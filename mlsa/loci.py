"""Candidate locus definitions and extraction from the GTDB reference
genomes symlinked under data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii/.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

from Bio import SeqIO
from Bio.Seq import Seq

from .protocols import get_hsp65_primers

# Genes extracted unambiguously by their GFF gene= symbol (single copy,
# confirmed on the GCF_000157895.3 type-strain annotation as gene=secA and
# gene=tuf). A couple of synonyms are included per locus since the 72
# reference genomes come from different NCBI annotation pipeline runs that
# don't always agree on gene naming (e.g. secA vs secA1, tuf vs tufA).
GENE_SYMBOL_LOCI: dict[str, list[str]] = {
    "rpoB": ["rpoB"],
    "gyrA": ["gyrA"],
    "gyrB": ["gyrB"],
    "recA": ["recA"],
    "secA1": ["secA", "secA1"],
    "tuf": ["tuf", "tufA"],
}


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


@dataclass
class GenomeRecord:
    accession: str
    species: str
    fna_path: Path
    gff_path: Path
    _contigs_cache: dict[str, str] | None = field(default=None, repr=False, compare=False)
    _features_cache: list[tuple] | None = field(default=None, repr=False, compare=False)

    def contigs(self) -> dict[str, str]:
        # Every candidate locus is extracted separately for every genome, so
        # this (and gff_features below) get called many times per genome;
        # cache on first parse rather than re-reading the .fna/.gff each time.
        if self._contigs_cache is None:
            self._contigs_cache = {rec.id: str(rec.seq) for rec in SeqIO.parse(self.fna_path, "fasta")}
        return self._contigs_cache

    def gff_features(self) -> Iterator[tuple[str, str, str, int, int, str, dict]]:
        """Yield (seqid, source, ftype, start, end, strand, attrs) for every
        9-column feature line in the GFF."""
        if self._features_cache is None:
            features = []
            with open(self.gff_path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 9:
                        continue
                    seqid, source, ftype, start, end, score, strand, phase, attrs_str = fields
                    features.append((seqid, source, ftype, int(start), int(end), strand, _parse_attrs(attrs_str)))
            self._features_cache = features
        yield from self._features_cache


def _parse_attrs(attr_field: str) -> dict[str, str]:
    attrs = {}
    for part in attr_field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def discover_genomes(mlsa_root: Path, species_list: Iterable[str]) -> list[GenomeRecord]:
    records = []
    for species in species_list:
        sp_dir = mlsa_root / species
        if not sp_dir.is_dir():
            continue
        for fna in sorted(sp_dir.glob("*.fna")):
            accession = fna.stem
            gff = sp_dir / f"{accession}.gff"
            if gff.exists():
                records.append(GenomeRecord(accession, species, fna, gff))
    return records


def extract_gene_by_symbol(genome: GenomeRecord, gene_symbols: list[str]) -> str | None:
    """Extract the nucleotide sequence of a single-copy gene by exact GFF
    gene= symbol match (case-insensitive). If more than one feature matches
    (shouldn't happen for the symbols in GENE_SYMBOL_LOCI, but genomes are
    inconsistently annotated), the longest is kept."""
    wanted = {g.lower() for g in gene_symbols}
    hits = []
    for seqid, source, ftype, start, end, strand, attrs in genome.gff_features():
        if ftype != "gene":
            continue
        if attrs.get("gene", "").lower() in wanted:
            hits.append((seqid, start, end, strand))
    if not hits:
        return None
    hits.sort(key=lambda h: h[2] - h[1], reverse=True)
    seqid, start, end, strand = hits[0]
    contig_seq = genome.contigs().get(seqid)
    if contig_seq is None:
        return None
    frag = contig_seq[start - 1:end]
    return revcomp(frag) if strand == "-" else frag


def extract_16s(genome: GenomeRecord) -> str | None:
    """Extract the full-length 16S rRNA gene via its rRNA/16S ribosomal RNA
    product annotation (single copy, unambiguous)."""
    for seqid, source, ftype, start, end, strand, attrs in genome.gff_features():
        if ftype != "rRNA":
            continue
        if "16S ribosomal RNA" in attrs.get("product", ""):
            contig_seq = genome.contigs().get(seqid)
            if contig_seq is None:
                continue
            frag = contig_seq[start - 1:end]
            return revcomp(frag) if strand == "-" else frag
    return None


def extract_its(genome: GenomeRecord) -> str | None:
    """Extract the 16S-23S intergenic spacer (ITS): the gap between the 16S
    and 23S rRNA features on the same contig."""
    by_contig: dict[str, list[tuple[int, int, str, str]]] = {}
    for seqid, source, ftype, start, end, strand, attrs in genome.gff_features():
        if ftype != "rRNA":
            continue
        product = attrs.get("product", "")
        if "16S ribosomal RNA" in product or "23S ribosomal RNA" in product:
            by_contig.setdefault(seqid, []).append((start, end, strand, product))

    contigs = genome.contigs()
    for seqid, feats in by_contig.items():
        s16 = next((f for f in feats if "16S" in f[3]), None)
        s23 = next((f for f in feats if "23S" in f[3]), None)
        if not s16 or not s23:
            continue
        strand = s16[2]
        if strand == "+":
            gap_start, gap_end = s16[1], s23[0] - 1
        else:
            gap_start, gap_end = s23[1], s16[0] - 1
        if gap_end <= gap_start:
            continue
        contig_seq = contigs.get(seqid)
        if contig_seq is None:
            continue
        frag = contig_seq[gap_start:gap_end]
        return revcomp(frag) if strand == "-" else frag
    return None


_BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
_IUPAC_EXPAND = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _seq_to_codes(seq: str) -> "np.ndarray":
    import numpy as np
    arr = np.frombuffer(seq.upper().encode("ascii", errors="replace"), dtype=np.uint8)
    codes = np.full(arr.shape, -1, dtype=np.int8)
    for base, code in _BASE_CODE.items():
        codes[arr == ord(base)] = code
    return codes


def find_primer_sites(seq: str, primer: str, max_mismatches: int) -> list[int]:
    """Sliding-window search (vectorized with numpy) for approximate matches
    of `primer` (IUPAC degenerate bases honoured) in `seq`, allowing up to
    max_mismatches substitutions. Real PCR primers routinely tolerate a
    handful of mismatches (especially away from the 3' end), so an
    exact-match search misses most true binding sites in genomes that are
    not identical to whatever strain the primer was originally designed
    against. Returns the 0-based start positions of every matching window."""
    import numpy as np
    codes = _seq_to_codes(seq)
    L = len(primer)
    n = len(codes) - L + 1
    if n <= 0:
        return []
    allowed = np.zeros((L, 4), dtype=bool)
    for i, base in enumerate(primer.upper()):
        for allowed_base in _IUPAC_EXPAND.get(base, "ACGT"):
            allowed[i, _BASE_CODE[allowed_base]] = True

    windows = np.lib.stride_tricks.sliding_window_view(codes, L)
    # Per-base broadcast instead of fancy-indexing every (window, position)
    # pair: the latter materializes O(n * L) index arrays and blows up
    # memory on multi-Mb genomes. Ambiguous/N bases in the genome (code -1)
    # never equal 0..3, so they fall out as mismatches automatically.
    matches = np.zeros(windows.shape, dtype=bool)
    for base_code in range(4):
        matches |= (windows == base_code) & allowed[:, base_code][None, :]
    mismatches = L - matches.sum(axis=1)
    return np.where(mismatches <= max_mismatches)[0].tolist()


def in_silico_pcr(contigs: dict[str, str], fwd_primer: str, rev_primer_options: list[str],
                   max_amplicon: int = 2000, max_mismatches: int = 3) -> str | None:
    """Find a primer-pair binding site in any contig/orientation and return
    the amplicon (inclusive of both primers), allowing up to max_mismatches
    per primer (see find_primer_sites). Tries each reverse-primer option in
    order (e.g. TB-12w before the TB-12 fallback) and returns the first
    valid (fwd upstream of rev, within max_amplicon) pairing found."""
    for rev_primer in rev_primer_options:
        rev_rc = revcomp(rev_primer)
        for contig_seq in contigs.values():
            for strand_seq in (contig_seq, revcomp(contig_seq)):
                fwd_sites = find_primer_sites(strand_seq, fwd_primer, max_mismatches)
                if not fwd_sites:
                    continue
                rev_sites = find_primer_sites(strand_seq, rev_rc, max_mismatches)
                if not rev_sites:
                    continue
                for fpos in sorted(fwd_sites):
                    for rpos in sorted(rev_sites):
                        rend = rpos + len(rev_rc)
                        if fpos < rend <= fpos + max_amplicon:
                            return strand_seq[fpos:rend]
    return None


def extract_hsp65(genome: GenomeRecord) -> str | None:
    fwd, rev_options = get_hsp65_primers()
    return in_silico_pcr(genome.contigs(), fwd, rev_options)


CANDIDATE_LOCI = list(GENE_SYMBOL_LOCI.keys()) + ["16S", "ITS", "hsp65"]


def extract_all_loci(genome: GenomeRecord) -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    for locus, symbols in GENE_SYMBOL_LOCI.items():
        result[locus] = extract_gene_by_symbol(genome, symbols)
    result["16S"] = extract_16s(genome)
    result["ITS"] = extract_its(genome)
    result["hsp65"] = extract_hsp65(genome)
    return result
