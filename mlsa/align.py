"""Alignment, tree-building, and distance/discrimination helpers.

Heavy binaries (mafft, iqtree) run inside the shared Singularity containers
used by immensekansasii rather than being installed into env_mlsa, per the
project's "keep the login-node conda env light" requirement. blast/iqtree
are reused as-is from immensekansasii's container set; mafft is pulled by
Singularity/pull_singularity_img.sh. See that script for exact image names.
"""
from __future__ import annotations

import itertools
import subprocess
from pathlib import Path

import pandas as pd
from Bio import Align, SeqIO
from Bio.Seq import Seq

CONTAINER_ROOT = Path(
    "/shares/sander.imm.uzh/software/pipelines/IMMense/IMMense_dependencies/containers"
)
MAFFT_IMG = CONTAINER_ROOT / "quay.io-biocontainers-mafft-7.525--h031d066_1.img"
IQTREE_IMG = CONTAINER_ROOT / "quay.io-biocontainers-iqtree-3.1.3--h8471819_0.img"


def _singularity_exec(img: Path, cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    full_cmd = ["singularity", "exec", "--bind", "/shares", str(img)] + cmd
    return subprocess.run(full_cmd, check=True, **kwargs)


def write_fasta(sequences: dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n{seq}\n")


def read_fasta(path: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}


def run_mafft(input_fasta: Path, output_fasta: Path) -> Path:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as out:
        _singularity_exec(MAFFT_IMG, ["mafft", "--auto", "--quiet", str(input_fasta)], stdout=out)
    return output_fasta


def run_iqtree(alignment_fasta: Path, prefix: Path, model: str = "GTR+G") -> Path:
    """Build a maximum-likelihood tree with IQ-TREE. Returns the .treefile
    path. Overwrites any previous run at the same prefix (-redo)."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    _singularity_exec(
        IQTREE_IMG,
        ["iqtree", "-s", str(alignment_fasta), "-m", model, "-T", "AUTO",
         "--prefix", str(prefix), "--quiet", "--redo"],
    )
    # IQ-TREE appends ".treefile" to the prefix string as-is; Path.with_suffix
    # would instead replace any existing suffix in `prefix` (e.g. turning
    # ".../16S.tree" into ".../16S.treefile" instead of
    # ".../16S.tree.treefile"), so build the path by straight concatenation.
    return Path(f"{prefix}.treefile")


def orient_to_reference(seq: str, reference: str) -> str:
    """Return seq in whichever orientation (as-is or reverse-complemented)
    aligns better to reference, using a quick local alignment score. Used to
    fix the inconsistent historical fwd/rev naming in the lab's Sanger
    filenames without relying on fragile filename parsing."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    fwd_score = aligner.score(reference, seq)
    rc = str(Seq(seq).reverse_complement())
    rev_score = aligner.score(reference, rc)
    return rc if rev_score > fwd_score else seq


def p_distance_matrix(alignment: dict[str, str]) -> pd.DataFrame:
    """Pairwise uncorrected p-distance (fraction mismatched over aligned
    columns where neither sequence has a gap) for an already-aligned set of
    equal-length sequences."""
    names = list(alignment.keys())
    dist = pd.DataFrame(0.0, index=names, columns=names)
    for a, b in itertools.combinations(names, 2):
        sa, sb = alignment[a], alignment[b]
        compared = mismatches = 0
        for ca, cb in zip(sa, sb):
            if ca == "-" or cb == "-":
                continue
            compared += 1
            if ca.upper() != cb.upper():
                mismatches += 1
        d = mismatches / compared if compared else float("nan")
        dist.loc[a, b] = dist.loc[b, a] = d
    return dist


def barcoding_gap_check(dist: pd.DataFrame, species_of: dict[str, str]) -> pd.DataFrame:
    """For every pair of species present in dist, compare the minimum
    inter-species distance to the maximum intra-species distance of either
    species (the DNA-barcoding "gap" criterion). A pair is "separated" when
    the closest pair of sequences from the two different species is still
    farther apart than the most divergent pair of sequences within either
    species. Returns one row per species pair."""
    species_seqs: dict[str, list[str]] = {}
    for name, sp in species_of.items():
        if name in dist.index:
            species_seqs.setdefault(sp, []).append(name)

    def max_intra(sp: str) -> float:
        seqs = species_seqs.get(sp, [])
        if len(seqs) < 2:
            return 0.0
        return max(dist.loc[a, b] for a, b in itertools.combinations(seqs, 2))

    def min_inter(sp_a: str, sp_b: str) -> float:
        seqs_a, seqs_b = species_seqs.get(sp_a, []), species_seqs.get(sp_b, [])
        if not seqs_a or not seqs_b:
            return float("nan")
        return min(dist.loc[a, b] for a in seqs_a for b in seqs_b)

    rows = []
    species = sorted(species_seqs.keys())
    for sp_a, sp_b in itertools.combinations(species, 2):
        intra = max(max_intra(sp_a), max_intra(sp_b))
        inter = min_inter(sp_a, sp_b)
        rows.append({
            "species_a": sp_a,
            "species_b": sp_b,
            "n_a": len(species_seqs.get(sp_a, [])),
            "n_b": len(species_seqs.get(sp_b, [])),
            "max_intra": intra,
            "min_inter": inter,
            "gap_margin": inter - intra,
            "separated": bool(inter > intra),
        })
    return pd.DataFrame(rows)


def resolves_all_species(dist: pd.DataFrame, species_of: dict[str, str]) -> tuple[bool, pd.DataFrame]:
    gap_df = barcoding_gap_check(dist, species_of)
    return bool(gap_df["separated"].all()), gap_df
