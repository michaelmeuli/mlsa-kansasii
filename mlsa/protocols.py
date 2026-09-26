"""Parse PCR primer sequences out of the lab's Sanger sequencing protocol
.docx files, so genome-derived loci can be extracted from the exact same
amplicon the wet-lab assay targets.

Only hsp65 needs primer-based in-silico PCR here: the two GroEL paralogs in
Mycobacterium (groEL1/groEL2) are annotated identically in NCBI GFFs
(gene=groL, product="chaperonin GroEL"), so gene-symbol extraction cannot
tell them apart, while the TB-11/TB-12(w) primers only bind the groEL2 copy
the hsp65 assay actually targets. 16S is extracted directly from the GFF
rRNA/16S ribosomal RNA feature instead (see mlsa.loci.extract_16s, which
picks the genome's own copy where contaminant contigs add extra ones), and
using the full-length gene
is a strict superset of whatever fragment any particular 16S sequencing
primer pair amplifies, so it aligns correctly against the lab's 16S Sanger
reads regardless of exactly which primers were used for a given read.
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path

HSP65_PROTOCOL = Path(
    "/shares/sander.imm.uzh/MM/kansasii/data/protocols/"
    "32.43505-4-Mykobakterien NTM 65kDa.docx"
)

# Confirmed by reading 32.43505-4 (Telenti et al. 1993, J Clin Microbiol
# 31:175-178): routine primer pair is TB-11 / TB-12w (fallback TB-12 if
# TB-12w gives no amplicon), 439 bp amplicon. Kept here as a fallback in case
# docx text extraction ever fails (e.g. the file is replaced/reformatted).
_FALLBACK_HSP65_FWD = "ACCAACGATGGTGTGTCCAT"          # TB-11
_FALLBACK_HSP65_REV = [
    "CTTGTCGAASCGCATACCCT",                            # TB-12w (S = C/G)
    "CTTGTCGAACCGCATACCCT",                             # TB-12
]

_PRIMER_SEQ_RE = r"([ACGTRYSWKMBDHVN ]+?)-?3[’']"


def _docx_text(path: Path) -> str:
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8", "ignore")
    text = re.sub(r"<[^>]+>", "", xml)
    return re.sub(r"\s+", " ", text)


def get_hsp65_primers(protocol_path: Path = HSP65_PROTOCOL) -> tuple[str, list[str]]:
    """Return (forward_primer, [reverse_primer_options]) for the hsp65 assay,
    parsed live from the protocol .docx (TB-11 forward; TB-12w tried before
    TB-12 as reverse, matching the lab's routine-then-fallback order). Falls
    back to hard-coded values confirmed against the same document if parsing
    fails for any reason.
    """
    try:
        text = _docx_text(protocol_path)
        fwd_m = re.search(r"TB-11\D*?5[’']\s*" + _PRIMER_SEQ_RE, text)
        revw_m = re.search(r"TB-12w\D*?5[’']\s*" + _PRIMER_SEQ_RE, text)
        rev_m = re.search(r"TB-12(?!w)\D*?5[’']\s*" + _PRIMER_SEQ_RE, text)

        fwd = fwd_m.group(1).replace(" ", "") if fwd_m else _FALLBACK_HSP65_FWD
        revs = []
        if revw_m:
            revs.append(revw_m.group(1).replace(" ", ""))
        if rev_m:
            revs.append(rev_m.group(1).replace(" ", ""))
        if not revs:
            revs = list(_FALLBACK_HSP65_REV)
        return fwd, revs
    except Exception:
        return _FALLBACK_HSP65_FWD, list(_FALLBACK_HSP65_REV)
