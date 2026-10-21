#!/usr/bin/env python
"""Write one IUPAC consensus hsp65 reference per species (main1's alignments/hsp65.aligned.fasta): a column
with more than one base among the species' genomes becomes the IUPAC code of those bases. Experiment for
`run_reference_alignment.py --iupac`; the genome-level references have no ambiguity codes.

Usage: python scripts/05_make_iupac_hsp65_references.py MIN_MINOR OUT_DIR
  MIN_MINOR  a base counts only if at least this many genomes of the species carry it (1 = every polymorphic
             site, 2 = ignore singleton variants)
Header ">REFSP_<species> Mycobacterium <species> strain consensus" so mlsa.refalign.load_references reads the species.
"""
from __future__ import annotations

import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mlsa.refalign import BASES_IUPAC  # noqa: E402

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii"))
SRC = ROOT / "output" / "mlsa" / "main1_locus_discovery" / "alignments" / "hsp65.aligned.fasta"


def consensus(seqs: list[str], min_minor: int) -> str:
    out = []
    for col in zip(*seqs):
        cnt = Counter(c for c in col if c in "ACGT")
        if not cnt:
            out.append("-")
            continue
        keep = {b for b, n in cnt.items() if n >= min_minor} or {cnt.most_common(1)[0][0]}
        out.append(BASES_IUPAC[frozenset(keep)])
    return "".join(out)


def main() -> None:
    min_minor, out = int(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    by_sp: dict[str, list[str]] = defaultdict(list)
    for rec in SeqIO.parse(SRC, "fasta"):
        by_sp[rec.id.split("__", 1)[0]].append(str(rec.seq).upper())
    for sp, seqs in by_sp.items():
        cons = consensus(seqs, min_minor).replace("-", "")
        amb = sum(c not in "ACGT" for c in cons)
        (out / f"REFSP_{sp}_hsp65.fasta").write_text(
            f">REFSP_{sp} Mycobacterium {sp} strain consensus\n{cons}\n", encoding="utf-8")
        print(f"{sp}: {len(seqs)} genomes, {len(cons)} bp, {amb} ambiguous sites")


if __name__ == "__main__":
    main()
