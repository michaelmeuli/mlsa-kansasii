#!/usr/bin/env python
"""Write one FASTA per reference genome with its in-silico hsp65 amplicon (locus_discovery's
alignments/hsp65.raw.fasta), the reference set (HSP65_AMPLICONS) for hsp65 Sanger reads in 09_sanger_refalign,
sanger-microsynth, hsp65-refalign and the public hsp65 check. File name <accession>_hsp65.fasta, header ">acc Mycobacterium <species> strain <acc>"
so mlsa.refalign.load_references reads the species from it.

Usage: python scripts/06_make_hsp65_references.py [--out HSP65_AMPLICONS_DIR] [--src hsp65.raw.fasta]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the mlsa package
from mlsa.kansasii_paths import HSP65_AMPLICONS, OUTPUT_ROOT  # noqa: E402

SRC = OUTPUT_ROOT / "mlsa" / "locus_discovery" / "alignments" / "hsp65.raw.fasta"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=HSP65_AMPLICONS)
    ap.add_argument("--src", type=Path, default=SRC)
    args = ap.parse_args()
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for rec in SeqIO.parse(args.src, "fasta"):
        species, acc = rec.id.split("__", 1)
        (out / f"{acc}_hsp65.fasta").write_text(
            f">{acc} Mycobacterium {species} strain {acc}\n{str(rec.seq).upper()}\n", encoding="utf-8")
        n += 1
    print(f"wrote {n} amplicons to {out}")


if __name__ == "__main__":
    main()
