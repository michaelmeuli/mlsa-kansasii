#!/usr/bin/env python
"""Write one FASTA per reference genome with its in-silico hsp65 amplicon (main1's
alignments/hsp65.raw.fasta), the reference set for hsp65 Sanger reads in main3 and
sanger-microsynth. File name <accession>_hsp65.fasta, header ">acc Mycobacterium <species> strain <acc>"
so mlsa.refalign.load_references reads the species from it.

Usage: python scripts/make_hsp65_references.py [OUT_DIR]
"""
from __future__ import annotations

import sys
from pathlib import Path

from Bio import SeqIO

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii")
SRC = ROOT / "output" / "mlsa" / "main1_locus_discovery" / "alignments" / "hsp65.raw.fasta"
DEFAULT_OUT = ROOT / "data" / "gtdb_genomes" / "Mycobacteriaceae" / "kansasii_complex_hsp65_amplicons"


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for rec in SeqIO.parse(SRC, "fasta"):
        species, acc = rec.id.split("__", 1)
        (out / f"{acc}_hsp65.fasta").write_text(
            f">{acc} Mycobacterium {species} strain {acc}\n{str(rec.seq).upper()}\n", encoding="utf-8")
        n += 1
    print(f"wrote {n} amplicons to {out}")


if __name__ == "__main__":
    main()
