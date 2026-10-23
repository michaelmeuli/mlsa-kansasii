#!/usr/bin/env python
"""Main program 2 re-run against the reduced reference set of
locus_discovery_excluded.

Same analysis as sanger_mlsa/run_sanger_mlsa.py
(which this imports and calls). Reference sequences and the per-locus
pair-separation table are read from output/mlsa/locus_discovery_excluded
(GCF_002705785.1, GCF_002705825.1 and GCF_002705865.1 removed), and results
are written to output/mlsa/sanger_mlsa_excluded. Run after
locus_discovery_excluded.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "sanger_mlsa"))

import run_sanger_mlsa  # noqa: E402

EXCLUDED_ACCESSIONS = ("GCF_002705785.1", "GCF_002705825.1", "GCF_002705865.1")
OUTPUT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii") + "/output") / "mlsa"

if __name__ == "__main__":
    run_sanger_mlsa.main(
        results_dir=OUTPUT / "sanger_mlsa_excluded",
        locus_discovery_results=OUTPUT / "locus_discovery_excluded",
        exclude_accessions=EXCLUDED_ACCESSIONS,
    )
