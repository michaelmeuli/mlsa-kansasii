#!/usr/bin/env python
"""Main program 2 re-run against the reduced reference set of
main1_locus_discovery_excluded.

Same analysis as main2_sanger_differentiation/run_sanger_differentiation.py
(which this imports and calls). Reference sequences and the per-locus
pair-separation table are read from output/mlsa/main1_locus_discovery_excluded
(GCF_002705785.1, GCF_002705825.1 and GCF_002705865.1 removed), and results
are written to output/mlsa/main2_sanger_differentiation_excluded. Run after
main1_locus_discovery_excluded.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "main2_sanger_differentiation"))

import run_sanger_differentiation  # noqa: E402

EXCLUDED_ACCESSIONS = ("GCF_002705785.1", "GCF_002705825.1", "GCF_002705865.1")
OUTPUT = Path("/shares/sander.imm.uzh/MM/kansasii/output") / "mlsa"

if __name__ == "__main__":
    run_sanger_differentiation.main(
        results_dir=OUTPUT / "main2_sanger_differentiation_excluded",
        main1_results=OUTPUT / "main1_locus_discovery_excluded",
        exclude_accessions=EXCLUDED_ACCESSIONS,
    )
