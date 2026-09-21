#!/usr/bin/env python
"""Main program 1 re-run with three reference genomes excluded.

Same analysis as main1_locus_discovery/run_locus_discovery.py (which this
imports and calls), but GCF_002705785.1, GCF_002705825.1 and
GCF_002705865.1 are dropped before locus extraction, and everything is
written to output/mlsa/main1_locus_discovery_excluded instead.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "main1_locus_discovery"))

import run_locus_discovery  # noqa: E402

EXCLUDED_ACCESSIONS = ("GCF_002705785.1", "GCF_002705825.1", "GCF_002705865.1")
RESULTS_EXCLUDED = (
    Path("/shares/sander.imm.uzh/MM/kansasii/output") / "mlsa" / "main1_locus_discovery_excluded"
)

if __name__ == "__main__":
    run_locus_discovery.main(results_dir=RESULTS_EXCLUDED, exclude_accessions=EXCLUDED_ACCESSIONS)
