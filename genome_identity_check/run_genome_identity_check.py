#!/usr/bin/env python
"""Genome-wide identity check: all-vs-all ANI (skani) between the GTDB r232
reference genomes of the 7 M. kansasii complex species, i.e. the same
genomes main1/main2 use. It answers whether a genome whose marker genes
look like another species (e.g. the persicum-like hsp65 of
GCF_002705785.1/825.1/865.1, or the persicum-type 16S of GCF_002086895.1)
is really a different species genome-wide. See LIT.md.

Outputs (in RESULTS):
- genomes.txt: the .fna paths passed to skani.
- skani_ani.tsv: raw skani triangle -E output, one row per genome pair.
- skani_long.tsv: every pair in both directions (a, sa = accession and
  species of one genome; b, sb = the other), with ANI and AF, the smaller of
  the two genomes' aligned fractions (a conservative, symmetric value).
- species_ani_summary.tsv: per species pair, the ANI range over all
  genome pairs.
- skani_version.txt

Run via submit_genome_identity_check.sbatch (requires env_mlsa and
apptainer).
"""
from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

from mlsa import GTDB_MLSA_ROOT, SPECIES  # noqa: E402
from mlsa.align import run_skani_triangle, skani_version  # noqa: E402
from mlsa.loci import discover_genomes  # noqa: E402

RESULTS = Path("/shares/sander.imm.uzh/MM/kansasii/output") / "mlsa" / "genome_identity_check"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("genome_identity_check")


def long_table(raw: pd.DataFrame) -> pd.DataFrame:
    """One row per ordered genome pair (both directions). AF is the smaller
    of the two aligned fractions, so both directions get the same value."""
    def ids(col: str) -> tuple[pd.Series, pd.Series]:
        paths = raw[col].map(Path)
        return paths.map(lambda p: p.stem), paths.map(lambda p: p.parent.name)

    ref_acc, ref_sp = ids("Ref_file")
    qry_acc, qry_sp = ids("Query_file")
    af = raw[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
    fwd = pd.DataFrame({"a": ref_acc, "sa": ref_sp, "b": qry_acc, "sb": qry_sp,
                        "ANI": raw["ANI"], "AF": af})
    rev = pd.DataFrame({"a": qry_acc, "sa": qry_sp, "b": ref_acc, "sb": ref_sp,
                        "ANI": raw["ANI"], "AF": af})
    return pd.concat([fwd, rev], ignore_index=True)


def species_summary(long: pd.DataFrame) -> pd.DataFrame:
    """ANI range per species pair, each genome pair counted once, species in
    SPECIES order."""
    rank = {sp: i for i, sp in enumerate(SPECIES)}
    ra, rb = long["sa"].map(rank), long["sb"].map(rank)
    pairs = long[(ra < rb) | ((ra == rb) & (long["a"] < long["b"]))]
    summary = (pairs.groupby(["sa", "sb"])["ANI"]
               .agg(n_pairs="count", min_ani="min", max_ani="max")
               .reset_index()
               .rename(columns={"sa": "species_a", "sb": "species_b"}))
    return summary.sort_values(["species_a", "species_b"], key=lambda c: c.map(rank))


def main(results_dir: Path = RESULTS) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    genomes = discover_genomes(GTDB_MLSA_ROOT, SPECIES)
    genome_list = results_dir / "genomes.txt"
    genome_list.write_text("".join(f"{p}\n" for p in sorted(str(g.fna_path) for g in genomes)))
    log.info("Running skani on %d genomes", len(genomes))

    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
    raw_tsv = run_skani_triangle(genome_list, results_dir / "skani_ani.tsv", threads)
    (results_dir / "skani_version.txt").write_text(skani_version() + "\n")

    raw = pd.read_csv(raw_tsv, sep="\t")
    n_expected = len(genomes) * (len(genomes) - 1) // 2
    if len(raw) != n_expected:
        log.warning("skani reported %d of %d genome pairs (pairs below the screening cutoff are omitted)",
                    len(raw), n_expected)
    long = long_table(raw)
    long.to_csv(results_dir / "skani_long.tsv", sep="\t", index=False)
    species_summary(long).to_csv(results_dir / "species_ani_summary.tsv", sep="\t", index=False,
                                 float_format="%.2f")
    log.info("Wrote %s", results_dir)


if __name__ == "__main__":
    main()
