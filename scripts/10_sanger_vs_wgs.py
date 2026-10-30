#!/usr/bin/env python
"""Is hsp65 alone enough? Benchmark hsp65 species calls against the
whole-genome (GTDB-Tk) species of the same isolates.

Two questions, both against the WGS species in output/immensekansasii/results/screening_map_results.csv
(complex species only):

  1. Sanger vs WGS: how often do the existing Sanger hsp65 calls (sanger_mlsa nearest species,
     sanger_mlsa_excluded, sanger_refalign plain alignment) match the genome species?
  2. Ceiling: extract the hsp65 amplicon in silico from each isolate's own assembly (TB-11/TB-12(w)
     primers) and classify it with sanger_refalign's logic. This is what a perfect hsp65 read would say, so
     it separates failed Sanger reads (assay) from the limits of the marker itself.

Needs the outputs of 08_sanger_mlsa (+ _excluded) and 09_sanger_refalign. Writes to
output/mlsa-kansasii/results/sanger_vs_wgs/: sanger_vs_wgs.tsv, insilico_vs_wgs.tsv, SUMMARY.txt.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

from mlsa.kansasii_paths import (  # noqa: E402
    ASSEMBLIES,
    GTDB_REPRESENTATIVES,
    OUT_MLSA,
    REFERENCE_STRAINS,
    SCREENING_MAP_RESULTS,
    SPECIES,
)
from mlsa.loci import in_silico_pcr  # noqa: E402
from mlsa.protocols import get_hsp65_primers  # noqa: E402
from mlsa.refalign import identify_read, load_references  # noqa: E402
from Bio import SeqIO  # noqa: E402

MLSA = OUT_MLSA / "results"
OUTPUT = MLSA / "sanger_vs_wgs"


def stem(read: str) -> str:
    """Read names in sanger_refalign are .ab1 file stems (they may contain dots, so no Path.stem)."""
    return read[:-4] if read.endswith(".ab1") else read


def sanger_vs_wgs(truth: pd.DataFrame, mlsa: Path) -> pd.DataFrame:
    reps = pd.read_csv(mlsa / "sanger_mlsa" / "representative_reads.tsv", sep="\t")
    reps = reps[reps["locus"] == "hsp65"].set_index("probennummer")
    ra = pd.read_csv(mlsa / "sanger_refalign" / "reference_alignment.tsv", sep="\t").drop_duplicates("read").set_index("read")
    df = truth[truth["PROBENNUMMER"].isin(reps.index)].copy()
    df["read"] = df["PROBENNUMMER"].map(reps["read"]).map(stem)
    for col in ("status", "closest_species", "closest_identity", "closest_diffs", "margin", "read_bp", "mixed_fraction"):
        df["sanger_refalign_" + col] = df["read"].map(ra[col])
    for tag, d in (("sanger_mlsa", "sanger_mlsa"), ("sanger_mlsa_excl", "sanger_mlsa_excluded")):
        c = pd.read_csv(mlsa / d / "isolate_classification.tsv", sep="\t")
        c = c[c["locus"] == "hsp65"].set_index("probennummer")
        df[tag + "_nearest"] = df["PROBENNUMMER"].map(c["nearest_species"])
        df[tag + "_unambiguous"] = df["PROBENNUMMER"].map(c["unambiguous"])
    return df


def insilico_vs_wgs(truth: pd.DataFrame, gtdb_representatives: Path, assemblies: Path) -> pd.DataFrame:
    refs = load_references(gtdb_representatives, SPECIES)
    fwd, rev_options = get_hsp65_primers()
    rows = []
    for _, t in truth.iterrows():
        f = assemblies / t["PROBENNUMMER"] / "1_unicycler" / f"{t['PROBENNUMMER']}.fasta"
        if not f.exists():
            continue
        contigs = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}
        amp = in_silico_pcr(contigs, fwd, rev_options)
        row: dict[str, object] = {"PROBENNUMMER": t["PROBENNUMMER"], "wgs_species": t["species"],
                                  "contaminated": t["contamination_flag"] != "", "amplicon_bp": len(amp) if amp else 0}
        if amp:
            res = identify_read(refs, t["PROBENNUMMER"], "hsp65", amp)
            b, s = res.best, res.runner_up
            row.update(status=res.status, closest_species=res.closest_species)
            if b is not None:
                row.update(identity=round(100 * res.identity(b), 2), diffs=res.diffs[b])
                if s is not None:
                    row.update(second_species=res.hits[s].ref.species, margin=res.diffs[s] - res.diffs[b])
        rows.append(row)
        print(row, flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtdb-representatives", type=Path, default=GTDB_REPRESENTATIVES)
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
    ap.add_argument("--assemblies", type=Path, default=ASSEMBLIES)
    ap.add_argument("--mlsa-dir", type=Path, default=MLSA, help="holds the sanger_mlsa*, sanger_refalign outputs")
    ap.add_argument("--out", type=Path, default=OUTPUT)
    args = ap.parse_args()
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(args.screening_map_results, dtype=str).fillna("")
    # Mkan329-126 to -133 are reference strains, not patient isolates
    truth = s[s["species"].isin(SPECIES) & ~s["PROBENNUMMER"].isin(REFERENCE_STRAINS)]
    sv = sanger_vs_wgs(truth, args.mlsa_dir)
    sv.to_csv(out / "sanger_vs_wgs.tsv", sep="\t", index=False)
    iv = insilico_vs_wgs(truth, args.gtdb_representatives, args.assemblies)
    iv.to_csv(out / "insilico_vs_wgs.tsv", sep="\t", index=False)

    lines = [f"WGS complex isolates: {len(truth)} ({(truth['contamination_flag'] != '').sum()} flagged contaminated)",
             f"with a Sanger hsp65 read: {len(sv)}", ""]
    for tag, call, called in (("sanger_mlsa (nearest species)", "sanger_mlsa_nearest", sv["sanger_mlsa_nearest"].notna()),
                              ("sanger_mlsa_excluded (nearest species)", "sanger_mlsa_excl_nearest", sv["sanger_mlsa_excl_nearest"].notna()),
                              ("sanger_refalign (status ok only)", "sanger_refalign_closest_species", sv["sanger_refalign_status"] == "ok")):
        d = sv[called]
        lines.append(f"{tag}: {len(d)} called, {(d[call] == d['species']).sum()} match WGS; {len(sv) - len(d)} no call")
    ok = iv[iv["status"] == "ok"] if "status" in iv else iv.iloc[0:0]
    lines += ["", f"in-silico hsp65 from own assembly: {len(iv)} isolates, amplicon found in {(iv['amplicon_bp'] > 0).sum()}",
              f"  sanger_refalign status ok: {len(ok)}, match WGS: {(ok['closest_species'] == ok['wgs_species']).sum()}"]
    if "status" in iv:
        lines.append("  status counts: " + ", ".join(f"{k} {v}" for k, v in iv["status"].value_counts().items()))
    (out / "SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
