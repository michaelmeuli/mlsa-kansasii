#!/usr/bin/env python
"""Main program 4: is hsp65 alone enough? Benchmark hsp65 species calls against the
whole-genome (GTDB-Tk) species of the same isolates.

Two questions, both against the WGS species in output/screening_map_results.csv
(complex species only):

  1. Sanger vs WGS: how often do the existing Sanger hsp65 calls (main2 nearest species,
     main2_excluded, main3 plain alignment) match the genome species?
  2. Ceiling: extract the hsp65 amplicon in silico from each isolate's own assembly (TB-11/TB-12(w)
     primers) and classify it with main3's logic. This is what a perfect hsp65 read would say, so
     it separates failed Sanger reads (assay) from the limits of the marker itself.

Needs the outputs of main2 (+ _excluded) and main3. Writes to
output/mlsa/main4_hsp65_vs_wgs/: sanger_vs_wgs.tsv, insilico_vs_wgs.tsv, SUMMARY.txt.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

from mlsa import DATA_ROOT, GTDB_REPRESENTATIVES, SPECIES  # noqa: E402
from mlsa.loci import in_silico_pcr  # noqa: E402
from mlsa.protocols import get_hsp65_primers  # noqa: E402
from mlsa.refalign import identify_read, load_references  # noqa: E402
from Bio import SeqIO  # noqa: E402

OUT_ROOT = DATA_ROOT.parent / "output"
MLSA = OUT_ROOT / "mlsa"
OUTPUT = MLSA / "main4_hsp65_vs_wgs"
GENOMES = DATA_ROOT.parent / "runs" / "mkan329" / "mkan329_transfer_result" / "genomes"


def stem(read: str) -> str:
    """Read names in main3 are .ab1 file stems (they may contain dots, so no Path.stem)."""
    return read[:-4] if read.endswith(".ab1") else read


def sanger_vs_wgs(truth: pd.DataFrame) -> pd.DataFrame:
    reps = pd.read_csv(MLSA / "main2_sanger_differentiation" / "representative_reads.tsv", sep="\t")
    reps = reps[reps["locus"] == "hsp65"].set_index("probennummer")
    ra = pd.read_csv(MLSA / "main3_reference_alignment" / "reference_alignment.tsv", sep="\t").drop_duplicates("read").set_index("read")
    df = truth[truth["PROBENNUMMER"].isin(reps.index)].copy()
    df["read"] = df["PROBENNUMMER"].map(reps["read"]).map(stem)
    for col in ("status", "closest_species", "closest_identity", "closest_diffs", "margin", "read_bp", "mixed_fraction"):
        df["main3_" + col] = df["read"].map(ra[col])
    for tag, d in (("main2", "main2_sanger_differentiation"), ("main2_excl", "main2_sanger_differentiation_excluded")):
        c = pd.read_csv(MLSA / d / "isolate_classification.tsv", sep="\t")
        c = c[c["locus"] == "hsp65"].set_index("probennummer")
        df[tag + "_nearest"] = df["PROBENNUMMER"].map(c["nearest_species"])
        df[tag + "_unambiguous"] = df["PROBENNUMMER"].map(c["unambiguous"])
    return df


def insilico_vs_wgs(truth: pd.DataFrame) -> pd.DataFrame:
    refs = load_references(GTDB_REPRESENTATIVES, SPECIES)
    fwd, rev_options = get_hsp65_primers()
    rows = []
    for _, t in truth.iterrows():
        f = next(iter(GENOMES.glob(f"*/{t['PROBENNUMMER']}_*.fna")), None)
        if f is None:
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
    OUTPUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(OUT_ROOT / "screening_map_results.csv", dtype=str).fillna("")
    truth = s[s["species"].isin(SPECIES)]
    sv = sanger_vs_wgs(truth)
    sv.to_csv(OUTPUT / "sanger_vs_wgs.tsv", sep="\t", index=False)
    iv = insilico_vs_wgs(truth)
    iv.to_csv(OUTPUT / "insilico_vs_wgs.tsv", sep="\t", index=False)

    lines = [f"WGS complex isolates: {len(truth)} ({(truth['contamination_flag'] != '').sum()} flagged contaminated)",
             f"with a Sanger hsp65 read: {len(sv)}", ""]
    for tag, call, called in (("main2 (nearest species)", "main2_nearest", sv["main2_nearest"].notna()),
                              ("main2_excluded (nearest species)", "main2_excl_nearest", sv["main2_excl_nearest"].notna()),
                              ("main3 (status ok only)", "main3_closest_species", sv["main3_status"] == "ok")):
        d = sv[called]
        lines.append(f"{tag}: {len(d)} called, {(d[call] == d['species']).sum()} match WGS; {len(sv) - len(d)} no call")
    ok = iv[iv["status"] == "ok"] if "status" in iv else iv.iloc[0:0]
    lines += ["", f"in-silico hsp65 from own assembly: {len(iv)} isolates, amplicon found in {(iv['amplicon_bp'] > 0).sum()}",
              f"  main3 status ok: {len(ok)}, match WGS: {(ok['closest_species'] == ok['wgs_species']).sum()}"]
    if "status" in iv:
        lines.append("  status counts: " + ", ".join(f"{k} {v}" for k, v in iv["status"].value_counts().items()))
    (OUTPUT / "SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
