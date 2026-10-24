#!/usr/bin/env python
"""Which GTDB-kansasii assemblies of the public BioProjects carry a persicum-like hsp65?

Same ceiling analysis as insilico_vs_wgs() in run_sanger_vs_wgs.py (hsp65 amplicon in silico from the
isolate's own assembly with the TB-11/TB-12(w) primers, classified by sanger_refalign's logic), but for the public
runs (runs/public_<BioProject>/, immensekansasii) and against the 72 per-genome hsp65 amplicons, so that the three
kansasii genomes with an atypical (persicum-like) hsp65 (mlsa.refalign.ATYPICAL_HSP65) are separate references.

Only samples that GTDB-Tk calls Mycobacterium kansasii are analysed. Per sample:
  closest_species / closest_accession / atypical  as in sanger_refalign (atypical = closest reference is an atypical-hsp65 kansasii genome)
  persicum_like = closest species is persicum or the closest reference is an atypical-hsp65 kansasii genome
Writes output/mlsa/sanger_vs_wgs_public/insilico_public.tsv and SUMMARY.txt.

Usage: run_insilico_public.py [--limit N] [--out DIR] [BioProject ...]   (default: the four public projects)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402
from Bio import SeqIO  # noqa: E402

from mlsa import DATA_ROOT, SPECIES  # noqa: E402
from mlsa.loci import in_silico_pcr  # noqa: E402
from mlsa.protocols import get_hsp65_primers  # noqa: E402
from mlsa.refalign import identify_read, load_references  # noqa: E402

ROOT = DATA_ROOT.parent
PROJECTS = ["PRJNA1078897", "PRJNA1048499", "PRJNA845115", "PRJNA323639"]
AMPLICONS = DATA_ROOT / "gtdb_genomes" / "Mycobacteriaceae" / "kansasii_complex_hsp65_amplicons"
GTDB_KANSASII = "s__Mycobacterium kansasii"


def gtdb_species(project: str) -> dict[str, str]:
    """sample -> GTDB-Tk species of every sample of the project that has a call."""
    out: dict[str, str] = {}
    results = ROOT / "runs" / f"public_{project}" / "assembly" / "results"
    for f in sorted(results.glob("*/3_quality/GTDB/*_gtdb_summary.tsv")):
        t = pd.read_csv(f, sep="\t", dtype=str)
        for _, r in t.iterrows():
            out[r["user_genome"]] = r["classification"].split(";")[-1]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("projects", nargs="*", default=PROJECTS)
    ap.add_argument("--limit", type=int, help="only the first N kansasii samples per project (test)")
    ap.add_argument("--species", default=GTDB_KANSASII, help="GTDB species to analyse (default kansasii; persicum as positive control)")
    ap.add_argument("--out", type=Path, default=ROOT / "output" / "mlsa" / "sanger_vs_wgs_public")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    refs = load_references(AMPLICONS, SPECIES)
    fwd, rev_options = get_hsp65_primers()
    rows = []
    for proj in args.projects:
        species = gtdb_species(proj)
        todo = [s for s, sp in species.items() if sp == args.species]
        print(f"{proj}: {len(species)} with GTDB call, {len(todo)} {args.species}", flush=True)
        for sid in todo[:args.limit]:
            f = ROOT / "runs" / f"public_{proj}" / "assembly" / "results" / sid / "1_unicycler" / f"{sid}.fasta"
            row: dict[str, object] = {"sample": sid, "bioproject": proj, "gtdb_species": species[sid], "amplicon_bp": 0}
            if f.exists():
                amp = in_silico_pcr({r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}, fwd, rev_options)
                if amp:
                    row["amplicon_bp"] = len(amp)
                    res = identify_read(refs, sid, "hsp65", amp)
                    b, s = res.best, res.runner_up
                    row.update(status=res.status, closest_species=res.closest_species, atypical=res.atypical)
                    if b is not None:
                        row.update(closest_accession=res.hits[b].ref.accession, diffs=res.diffs[b],
                                   identity=round(100 * res.identity(b), 2))
                    if s is not None:
                        row.update(second_species=res.hits[s].ref.species, margin=res.diffs[s] - res.diffs[b])
            rows.append(row)
            print(row, flush=True)

    df = pd.DataFrame(rows)
    for c in ("status", "closest_species", "atypical", "second_species"):
        if c not in df:
            df[c] = None
    df["persicum_like"] = (df["closest_species"] == "persicum") | (df["atypical"] == True)  # noqa: E712
    df.to_csv(args.out / "insilico_public.tsv", sep="\t", index=False)

    found = df[df["amplicon_bp"] > 0]
    lines = [f"GTDB-kansasii samples analysed: {len(df)}; hsp65 amplicon found in {len(found)}, no amplicon in {len(df) - len(found)}", ""]
    lines.append("closest hsp65 reference (all samples with amplicon):")
    for label, mask in (("kansasii (typical)", (found["closest_species"] == "kansasii") & (found["atypical"] != True)),  # noqa: E712
                        ("kansasii, atypical (persicum-like) hsp65 genome", found["atypical"] == True),  # noqa: E712
                        ("persicum", found["closest_species"] == "persicum")):
        lines.append(f"  {label}: {mask.sum()}")
    other = found[~found["closest_species"].isin(["kansasii", "persicum"]) & (found["atypical"] != True)]  # noqa: E712
    lines.append(f"  other species: {len(other)} " + str(other["closest_species"].value_counts().to_dict()))
    pl = found[found["persicum_like"]]
    lines += ["", f"persicum-like hsp65 but GTDB kansasii: {len(pl)} of {len(found)}",
              "  by status: " + str(pl["status"].value_counts().to_dict()),
              "  by project: " + str(pl["bioproject"].value_counts().to_dict()),
              f"  kansasii closest but persicum runner-up within 2 diffs (ambiguous): "
              f"{((found['closest_species'] == 'kansasii') & (found['second_species'] == 'persicum') & (found['status'] == 'ambiguous')).sum()}"]
    (args.out / "SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
