#!/usr/bin/env python
"""Main program 2: using the Sanger hsp65 + 16S reads the lab already has for
the TNR isolates in data/imm/screening_map_link.csv (symlinked into
data/sanger/seq_kansasii by scripts/link_sanger_kansasii.sh), assess how far
species-level differentiation actually gets today, and which additional
locus/loci (from main1_locus_discovery's candidate set) would need to be
sequenced to resolve the isolates that remain ambiguous.

Isolates are keyed by PROBENNUMMER: reads under any of an isolate's TNRs
(TNR, TNR_NGS, TNR3..TNR6) are pooled. For each isolate + locus, all
matching .ab1 reads are quality-trimmed and the
best (longest, then highest quality) read is taken as that isolate's
representative sequence (no fwd/rev consensus assembly in this first pass).
Representative sequences are oriented against a reference sequence, aligned
together with the 7-species reference sequences from main1 (reused from
output/mlsa/main1_locus_discovery/alignments/*.raw.fasta when available, else
re-extracted directly), and each isolate is classified as unambiguously
nested in one species' cluster, or ambiguous, using the same DNA-barcoding
gap logic as main1: the gap between the isolate's nearest and second-nearest
species must exceed the larger within-species diversity of the two among
the references.

Run via submit_sanger_differentiation.sbatch (requires env_mlsa active and
the mafft Singularity container pulled). Best run after main1, but falls
back to extracting reference sequences itself if main1's output is missing.
"""
from __future__ import annotations

import csv
import logging
import sys
from collections import defaultdict
from typing import Iterable
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

from mlsa import (  # noqa: E402
    GTDB_MLSA_ROOT,
    SANGER_LINK_DIR,
    SCREENING_MAP,
    SCREENING_MAP_TNR_COLUMNS,
    SPECIES,
)
from mlsa.align import (  # noqa: E402
    orient_to_reference,
    p_distance_matrix,
    read_fasta,
    run_mafft,
    write_fasta,
)
from mlsa.loci import discover_genomes, extract_16s, extract_hsp65  # noqa: E402
from mlsa.plotting import plot_alignment_heatmap, plot_tree  # noqa: E402
from mlsa.sanger_io import classify_locus, extract_tnr, load_trimmed_ab1  # noqa: E402

RESULTS = Path("/shares/sander.imm.uzh/MM/kansasii/output") / "mlsa" / "main2_sanger_differentiation"
MAIN1_RESULTS = Path("/shares/sander.imm.uzh/MM/kansasii/output") / "mlsa" / "main1_locus_discovery"
LOCI = ["hsp65", "16S"]

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("sanger_differentiation")


def load_tnr_to_probennummer() -> dict[str, str]:
    """Map every TNR in any of SCREENING_MAP_TNR_COLUMNS to its row's
    PROBENNUMMER. Reference strains have no TNR and are skipped."""
    tnr_to_pnr: dict[str, str] = {}
    with open(SCREENING_MAP) as fh:
        for row in csv.DictReader(fh):
            for col in SCREENING_MAP_TNR_COLUMNS:
                tnr = row[col]
                if not tnr:
                    continue
                if tnr_to_pnr.get(tnr, row["PROBENNUMMER"]) != row["PROBENNUMMER"]:
                    raise ValueError(f"TNR {tnr} listed for both {tnr_to_pnr[tnr]} "
                                     f"and {row['PROBENNUMMER']} in {SCREENING_MAP}")
                tnr_to_pnr[tnr] = row["PROBENNUMMER"]
    return tnr_to_pnr


def collect_isolate_reads(tnr_to_pnr: dict[str, str]) -> dict[str, dict[str, list[Path]]]:
    """locus -> {probennummer: [ab1 paths]}"""
    by_locus_pnr: dict[str, dict[str, list[Path]]] = {locus: defaultdict(list) for locus in LOCI}
    for path in sorted(SANGER_LINK_DIR.glob("*.ab1")):
        locus = classify_locus(path.name)
        if locus is None:
            continue
        tnr = extract_tnr(path.name, tnr_to_pnr)
        if tnr is None:
            continue
        by_locus_pnr[locus][tnr_to_pnr[tnr]].append(path)
    return by_locus_pnr


def pick_representative(paths: list[Path]) -> tuple[str, float, Path] | None:
    """Quality-trim every candidate read for a (probennummer, locus) and keep the
    longest (ties broken by mean quality) as the representative sequence."""
    candidates = []
    for path in paths:
        trimmed = load_trimmed_ab1(path)
        if trimmed is None:
            continue
        seq, mean_q = trimmed
        candidates.append((seq, mean_q, path))
    if not candidates:
        return None
    candidates.sort(key=lambda c: (len(c[0]), c[1]), reverse=True)
    return candidates[0]


def load_reference_sequences(locus: str, main1_results: Path = MAIN1_RESULTS,
                              exclude_accessions: frozenset[str] = frozenset()) -> dict[str, str]:
    """Reuse main1's already-extracted reference sequences when available;
    otherwise extract them directly so main2 can run standalone."""
    raw_fasta = main1_results / "alignments" / f"{locus}.raw.fasta"
    if raw_fasta.exists():
        log.info("Reusing main1 reference sequences for %s from %s", locus, raw_fasta)
        return read_fasta(raw_fasta)

    log.warning("main1 output for %s not found at %s; extracting references directly", locus, raw_fasta)
    genomes = [g for g in discover_genomes(GTDB_MLSA_ROOT, SPECIES)
               if g.accession not in exclude_accessions]
    extractor = extract_hsp65 if locus == "hsp65" else extract_16s
    refs = {}
    for g in genomes:
        seq = extractor(g)
        if seq:
            refs[f"{g.species}__{g.accession}"] = seq
    return refs


def species_of_name(name: str) -> str:
    prefix = name.split("__", 1)[0]
    return prefix if prefix in SPECIES else "isolate"


def max_intra_species(dist: pd.DataFrame, species_of: dict[str, str]) -> dict[str, float]:
    by_species: dict[str, list[str]] = defaultdict(list)
    for name, sp in species_of.items():
        if sp in SPECIES and name in dist.index:
            by_species[sp].append(name)
    result = {}
    for sp, names in by_species.items():
        if len(names) < 2:
            result[sp] = 0.0
            continue
        result[sp] = max(dist.loc[a, b] for i, a in enumerate(names) for b in names[i + 1:])
    return result


def classify_isolate(isolate_name: str, dist: pd.DataFrame, species_of: dict[str, str],
                      intra_max: dict[str, float]) -> dict | None:
    """Nearest and second-nearest species by the isolate's distance to each
    species' closest reference. The tolerance is the larger within-species
    diversity of the two species, as in main1's barcoding-gap check
    (mlsa.align.barcoding_gap_check). References whose distance is NaN (no
    overlapping, gap-free columns with the isolate's read) are skipped."""
    per_species_min: dict[str, float] = {}
    for ref_name, sp in species_of.items():
        if sp not in SPECIES or ref_name not in dist.columns:
            continue
        d = dist.loc[isolate_name, ref_name]
        if pd.isna(d):
            continue
        per_species_min[sp] = min(per_species_min.get(sp, float("inf")), d)
    if not per_species_min:
        log.warning("%s shares no aligned columns with any reference; not classified", isolate_name)
        return None
    ranked = sorted(per_species_min.items(), key=lambda kv: kv[1])
    best_sp, best_d = ranked[0]
    second_sp, second_d = ranked[1] if len(ranked) > 1 else (None, float("inf"))
    tolerance = max(intra_max.get(best_sp, 0.0), intra_max.get(second_sp, 0.0))
    margin = second_d - best_d
    return {
        "nearest_species": best_sp,
        "nearest_dist": best_d,
        "second_species": second_sp,
        "second_dist": second_d,
        "margin": margin,
        "tolerance": tolerance,
        "unambiguous": bool(margin > tolerance),
    }


def recommend_loci(ambiguous_df: pd.DataFrame, pair_sep_path: Path) -> pd.DataFrame:
    if not pair_sep_path.exists():
        log.warning("main1 pair-separation table not found at %s; skipping locus recommendations", pair_sep_path)
        return pd.DataFrame()
    pair_sep = pd.read_csv(pair_sep_path, sep="\t")
    pair_sep = pair_sep[~pair_sep["locus"].isin(LOCI)]

    rows = []
    for _, row in ambiguous_df.iterrows():
        if pd.isna(row["second_species"]):
            continue
        sp_a, sp_b = sorted([row["nearest_species"], row["second_species"]])
        candidates = pair_sep[
            ((pair_sep["species_a"] == sp_a) & (pair_sep["species_b"] == sp_b))
            & pair_sep["separated"]
        ].sort_values("gap_margin", ascending=False)
        recommended = candidates.iloc[0]["locus"] if len(candidates) else "none of the candidate loci"
        rows.append({"probennummer": row["probennummer"], "species_a": sp_a, "species_b": sp_b, "recommended_locus": recommended})
    return pd.DataFrame(rows)


def run_locus_or_combo(label: str, combined_seqs: dict[str, str], species_of: dict[str, str],
                        isolate_names: list[str], results_dir: Path = RESULTS) -> pd.DataFrame:
    aln_fasta = results_dir / "alignments" / f"{label}.raw.fasta"
    aligned_fasta = results_dir / "alignments" / f"{label}.aligned.fasta"
    write_fasta(combined_seqs, aln_fasta)
    run_mafft(aln_fasta, aligned_fasta)
    aligned = read_fasta(aligned_fasta)

    dist = p_distance_matrix(aligned)
    intra_max = max_intra_species(dist, species_of)

    rows = []
    for name in isolate_names:
        if name not in aligned:
            continue
        result = classify_isolate(name, dist, species_of, intra_max)
        if result is None:
            continue
        result["probennummer"] = name.split("__", 1)[1]
        result["locus"] = label
        rows.append(result)
    df = pd.DataFrame(rows)

    figures_dir = results_dir / "figures"
    plot_alignment_heatmap(aligned, species_of, figures_dir / f"heatmap_{label}",
                            f"Isolates vs. reference species — {label}")
    n_tips = len(aligned)
    if n_tips <= 150:
        from mlsa.align import run_iqtree
        treefile = run_iqtree(aligned_fasta, results_dir / "alignments" / f"{label}.tree")
        plot_tree(treefile, species_of, figures_dir / f"tree_{label}", f"Tree — {label}")
    else:
        log.info("Skipping tree image for %s (%d sequences, too many to render legibly)", label, n_tips)
    return df


def main(results_dir: Path = RESULTS, main1_results: Path = MAIN1_RESULTS,
         exclude_accessions: Iterable[str] = ()) -> None:
    exclude_accessions = frozenset(exclude_accessions)
    results_dir.mkdir(parents=True, exist_ok=True)
    tnr_to_pnr = load_tnr_to_probennummer()
    log.info("Known TNRs from %s: %d (for %d isolates)",
             SCREENING_MAP.name, len(tnr_to_pnr), len(set(tnr_to_pnr.values())))

    by_locus_pnr = collect_isolate_reads(tnr_to_pnr)

    isolate_seqs: dict[str, dict[str, str]] = {}
    for locus in LOCI:
        refs = load_reference_sequences(locus, main1_results, exclude_accessions)
        ref_example = next(iter(refs.values())) if refs else None

        seqs = dict(refs)
        n_ok = n_fail = 0
        for pnr, paths in by_locus_pnr[locus].items():
            rep = pick_representative(paths)
            if rep is None:
                n_fail += 1
                continue
            seq, mean_q, source = rep
            if ref_example:
                seq = orient_to_reference(seq, ref_example)
            seqs[f"isolate__{pnr}"] = seq
            n_ok += 1
        log.info("%s: %d isolates with a usable representative read, %d with none (of %d isolates seen)",
                  locus, n_ok, n_fail, len(by_locus_pnr[locus]))
        isolate_seqs[locus] = seqs

    species_of = {name: species_of_name(name) for seqs in isolate_seqs.values() for name in seqs}

    all_result_dfs = []
    for locus in LOCI:
        isolate_names = [n for n in isolate_seqs[locus] if species_of[n] == "isolate"]
        df = run_locus_or_combo(locus, isolate_seqs[locus], species_of, isolate_names, results_dir)
        all_result_dfs.append(df)

    # hsp65 + 16S concatenated, for isolates (and reference genomes) that have both.
    both_isolate_pnrs = (
        {n.split("__", 1)[1] for n in isolate_seqs["hsp65"] if species_of[n] == "isolate"}
        & {n.split("__", 1)[1] for n in isolate_seqs["16S"] if species_of[n] == "isolate"}
    )
    if both_isolate_pnrs:
        # Re-align each locus alone first (needed to get equal-length columns to concatenate),
        # reusing the alignments just written by run_locus_or_combo.
        hsp65_aligned = read_fasta(results_dir / "alignments" / "hsp65.aligned.fasta")
        s16_aligned = read_fasta(results_dir / "alignments" / "16S.aligned.fasta")
        common_refs = set(hsp65_aligned) & set(s16_aligned) - {f"isolate__{t}" for t in both_isolate_pnrs}
        combined = {}
        for name in common_refs | {f"isolate__{t}" for t in both_isolate_pnrs}:
            if name in hsp65_aligned and name in s16_aligned:
                combined[name] = hsp65_aligned[name] + s16_aligned[name]
        isolate_names = [f"isolate__{t}" for t in both_isolate_pnrs]
        df_both = run_locus_or_combo("hsp65+16S", combined, species_of, isolate_names, results_dir)
        all_result_dfs.append(df_both)
    else:
        log.info("No isolates have both hsp65 and 16S representative reads; skipping combined analysis")
        df_both = pd.DataFrame()

    results_df = pd.concat(all_result_dfs, ignore_index=True) if all_result_dfs else pd.DataFrame()
    results_df.to_csv(results_dir / "isolate_classification.tsv", sep="\t", index=False)

    summary_lines = ["Differentiation summary (unambiguous / total isolates):"]
    for label, df in [("hsp65", all_result_dfs[0]), ("16S", all_result_dfs[1]), ("hsp65+16S", df_both)]:
        if len(df):
            summary_lines.append(f"  {label}: {df['unambiguous'].sum()} / {len(df)}")
        else:
            summary_lines.append(f"  {label}: no data")

    if len(df_both):
        ambiguous_both = df_both[~df_both["unambiguous"]]
        recs = recommend_loci(ambiguous_both, main1_results / "single_locus_pair_separation.tsv")
        if len(recs):
            recs.to_csv(results_dir / "recommended_additional_loci.tsv", sep="\t", index=False)
            top = recs["recommended_locus"].value_counts()
            summary_lines.append("")
            summary_lines.append("Most commonly recommended additional locus for isolates still "
                                  "ambiguous with hsp65+16S:")
            for locus, n in top.items():
                summary_lines.append(f"  {locus}: {n} isolate(s)")

    summary_text = "\n".join(summary_lines)
    log.info("\n%s", summary_text)
    if exclude_accessions:
        summary_text += "\n\nReference genomes excluded: " + ", ".join(sorted(exclude_accessions))
    with open(results_dir / "SUMMARY.txt", "w") as fh:
        fh.write(summary_text + "\n")


if __name__ == "__main__":
    main()
