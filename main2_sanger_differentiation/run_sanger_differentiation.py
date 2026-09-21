#!/usr/bin/env python
"""Main program 2: using the Sanger hsp65 + 16S reads the lab already has for
the TNR isolates in data/imm/screening_map.csv (symlinked into
data/sanger/seq_kansasii by scripts/link_sanger_kansasii.sh), assess how far
species-level differentiation actually gets today, and which additional
locus/loci (from main1_locus_discovery's candidate set) would need to be
sequenced to resolve the isolates that remain ambiguous.

For each isolate + locus, all matching .ab1 reads are quality-trimmed and the
best (longest, then highest quality) read is taken as that isolate's
representative sequence (no fwd/rev consensus assembly in this first pass).
Representative sequences are oriented against a reference sequence, aligned
together with the 7-species reference sequences from main1 (reused from
results/main1_locus_discovery/alignments/*.raw.fasta when available, else
re-extracted directly), and each isolate is classified as unambiguously
nested in one species' cluster, or ambiguous, using the same DNA-barcoding
gap logic as main1 (isolate-to-species distance vs. that species' own
observed diversity among references).

Run via submit_sanger_differentiation.sbatch (requires env_mlsa active and
the mafft Singularity container pulled). Best run after main1, but falls
back to extracting reference sequences itself if main1's output is missing.
"""
from __future__ import annotations

import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

from mlsa import SANGER_LINK_DIR, SCREENING_MAP, SPECIES, GTDB_MLSA_ROOT  # noqa: E402
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


def load_known_tnrs() -> set[str]:
    with open(SCREENING_MAP) as fh:
        reader = csv.DictReader(fh)
        return {row["TNR"] for row in reader if row["TNR"]}


def collect_isolate_reads(known_tnrs: set[str]) -> dict[str, dict[str, list[Path]]]:
    """locus -> {tnr: [ab1 paths]}"""
    by_locus_tnr: dict[str, dict[str, list[Path]]] = {locus: defaultdict(list) for locus in LOCI}
    for path in sorted(SANGER_LINK_DIR.glob("*.ab1")):
        locus = classify_locus(path.name)
        if locus is None:
            continue
        tnr = extract_tnr(path.name, known_tnrs)
        if tnr is None:
            continue
        by_locus_tnr[locus][tnr].append(path)
    return by_locus_tnr


def pick_representative(paths: list[Path]) -> tuple[str, float, Path] | None:
    """Quality-trim every candidate read for a (tnr, locus) and keep the
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


def load_reference_sequences(locus: str) -> dict[str, str]:
    """Reuse main1's already-extracted reference sequences when available;
    otherwise extract them directly so main2 can run standalone."""
    raw_fasta = MAIN1_RESULTS / "alignments" / f"{locus}.raw.fasta"
    if raw_fasta.exists():
        log.info("Reusing main1 reference sequences for %s from %s", locus, raw_fasta)
        return read_fasta(raw_fasta)

    log.warning("main1 output for %s not found at %s; extracting references directly", locus, raw_fasta)
    genomes = discover_genomes(GTDB_MLSA_ROOT, SPECIES)
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
                      intra_max: dict[str, float]) -> dict:
    per_species_min: dict[str, float] = {}
    for ref_name, sp in species_of.items():
        if sp not in SPECIES or ref_name not in dist.columns:
            continue
        d = dist.loc[isolate_name, ref_name]
        per_species_min[sp] = min(per_species_min.get(sp, float("inf")), d)
    ranked = sorted(per_species_min.items(), key=lambda kv: kv[1])
    best_sp, best_d = ranked[0]
    second_sp, second_d = ranked[1] if len(ranked) > 1 else (None, float("inf"))
    tolerance = intra_max.get(best_sp, 0.0)
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
        sp_a, sp_b = sorted([row["nearest_species"], row["second_species"]])
        candidates = pair_sep[
            ((pair_sep["species_a"] == sp_a) & (pair_sep["species_b"] == sp_b))
            & pair_sep["separated"]
        ].sort_values("gap_margin", ascending=False)
        recommended = candidates.iloc[0]["locus"] if len(candidates) else "none of the candidate loci"
        rows.append({"tnr": row["tnr"], "species_a": sp_a, "species_b": sp_b, "recommended_locus": recommended})
    return pd.DataFrame(rows)


def run_locus_or_combo(label: str, combined_seqs: dict[str, str], species_of: dict[str, str],
                        isolate_names: list[str]) -> pd.DataFrame:
    aln_fasta = RESULTS / "alignments" / f"{label}.raw.fasta"
    aligned_fasta = RESULTS / "alignments" / f"{label}.aligned.fasta"
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
        result["tnr"] = name.split("__", 1)[1]
        result["locus"] = label
        rows.append(result)
    df = pd.DataFrame(rows)

    figures_dir = RESULTS / "figures"
    plot_alignment_heatmap(aligned, species_of, figures_dir / f"heatmap_{label}",
                            f"Isolates vs. reference species — {label}")
    n_tips = len(aligned)
    if n_tips <= 150:
        from mlsa.align import run_iqtree
        treefile = run_iqtree(aligned_fasta, RESULTS / "alignments" / f"{label}.tree")
        plot_tree(treefile, species_of, figures_dir / f"tree_{label}", f"Tree — {label}")
    else:
        log.info("Skipping tree image for %s (%d sequences, too many to render legibly)", label, n_tips)
    return df


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    known_tnrs = load_known_tnrs()
    log.info("Known TNRs from screening_map.csv: %d", len(known_tnrs))

    by_locus_tnr = collect_isolate_reads(known_tnrs)

    isolate_seqs: dict[str, dict[str, str]] = {}
    for locus in LOCI:
        refs = load_reference_sequences(locus)
        ref_example = next(iter(refs.values())) if refs else None

        seqs = dict(refs)
        n_ok = n_fail = 0
        for tnr, paths in by_locus_tnr[locus].items():
            rep = pick_representative(paths)
            if rep is None:
                n_fail += 1
                continue
            seq, mean_q, source = rep
            if ref_example:
                seq = orient_to_reference(seq, ref_example)
            seqs[f"isolate__{tnr}"] = seq
            n_ok += 1
        log.info("%s: %d isolates with a usable representative read, %d with none (of %d TNRs seen)",
                  locus, n_ok, n_fail, len(by_locus_tnr[locus]))
        isolate_seqs[locus] = seqs

    species_of = {name: species_of_name(name) for seqs in isolate_seqs.values() for name in seqs}

    all_result_dfs = []
    for locus in LOCI:
        isolate_names = [n for n in isolate_seqs[locus] if species_of[n] == "isolate"]
        df = run_locus_or_combo(locus, isolate_seqs[locus], species_of, isolate_names)
        all_result_dfs.append(df)

    # hsp65 + 16S concatenated, for isolates (and reference genomes) that have both.
    both_isolate_tnrs = (
        {n.split("__", 1)[1] for n in isolate_seqs["hsp65"] if species_of[n] == "isolate"}
        & {n.split("__", 1)[1] for n in isolate_seqs["16S"] if species_of[n] == "isolate"}
    )
    if both_isolate_tnrs:
        # Re-align each locus alone first (needed to get equal-length columns to concatenate),
        # reusing the alignments just written by run_locus_or_combo.
        hsp65_aligned = read_fasta(RESULTS / "alignments" / "hsp65.aligned.fasta")
        s16_aligned = read_fasta(RESULTS / "alignments" / "16S.aligned.fasta")
        common_refs = set(hsp65_aligned) & set(s16_aligned) - {f"isolate__{t}" for t in both_isolate_tnrs}
        combined = {}
        for name in common_refs | {f"isolate__{t}" for t in both_isolate_tnrs}:
            if name in hsp65_aligned and name in s16_aligned:
                combined[name] = hsp65_aligned[name] + s16_aligned[name]
        isolate_names = [f"isolate__{t}" for t in both_isolate_tnrs]
        df_both = run_locus_or_combo("hsp65+16S", combined, species_of, isolate_names)
        all_result_dfs.append(df_both)
    else:
        log.info("No isolates have both hsp65 and 16S representative reads; skipping combined analysis")
        df_both = pd.DataFrame()

    results_df = pd.concat(all_result_dfs, ignore_index=True) if all_result_dfs else pd.DataFrame()
    results_df.to_csv(RESULTS / "isolate_classification.tsv", sep="\t", index=False)

    summary_lines = ["Differentiation summary (unambiguous / total isolates):"]
    for label, df in [("hsp65", all_result_dfs[0]), ("16S", all_result_dfs[1]), ("hsp65+16S", df_both)]:
        if len(df):
            summary_lines.append(f"  {label}: {df['unambiguous'].sum()} / {len(df)}")
        else:
            summary_lines.append(f"  {label}: no data")

    if len(df_both):
        ambiguous_both = df_both[~df_both["unambiguous"]]
        recs = recommend_loci(ambiguous_both, MAIN1_RESULTS / "single_locus_pair_separation.tsv")
        if len(recs):
            recs.to_csv(RESULTS / "recommended_additional_loci.tsv", sep="\t", index=False)
            top = recs["recommended_locus"].value_counts()
            summary_lines.append("")
            summary_lines.append("Most commonly recommended additional locus for isolates still "
                                  "ambiguous with hsp65+16S:")
            for locus, n in top.items():
                summary_lines.append(f"  {locus}: {n} isolate(s)")

    summary_text = "\n".join(summary_lines)
    log.info("\n%s", summary_text)
    with open(RESULTS / "SUMMARY.txt", "w") as fh:
        fh.write(summary_text + "\n")


if __name__ == "__main__":
    main()
