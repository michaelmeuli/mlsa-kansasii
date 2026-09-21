#!/usr/bin/env python
"""Main program 1: find the smallest set of loci that, extracted from the
GTDB r232 reference genomes, fully separates the 7 species of the
M. kansasii complex (kansasii, persicum, pseudokansasii, innocens,
attenuatum, ostraviense, gastri).

Candidate loci: rpoB, gyrA, gyrB, recA, secA1, tuf (single-copy GFF gene
symbols), 16S (full-length rRNA feature), ITS (16S-23S spacer), hsp65
(primer-based in-silico PCR matching the lab's TB-11/TB-12(w) assay).

For each candidate locus and then for combinations of loci (smallest first),
species are considered fully separated when every pair of species shows a
DNA-barcoding "gap": the closest pair of sequences from two different
species is farther apart than the most divergent pair of sequences within
either species (mlsa.align.resolves_all_species). Exhaustive search is used
up to combination size 4 (choose(9,4) = 126, trivial); if nothing resolves
by then, greedy forward selection continues by always adding the locus that
most reduces the number of unseparated species pairs.

Run via submit_locus_discovery.sbatch on the cluster (requires env_mlsa
active and the mafft/iqtree Singularity containers pulled).
"""
from __future__ import annotations

import itertools
import logging
import sys
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

from mlsa import SPECIES, GTDB_MLSA_ROOT  # noqa: E402
from mlsa.align import (  # noqa: E402
    barcoding_gap_check,
    p_distance_matrix,
    read_fasta,
    resolves_all_species,
    run_iqtree,
    run_mafft,
    species_sort_key,
    write_fasta,
)
from mlsa.loci import CANDIDATE_LOCI, discover_genomes, extract_all_loci  # noqa: E402
from mlsa.plotting import plot_alignment_heatmap, plot_tree  # noqa: E402

RESULTS = Path("/shares/sander.imm.uzh/MM/kansasii/output") / "mlsa" / "main1_locus_discovery"
MAX_EXHAUSTIVE_SIZE = 4

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("locus_discovery")


def genome_name(genome) -> str:
    return f"{genome.species}__{genome.accession}"


def extract_loci_for_all_genomes(genomes) -> dict[str, dict[str, str]]:
    """locus -> {genome_name: raw (unaligned) sequence}, skipping genomes
    where extraction failed for that locus."""
    per_locus: dict[str, dict[str, str]] = {locus: {} for locus in CANDIDATE_LOCI}
    for g in genomes:
        seqs = extract_all_loci(g)
        for locus, seq in seqs.items():
            if seq:
                per_locus[locus][genome_name(g)] = seq
    return per_locus


def write_coverage_table(per_locus: dict[str, dict[str, str]], species_of: dict[str, str],
                          out_path: Path) -> None:
    rows = []
    for locus, seqs in per_locus.items():
        species_present = {species_of[n] for n in seqs}
        rows.append({
            "locus": locus,
            "n_genomes": len(seqs),
            "n_species": len(species_present),
            "species_missing": ",".join(sorted(set(SPECIES) - species_present)) or "-",
        })
    df = pd.DataFrame(rows).sort_values("locus")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    log.info("Wrote locus coverage table to %s", out_path)


def align_each_locus(per_locus: dict[str, dict[str, str]], usable_loci: list[str],
                      results_dir: Path = RESULTS) -> dict[str, dict[str, str]]:
    aligned: dict[str, dict[str, str]] = {}
    for locus in usable_loci:
        raw_fasta = results_dir / "alignments" / f"{locus}.raw.fasta"
        aln_fasta = results_dir / "alignments" / f"{locus}.aligned.fasta"
        write_fasta(per_locus[locus], raw_fasta)
        run_mafft(raw_fasta, aln_fasta)
        aligned[locus] = read_fasta(aln_fasta)
        log.info("Aligned locus %s: %d sequences, length %d", locus, len(aligned[locus]),
                  len(next(iter(aligned[locus].values()))))
    return aligned


def concat_combo(aligned: dict[str, dict[str, str]], combo: tuple[str, ...]) -> dict[str, str]:
    common_names = set.intersection(*(set(aligned[locus]) for locus in combo))
    return {name: "".join(aligned[locus][name] for locus in combo) for name in sorted(common_names, key=species_sort_key)}


def evaluate_combo(aligned: dict[str, dict[str, str]], combo: tuple[str, ...],
                    species_of: dict[str, str]) -> dict:
    concat = concat_combo(aligned, combo)
    species_present = {species_of[n] for n in concat}
    missing_species = set(SPECIES) - species_present
    result = {
        "loci": "+".join(combo),
        "size": len(combo),
        "n_genomes": len(concat),
        "species_missing": ",".join(sorted(missing_species)) or "-",
        "resolves_all": False,
        "min_gap_margin": float("nan"),
        "mean_gap_margin": float("nan"),
        "n_pairs_separated": 0,
        "n_pairs_total": 0,
    }
    if missing_species or len(concat) < 2:
        return result
    dist = p_distance_matrix(concat)
    resolved, gap_df = resolves_all_species(dist, {n: species_of[n] for n in concat})
    result.update({
        "resolves_all": resolved,
        "min_gap_margin": gap_df["gap_margin"].min(),
        "mean_gap_margin": gap_df["gap_margin"].mean(),
        "n_pairs_separated": int(gap_df["separated"].sum()),
        "n_pairs_total": len(gap_df),
    })
    return result


def write_single_locus_pair_table(aligned: dict[str, dict[str, str]], usable_loci: list[str],
                                   species_of: dict[str, str], out_path: Path) -> pd.DataFrame:
    """Per-locus, per-species-pair separation detail (used by
    main2_sanger_differentiation to recommend which single additional locus
    would resolve a given ambiguous species pair among real isolates)."""
    frames = []
    for locus in usable_loci:
        dist = p_distance_matrix(aligned[locus])
        gap_df = barcoding_gap_check(dist, {n: species_of[n] for n in aligned[locus]})
        gap_df.insert(0, "locus", locus)
        frames.append(gap_df)
    df = pd.concat(frames, ignore_index=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    log.info("Wrote per-locus species-pair separation table to %s", out_path)
    return df


def search_minimal_combo(aligned: dict[str, dict[str, str]], usable_loci: list[str],
                          species_of: dict[str, str]) -> tuple[list[dict], list[str] | None]:
    all_results = []
    max_size = min(MAX_EXHAUSTIVE_SIZE, len(usable_loci))
    for size in range(1, max_size + 1):
        size_results = []
        for combo in itertools.combinations(sorted(usable_loci), size):
            res = evaluate_combo(aligned, combo, species_of)
            size_results.append(res)
        all_results.extend(size_results)
        resolving = [r for r in size_results if r["resolves_all"]]
        if resolving:
            best = max(resolving, key=lambda r: r["mean_gap_margin"])
            log.info("Smallest resolving combination found at size %d: %s", size, best["loci"])
            return all_results, best["loci"].split("+")

    # Greedy forward selection fallback beyond MAX_EXHAUSTIVE_SIZE.
    log.info("No combination up to size %d resolves all species; falling back to greedy selection", max_size)
    chosen: list[str] = []
    remaining = list(usable_loci)
    best_result = None
    while remaining:
        candidates = []
        for locus in remaining:
            combo = tuple(sorted(chosen + [locus]))
            res = evaluate_combo(aligned, combo, species_of)
            candidates.append((res, locus))
        candidates.sort(key=lambda cr: (-cr[0]["n_pairs_separated"], -cr[0]["mean_gap_margin"]))
        best_res, best_locus = candidates[0]
        all_results.append(best_res)
        chosen.append(best_locus)
        remaining.remove(best_locus)
        best_result = best_res
        log.info("Greedy step: added %s -> %d/%d pairs separated", best_locus,
                  best_res["n_pairs_separated"], best_res["n_pairs_total"])
        if best_res["resolves_all"]:
            return all_results, chosen
    log.warning("No combination of the %d candidate loci fully resolves all 7 species; "
                "best effort uses all of them (%d/%d pairs separated)",
                len(usable_loci), best_result["n_pairs_separated"], best_result["n_pairs_total"])
    return all_results, None


def main(results_dir: Path = RESULTS, exclude_accessions: Iterable[str] = ()) -> None:
    exclude_accessions = frozenset(exclude_accessions)
    results_dir.mkdir(parents=True, exist_ok=True)
    genomes = discover_genomes(GTDB_MLSA_ROOT, SPECIES)
    log.info("Discovered %d reference genomes across %d species", len(genomes), len(SPECIES))
    if exclude_accessions:
        present = {g.accession for g in genomes}
        not_found = sorted(set(exclude_accessions) - present)
        if not_found:
            log.warning("Requested exclusions not found among discovered genomes: %s", not_found)
        genomes = [g for g in genomes if g.accession not in exclude_accessions]
        log.info("Excluded %d genome(s) (%s); %d remain",
                 len(present & exclude_accessions), ", ".join(sorted(exclude_accessions)), len(genomes))
    species_of = {genome_name(g): g.species for g in genomes}

    per_locus = extract_loci_for_all_genomes(genomes)
    write_coverage_table(per_locus, species_of, results_dir / "locus_coverage.tsv")

    usable_loci = [locus for locus, seqs in per_locus.items()
                   if {species_of[n] for n in seqs} == set(SPECIES)]
    dropped = sorted(set(CANDIDATE_LOCI) - set(usable_loci))
    if dropped:
        log.warning("Dropping candidate loci with incomplete species coverage: %s", dropped)
    log.info("Usable candidate loci (all 7 species represented): %s", usable_loci)

    aligned = align_each_locus(per_locus, usable_loci, results_dir)
    write_single_locus_pair_table(aligned, usable_loci, species_of,
                                   results_dir / "single_locus_pair_separation.tsv")

    all_results, winning_combo = search_minimal_combo(aligned, usable_loci, species_of)
    results_df = pd.DataFrame(all_results).sort_values(["size", "resolves_all"], ascending=[True, False])
    results_df.to_csv(results_dir / "combo_results.tsv", sep="\t", index=False)
    log.info("Wrote %d evaluated combinations to %s", len(results_df), results_dir / "combo_results.tsv")

    if winning_combo is None:
        log.warning("Finished without a fully resolving locus combination.")
        return

    log.info("Winning minimal locus set: %s", "+".join(winning_combo))
    concat = concat_combo(aligned, tuple(winning_combo))
    winner_dir = results_dir / "winning_combo"
    winner_fasta = winner_dir / f"{'+'.join(winning_combo)}.aligned.fasta"
    write_fasta(concat, winner_fasta)

    treefile = run_iqtree(winner_fasta, winner_dir / "tree")
    figures_dir = results_dir / "figures"
    combo_label = "+".join(winning_combo)
    plot_alignment_heatmap(concat, species_of, figures_dir / f"heatmap_{combo_label}",
                            f"Sequence differences across the M. kansasii complex — {combo_label}")
    plot_tree(treefile, species_of, figures_dir / f"tree_{combo_label}",
              f"Maximum-likelihood tree — {combo_label}")

    with open(results_dir / "SUMMARY.txt", "w") as fh:
        fh.write(f"Minimal locus set fully separating all 7 species: {combo_label}\n")
        fh.write(f"Genomes used: {len(concat)} / {len(genomes)}\n")
        if exclude_accessions:
            fh.write(f"Excluded genomes: {', '.join(sorted(exclude_accessions))}\n")
        fh.write("See combo_results.tsv for every tested combination and "
                  "locus_coverage.tsv for per-locus extraction coverage.\n")
    log.info("Done. Summary written to %s", results_dir / "SUMMARY.txt")


if __name__ == "__main__":
    main()
