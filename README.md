# mlsa-kansasii

Can the *Mycobacterium kansasii* species complex (kansasii, persicum,
pseudokansasii, innocens, attenuatum, ostraviense, gastri) be differentiated
by Sanger sequencing / MLSA with as few loci as possible, and how far does
the lab's existing hsp65 + 16S Sanger data already get?

## Layout

- `scripts/link_sanger_kansasii.sh` — symlinks the `.ab1` Sanger traces of
  the TNR isolates in `data/imm/screening_map_link.csv` into
  `data/sanger/seq_kansasii/`.
- `scripts/link_gtdb_kansasii_complex.sh` — downloads/filters GTDB r232
  metadata for the 7 complex species and symlinks their genomes
  (`.fna`/`genomic.gff`) into
  `data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii/<species>/`.
- `mlsa/` — shared Python package: locus extraction (GFF gene symbols; 16S
  and ITS via rRNA features; hsp65 via in-silico PCR with the lab's
  TB-11/TB-12(w) primers parsed from the protocol `.docx`), AB1 parsing +
  quality trimming, mafft/iqtree wrappers (via Singularity), and plotting.
- `main1_locus_discovery/` — finds the smallest locus combination that fully
  separates all 7 species in the GTDB reference genomes.
- `main2_sanger_differentiation/` — checks how far the isolates' existing
  hsp65/16S Sanger reads actually resolve species, and recommends which
  additional locus (from main1's candidates) would resolve the ambiguous
  ones.
- `genome_identity_check/` — all-vs-all ANI (skani) between the reference
  genomes, to check whether genomes with atypical marker genes are really
  another species (see LIT.md).

## Setup

```bash
conda env create -f environment.yml
conda activate env_mlsa
bash Singularity/pull_singularity_img.sh /shares/sander.imm.uzh/software/pipelines/IMMense/IMMense_dependencies/containers
```

## Run order

```bash
bash scripts/link_sanger_kansasii.sh
bash scripts/link_gtdb_kansasii_complex.sh
sbatch main1_locus_discovery/submit_locus_discovery.sbatch
sbatch main2_sanger_differentiation/submit_sanger_differentiation.sbatch   # after main1 finishes
sbatch genome_identity_check/submit_genome_identity_check.sbatch           # independent of main1/main2
```

Outputs are written to the shares, not into the repo:
`/shares/sander.imm.uzh/MM/kansasii/output/mlsa/main1_locus_discovery/` and
`/shares/sander.imm.uzh/MM/kansasii/output/mlsa/main2_sanger_differentiation/`
(tables + `figures/`), and `output/mlsa/genome_identity_check/`.

## How outputs are used downstream

| Step | Produces | Used by |
|---|---|---|
| `scripts/link_sanger_kansasii.sh` | `.ab1` symlinks in `data/sanger/seq_kansasii/` | main2 |
| `scripts/link_gtdb_kansasii_complex.sh` | genome `.fna`/`.gff` symlinks in `data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii/` | main1, main2 (only if main1's output is missing), genome_identity_check |
| main1 | `alignments/<locus>.raw.fasta` | main2 (hsp65 and 16S reference sequences) |
| main1 | `single_locus_pair_separation.tsv` | main2 (`recommended_additional_loci.tsv`) |
| main1 | `combo_results.tsv`, `winning_combo/`, `figures/`, `SUMMARY.txt` | final results |
| main2 | `isolate_classification.tsv`, `recommended_additional_loci.tsv`, `figures/`, `SUMMARY.txt` | final results |
| `*_excluded` variants | the same, in their own output folders | main2_excluded reads main1_excluded; otherwise final results |
| genome_identity_check | `skani_ani.tsv`, `skani_long.tsv`, `species_ani_summary.tsv` | no code; the ANI numbers are cited by hand in LIT.md |

The genome_identity_check results inform decisions about which reference
genomes to trust or exclude. No pipeline step reads them.

## How main2 classifies isolates (`isolate_classification.tsv`)

Implemented in `main2_sanger_differentiation/run_sanger_differentiation.py`.
The table has one row per isolate (`probennummer`) per `locus` (`hsp65`,
`16S`, `hsp65+16S`).

1. **One sequence per isolate and locus.** Reads under all of an isolate's
   TNR columns are pooled. Every `.ab1` read is quality-trimmed (Mott, reads
   shorter than 100 bp are dropped) and the longest one is kept, with ties
   going to higher mean quality. Forward and reverse reads are not merged
   into a consensus. The read is flipped if needed to match the reference
   strand.
2. **Alignment.** The isolate reads are aligned with MAFFT together with the
   reference sequences of the 7 species from main1
   (`main1_locus_discovery/alignments/<locus>.raw.fasta`, or extracted again
   if main1's output is missing). Distance is the uncorrected p-distance: the
   fraction of positions that differ, counting only positions where neither
   sequence has a gap.
3. **Tolerance per species.** A species' tolerance is the largest distance
   between any two of its own reference genomes. It is 0 if the species has
   only one reference.
4. **Classification.** For each species, take the isolate's distance to that
   species' closest reference. The closest and second-closest species give
   `nearest_species`/`nearest_dist` and `second_species`/`second_dist`.
   `margin = second_dist - nearest_dist`, `tolerance` is the nearest
   species' tolerance, and **`unambiguous = margin > tolerance`**. So the
   second-closest species has to be further away than the nearest species'
   own references are from each other.
5. **hsp65+16S.** Only isolates with both reads get this row. The two
   single-locus alignments are joined end to end, keeping only references
   that have both loci, and steps 3–4 are repeated.

Isolates still ambiguous with hsp65+16S get a suggested extra locus in
`recommended_additional_loci.tsv`: the main1 candidate locus that best
separates their nearest and second-closest species.

Things to keep in mind when reading the table:

- The tolerance comes from the reference genomes of the nearest species.
  Diverse species such as kansasii (hsp65 tolerance 0.025) make isolates
  ambiguous even when they sit right next to a reference.
- Gaps are skipped, so a short read that overlaps little of the alignment is
  judged on only a few positions.
- Reference 16S/ITS: some GTDB assemblies contain contaminant contigs with
  their own 16S. `mlsa.loci.extract_16s` skips partial copies and keeps the
  copy closest to the type-strain 16S (see LIT.md, "Extraction artifact").

The TNR isolates will also be Illumina-sequenced and speciated with the
`kansasii` branch of
[immensekansasii](https://gitlab.uzh.ch/appliedmicrobiologyresearch/immense.git)
(`../immensekansasii`); this repo answers the same species-complex question
computationally ahead of that NGS data being available.
