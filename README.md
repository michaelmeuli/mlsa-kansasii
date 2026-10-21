# mlsa-kansasii

Can the *Mycobacterium kansasii* species complex (kansasii, persicum,
pseudokansasii, innocens, attenuatum, ostraviense, gastri) be differentiated
by Sanger sequencing / MLSA with as few loci as possible, and how far does
the lab's existing hsp65 + 16S Sanger data already get?

## Layout

- `scripts/01_download_filter_gtdb.sh` — downloads GTDB r232 metadata (if missing) and lists the
  accessions of the 7 complex species.
- `scripts/03_link_sanger_kansasii.sh` — symlinks the `.ab1` Sanger traces of
  the TNR isolates in `data/imm/screening_map_link.csv` into
  `data/sanger/seq_kansasii/`.
- `scripts/02_link_gtdb_kansasii_complex.sh` — symlinks the genomes of that list
  (`.fna`/`genomic.gff`, downloaded by immensekansasii step 1) into
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
- `main3_reference_alignment/` — species ID of Sanger reads (`.ab1`/FASTA) by plain alignment to one
  reference genome per species (`kansasii_complex_gtdb_representatives/`), no MLSA. One row per read in
  `reference_alignment.tsv` (no PDFs: sanger-microsynth makes and emails the per-read alignment PDF).
  Logic in `mlsa/refalign.py` (copied to sanger-microsynth).
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

Run from the repo root. Inputs: `data/imm/screening_map_link.csv` (immensekansasii step 2) and the
GTDB genomes downloaded by immensekansasii step 1. `scripts/` is numbered in the order it is called;
the `mainN` programs are numbered by their own order:

| # | Script | Needs | Writes |
|---|---|---|---|
| 1 | `scripts/01_download_filter_gtdb.sh` | - | accession list of the 7 species (`data/lit/gtdb/gtdb232/`) |
| 2 | `scripts/02_link_gtdb_kansasii_complex.sh` | 1, genomes | `data/gtdb_genomes/.../mlsa-kansasii/<species>/` |
| 3 | `scripts/03_link_sanger_kansasii.sh` | link table | `data/sanger/seq_kansasii/` |
| 4 | `sbatch main1_locus_discovery/submit_locus_discovery.sbatch` (+ `main1_locus_discovery_excluded`) | 2 | `output/mlsa/main1_*` |
| 5 | `python scripts/04_make_hsp65_references.py` | main1 | `data/gtdb_genomes/.../kansasii_complex_hsp65_amplicons/` |
| 6 | `sbatch main2_sanger_differentiation/submit_sanger_differentiation.sbatch` (+ `_excluded`) | 3, main1 | `output/mlsa/main2_*` |
| 7 | `sbatch main3_reference_alignment/submit_reference_alignment.sbatch`, then `submit_hsp65_amplicon_alignment.sbatch` | 3 (and 5 for the second) | `output/mlsa/main3_*` |
| 8 | `sbatch genome_identity_check/submit_genome_identity_check.sbatch` | 2 | `output/mlsa/genome_identity_check/` |
| 9 | immensekansasii `06_screening_map_results.py` | 6, 7 | `data/imm/screening_map_results.csv` |
| 10 | `sbatch main4_hsp65_vs_wgs/submit_hsp65_vs_wgs.sbatch` | 6, 7, 9 | `output/mlsa/main4_hsp65_vs_wgs/` |

`immensekansasii/scripts/run_downstream.sh` submits all of this (and the other repos' downstream
steps) as one dependency chain.

Outputs are written to the shares, not into the repo:
`/shares/sander.imm.uzh/MM/kansasii/output/mlsa/main1_locus_discovery/` and
`/shares/sander.imm.uzh/MM/kansasii/output/mlsa/main2_sanger_differentiation/`
(tables + `figures/`), and `output/mlsa/genome_identity_check/`.

To download all results to a Windows machine (PowerShell):

```powershell
New-Item -ItemType Directory -Path "$env:USERPROFILE\kansasii_C\downloads\" -Force
scp -r mimeul@cluster.s3it.uzh.ch:/shares/sander.imm.uzh/MM/kansasii/output/* "$env:USERPROFILE\kansasii_C\downloads\"
```

## Plain alignment to references vs MLSA

`python main3_reference_alignment/run_reference_alignment.py <reads or dirs>` (about 3 s per read).
Call = closest reference by differences over the read: `ok` (&ge;98% identity, &ge;2
differences fewer than the closest other species), `ambiguous`, `divergent` (&lt;98%, e.g. new lineage, introgressed
locus, or bad read). By default the 7 representatives; for hsp65 reads use `--ref-dir` with the 72 per-genome
hsp65 amplicons (`scripts/04_make_hsp65_references.py`, job `submit_hsp65_amplicon_alignment.sbatch`), which also flags
reads matching the atypical-hsp65 kansasii genomes (`atypical_hsp65` column). main4 benchmarks this against
whole-genome species: 82 of 83 calls correct on 92 isolates (the miss is a mixed culture). Compared with main1/main2:

- Plain alignment needs no locus extraction, no annotation and no tolerance estimate, works on any read,
  and takes seconds. It also never sees the 72-genome diversity, so the 3 Korean *M. kansasii* genomes
  (LIT.md) cannot disturb it. But that is blindness, not robustness: an isolate carrying their
  persicum-like hsp65 would be called persicum, with the representatives giving no hint of the problem.
- main2 (MLSA-style) measures within-species diversity from all references and so says when a call is
  not supported (`unambiguous`), but is hsp65/16S-limited and needs main1's alignments.
- Neither fixes a single-locus limit: hsp65 and 16S alone cannot separate every species (16S reads
  tie between several references). More loci (gyrA, see main2 recommendations) or genome data are needed.

## IUPAC-aware alignment (opt-in, tested 2026-10-21: no gain)

Sanger reads carry IUPAC codes (Y, R, K, M, S, W, N) where the trace has mixed peaks. By default main3 *ignores*
them: `_diff` in `mlsa/refalign.py` compares only columns where read and reference both have A/C/G/T, and any
non-ACGT reference base becomes `N`. The references (genome amplicons, locus FASTAs) are pure ACGT, so
IUPAC-aware scoring on the read side alone has little to act on. The test therefore also built references that
carry ambiguity: one consensus per species from `output/mlsa/main1_locus_discovery/alignments/hsp65.aligned.fasta`
where every column with more than one base becomes the IUPAC code of those bases.

`run_reference_alignment.py --iupac` (default off, so the default output is unchanged):
- the aligner scores identical codes 2, overlapping codes (R vs A) 1, anything with N 0, disjoint codes -3;
- ambiguity codes in the reference are kept;
- a column is a difference only if the two codes share no base (R vs A matches, R vs C differs); N columns are not compared.

`scripts/05_make_iupac_hsp65_references.py MIN_MINOR OUT_DIR` builds the consensus references (MIN_MINOR 1 = every
polymorphic site, 2 = a base needs >= 2 genomes). Job `main3_reference_alignment/submit_hsp65_iupac_experiment.sbatch`
writes `output/mlsa/main3_reference_alignment_hsp65_iupac/`: A = baseline (the amplicon run), B = `--iupac` on the 72
amplicons, C/D = `--iupac` on the consensus (MIN_MINOR 1/2), E = old scoring on the same 7 consensus references.

**Result** (127 hsp65 reads with a whole-genome species): A, B, C, D and E are identical, 110 correct, 1 wrong
(the mixed culture), 16 without a call (15 divergent, 1 no hit), and no call flips. The no-calls are 77-85 % identity
reads (Rothia contaminant, failed traces), not ambiguity-code problems. The margin to the second species changed
for 101 of the 201 `ok` reads (see below), but no call. 16S and main2's distances were not tested.
The consensus has 11 ambiguous sites for kansasii (mostly from the 3 atypical-hsp65 genomes) and at most 3 for the others.

Independent of `--iupac`, `_COMP` (used to put a reverse-strand hit back in the read's orientation) now has the full
IUPAC complement; before, it only knew N, so R, Y, K, M, S and W were not complemented. Same fix in
sanger-microsynth `sanger_ms/refalign.py`; the `--iupac` code exists only in mlsa-kansasii. Known remaining
inconsistency: `p_distance_matrix` in `mlsa/align.py` (main2) counts IUPAC letters as mismatches, whereas
`pairwise_p_distance` skips them.

### The margin to the second-best species

`margin` = differences of the best hit of *another species* minus differences of the closest hit, over the read
(`ReadResult.runner_up`, `status` in `mlsa/refalign.py`). A call is `ok` only if identity is >= 98 % and the margin is
>= 2 (`--min-margin`). It is a count of differing bases, not a probability: the threshold is a heuristic and a margin
of 3 is 3 bases, not 3 sigma.

What drives it (all 201 `ok` reads of the Sanger set; truth subset = the 127 reads with a whole-genome species):
- **The 3 Korean *M. kansasii* genomes (GCF_002705785.1, -825.1, -865.1), not IUPAC handling.** Their hsp65 is
  persicum-like (11 SNPs from the other kansasii), so for a persicum read the "second species" is *kansasii*, via these
  three genomes, 3 differences away. That is why persicum reads sit at margin 2-3 (mean 2.85; 5 reads at margin 2, 24 at 3) while
  all other species have 6-16. Run F (consensus references built without the 3 genomes, `submit_hsp65_iupac_excluded.sbatch`)
  moves the persicum margin to 7.1-7.4 (second species is then gastri) and the 2 `ambiguous` reads become `ok`;
  the 111 `ok` / 110 correct on the truth subset are unchanged. The margin grows only because the nearest kansasii-like
  sequence is removed, so the larger number is not better evidence: a *Korean-type* kansasii isolate (kansasii by ANI
  and at every other locus, LIT.md) would still be called persicum, with no hint of the problem. Keeping the 3 genomes
  gives the honest margin.
- **No Sanger read matches these genomes.** `atypical_hsp65` is false for all reads in every run, so there is no
  Korean-like isolate in this set to test the case; "Korean" refers to the genomes, not to reads. Runs A-E (and the
  published amplicon run) include the 3 genomes; only F leaves them out.
- **IUPAC consensus (C/D) changes margins only by averaging genomes:** kansasii +0.6 (7.4 -> 8.0), pseudokansasii -0.6,
  attenuatum -2, innocens -3, persicum 0 (still 2.85, because the consensus kansasii codes the Korean alleles at 11
  sites and persicum reads match them). That is a property of 7 consensus references vs 72 genomes, not a gain in
  evidence. It also hides which genome is closest, and the atypical flag cannot be raised from a species consensus.

**Are the results more trustworthy?** No. Calls, correctness (110 of 111, one mixed culture) and the number of reads
with margin <= 3 on the truth subset (11) are the same in A-E. The only run with fewer weak margins (F, 0) removes
the information that makes them weak. The margin that matters for trust is the one against the nearest *real genome*
(run A); it says persicum vs kansasii in hsp65 alone is a 3-base call. A tie-breaking second locus (gyrA) is the fix,
not ambiguity handling.

## References used by main3 (and in the lab's alignment PDFs)

main3 aligns each read to one whole genome per species
(`data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_gtdb_representatives/`). The lab's own
hsp65 alignment reports (`data/imm/2026500522_hsp65_Tb-11_alignment.pdf` and
`..._alignment_comparison.pdf`, made with Apache FOP, not by this pipeline) use only three hsp65
reference sequences (373 bp). Comparison:

| Species | main3 reference (genome) | In the lab PDFs (hsp65 only) | Same strain? |
|---|---|---|---|
| kansasii | GCF_000157895.3, ATCC 12478 | ATCC 12478 | yes |
| gastri | GCF_002102175.1, DSM 43505 | CIP_104530 | yes: CIP 104530 = DSM 43505 = ATCC 15754 = CCUG 20995 = JCM 12407, the type strain |
| ostraviense | GCF_002705925.1, 241/15 (type strain, DSM 110538) | FDAARGOS_1613 | no, a different strain |
| persicum | GCF_002086675.1, AFPC-000227 | not used | n/a |
| pseudokansasii | GCF_900566075.1, MK142 | not used | n/a |
| innocens | GCF_900566055.1, MK13 | not used | n/a |
| attenuatum | GCF_900566085.1, MK41 | not used | n/a |

The PDFs therefore only separate kansasii, gastri and ostraviense. They cannot give persicum,
pseudokansasii, innocens or attenuatum, so a read from one of those is forced onto the nearest of
the three. For 2026500522 the isolate differs from gastri at 7 positions, ostraviense 13 and
kansasii 15 (98.12% identity to gastri at best), so the PDF call is not comparable with main3's.

### Strains behind the two extra references

- **CIP_104530** (*M. gastri*): the Institut Pasteur collection number of the species' type strain.
  The same strain is DSM 43505 (the genome GCF_002102175.1 used by main3), ATCC 15754, CCUG 20995
  and JCM 12407. It was isolated from gastric lavage. *M. gastri* shares an identical 16S rDNA
  sequence with *M. kansasii*; hsp65 and ITS tell them apart. So the PDF's gastri and main3's
  gastri are the same organism, just named by a different collection.
- **FDAARGOS_1613** (*M. ostraviense*): a complete, circular chromosome (GenBank CP089224,
  assembly GCA_021183725.1 / ASM2118372v1, 5,114 protein genes and 53 RNA genes, entered 2022)
  from the FDA-ARGOS database of quality-controlled reference genomes for diagnostic use. It is
  not the type strain (241/15, from sputum, Karviná, Czech Republic; Jagielski et al. 2019,
  Front Microbiol). Submitter and exact origin were not checked; the NCBI assembly page did not
  load. It is not among the genomes in this repo (no GCF accession is linked).
  Being a different ostraviense strain, it may differ from 241/15 at a few hsp65 positions, so
  the ostraviense distances in the PDF and in main3 are not directly comparable.

Sources: [Mycobacterium gastri, type strain designations](https://lpsn.dsmz.de/species/mycobacterium-gastri),
[DSM 43505](https://www.dsmz.de/collection/catalogue/details/culture/DSM-43505),
[KEGG Mycobacterium ostraviense FDAARGOS_1613](https://www.kegg.jp/kegg-bin/show_organism?org=mot),
[M. ostraviense DSM 110538](https://www.dsmz.de/collection/catalogue/details/culture/DSM-110538).

## How outputs are used downstream

| Step | Produces | Used by |
|---|---|---|
| `scripts/03_link_sanger_kansasii.sh` | `.ab1` symlinks in `data/sanger/seq_kansasii/` | main2 |
| `scripts/02_link_gtdb_kansasii_complex.sh` | genome `.fna`/`.gff` symlinks in `data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii/` | main1, main2 (only if main1's output is missing), genome_identity_check |
| main1 | `alignments/<locus>.raw.fasta` | main2 (hsp65 and 16S reference sequences) |
| main1 | `single_locus_pair_separation.tsv` | main2 (`recommended_additional_loci.tsv`) |
| main1 | `combo_results.tsv`, `winning_combo/`, `figures/`, `SUMMARY.txt` | final results |
| main2 | `isolate_classification.tsv`, `representative_reads.tsv`, `recommended_additional_loci.tsv`, `figures/`, `SUMMARY.txt` | final results |
| `*_excluded` variants | the same, in their own output folders | main2_excluded reads main1_excluded; otherwise final results |
| genome_identity_check | `skani_ani.tsv`, `skani_long.tsv`, `species_ani_summary.tsv` | no code; the ANI numbers are cited by hand in LIT.md |

The genome_identity_check results inform decisions about which reference
genomes to trust or exclude. No pipeline step reads them.

## How main2 classifies isolates (`isolate_classification.tsv`)

Implemented in `main2_sanger_differentiation/run_sanger_differentiation.py`.
The table has one row per isolate (`probennummer`) per `locus` (`hsp65`,
`16S`, `hsp65+16S`). `TNR` is the TNR of the read that was picked (for
hsp65+16S, both TNRs, hsp65 first, if they differ).

1. **One sequence per isolate and locus.** Reads under all of an isolate's
   TNR columns are pooled. The locus comes from the filename
   (`mlsa.sanger_io.classify_locus`): hsp65 by gene name or TB-11/TB-12(w)
   primer in their historical spellings, 16S by the MBAK-14, Mbak259r or
   Mbak264r primer. Every `.ab1` read is quality-trimmed (Mott, reads
   shorter than 100 bp are dropped) and flipped if needed to match the
   reference strand. Each read is compared to every reference by pairwise
   local alignment. The longest read (ties going to higher mean quality) that
   is within `max_dist` (step 4) of some reference is kept. If no read is,
   the longest read is kept anyway and ends up NA. So a failed first read
   (contaminant, mixed culture) doesn't hide a good repeat or re-extraction.
   Forward and reverse reads are not merged into a consensus. The chosen read
   per isolate and locus, with its screening distance and how many reads
   were available, is in `representative_reads.tsv`.
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
   `margin = second_dist - nearest_dist`, `tolerance` is the larger of the
   two species' tolerances (the same rule as main1's barcoding-gap check),
   and **`unambiguous = margin > tolerance`**. So the second-closest species
   has to be further away than the references of either species are from
   each other.
   **Outside the complex → NA.** `max_dist` is the largest distance between
   any two reference sequences at that locus (hsp65 0.054, 16S 0.008,
   hsp65+16S 0.018 with the full reference set). If `nearest_dist >
   max_dist`, the read is further from every reference than any two members
   of the complex are from each other (a contaminant, another genus' groEL,
   or a bad read). `nearest_species`, `second_species` and `unambiguous` are
   then NA, and the isolate gets no locus recommendation.
5. **hsp65+16S.** Only isolates with both reads get this row. The two
   single-locus alignments are joined end to end, keeping only references
   that have both loci, and steps 3–4 are repeated.

Isolates still ambiguous with hsp65+16S get a suggested extra locus in
`recommended_additional_loci.tsv`: the main1 candidate locus that best
separates their nearest and second-closest species.

Things to keep in mind when reading the table:

- The tolerance comes from the reference genomes of the two closest
  species. A diverse species such as kansasii (hsp65 tolerance 0.025) makes
  every isolate whose nearest or second-closest species is kansasii
  ambiguous, even when it sits right next to a reference.
- The tolerance comes from full-length reference sequences, while an
  isolate's distances only cover its read. For short reads (16S reads cover
  about a third of the gene) the tolerance can be too strict.
- Gaps are skipped, so a short read that overlaps little of the alignment is
  judged on only a few positions.
- Reference 16S/ITS: some GTDB assemblies contain contaminant contigs with
  their own 16S. `mlsa.loci.extract_16s` skips partial copies and keeps the
  copy closest to the type-strain 16S (see LIT.md, "Extraction artifact").
  Partial genes are also skipped for the gene-symbol loci (rpoB, gyrA, gyrB,
  recA, secA1, tuf).

The TNR isolates will also be Illumina-sequenced and speciated with the
`kansasii` branch of
[immensekansasii](https://gitlab.uzh.ch/appliedmicrobiologyresearch/immense.git)
(`../immensekansasii`); this repo answers the same species-complex question
computationally ahead of that NGS data being available.

## Pipeline runs (immensekansasii)

Start every run of the immensekansasii (IMMense) pipeline in its own
subdirectory of `/shares/sander.imm.uzh/MM/kansasii/runs/`, not in `output/`.
Nextflow writes its `work/` directory into the directory the run is started
from. `work/` is large and only temporary, so it must not end up in
`output/`, which gets downloaded to local `kansasii_C`. Once a run has
finished and been checked, copy its end results to
`/shares/sander.imm.uzh/MM/kansasii/output/<run_name>/` and delete `work/`.

```bash
mkdir -p /shares/sander.imm.uzh/MM/kansasii/runs/<run_name>
cd /shares/sander.imm.uzh/MM/kansasii/runs/<run_name>
bash /shares/sander.imm.uzh/MM/kansasii/repos/immensekansasii/run_IMMENSE.sh -j <job_name> -t <input_type> -r <run_name> -i <input_dir>

# after the run: collect results (without work/) in output/
rsync -a --exclude work --exclude .nextflow /shares/sander.imm.uzh/MM/kansasii/runs/<run_name>/ /shares/sander.imm.uzh/MM/kansasii/output/<run_name>/
```

## Windows checkout

- `.gitattributes` forces LF line endings, so a Windows checkout (even with
  `core.autocrlf=true`) keeps scripts runnable. Recommended: `git config core.autocrlf false`
  and `git config core.longpaths true`.
- Scripts default to the cluster root `/shares/sander.imm.uzh/MM/kansasii`; set the
  `KANSASII_ROOT` environment variable to point elsewhere (e.g. a mapped drive).
- Nextflow, Singularity and sbatch steps only run on the cluster (or WSL).


## Provenance

The sbatch files write `MANIFEST.txt` (commit, clean/dirty flag, script, conda env hash, inputs) into the output
directory when they finish; see `kansasii-docs/docs/PROVENANCE.md`. Project index: `kansasii-docs/README.md`.
