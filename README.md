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
- `main3_reference_alignment/` — species ID of Sanger reads (`.ab1`/FASTA) by plain alignment to one
  reference genome per species (`kansasii_complex_gtdb_representatives/`), no MLSA. One row per read in
  `reference_alignment.tsv` plus a PDF per read (difference matrix + alignment with all 7 references).
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

To download all results to a Windows machine (PowerShell):

```powershell
New-Item -ItemType Directory -Path "$env:USERPROFILE\kansasii_C\downloads\" -Force
scp -r mimeul@cluster.s3it.uzh.ch:/shares/sander.imm.uzh/MM/kansasii/output/* "$env:USERPROFILE\kansasii_C\downloads\"
```

## Plain alignment to references vs MLSA

`python main3_reference_alignment/run_reference_alignment.py <reads or dirs>` (about 3 s per read).
Call = closest of the 7 representatives by differences over the read: `ok` (&ge;99% identity, &ge;2
differences fewer than the next), `ambiguous`, `divergent` (&lt;99%, e.g. new lineage, introgressed
locus, or bad read). Compared with main1/main2:

- Plain alignment needs no locus extraction, no annotation and no tolerance estimate, works on any read,
  and takes seconds. It also never sees the 72-genome diversity, so the 3 Korean *M. kansasii* genomes
  (LIT.md) cannot disturb it. But that is blindness, not robustness: an isolate carrying their
  persicum-like hsp65 would be called persicum, with the representatives giving no hint of the problem.
- main2 (MLSA-style) measures within-species diversity from all references and so says when a call is
  not supported (`unambiguous`), but is hsp65/16S-limited and needs main1's alignments.
- Neither fixes a single-locus limit: hsp65 and 16S alone cannot separate every species (16S reads
  tie between several references). More loci (gyrA, see main2 recommendations) or genome data are needed.

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
| `scripts/link_sanger_kansasii.sh` | `.ab1` symlinks in `data/sanger/seq_kansasii/` | main2 |
| `scripts/link_gtdb_kansasii_complex.sh` | genome `.fna`/`.gff` symlinks in `data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii/` | main1, main2 (only if main1's output is missing), genome_identity_check |
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
