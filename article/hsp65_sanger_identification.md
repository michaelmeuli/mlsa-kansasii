# One Sanger read is enough: hsp65 identifies the *Mycobacterium kansasii* complex, and the lab can automate it

*Draft, updated 2026-10-09. Numbers come from `output/mlsa/main1`–`main4`, `output/mlsa/main3_reference_alignment_hsp65_amplicons` and `output/screening_map_results.csv`.*

## Summary

The *M. kansasii* complex contains at least seven species (kansasii, persicum, pseudokansasii, innocens, attenuatum, ostraviense, gastri) with different clinical and resistance profiles. Our routine identifies them by Sanger sequencing of a 441 bp hsp65 (groEL1) amplicon (primers TB-11/TB-12w) and, sometimes, 16S. We asked whether this is enough, using 126 clinical isolates whose species is known from whole-genome sequencing (GTDB-Tk), 92 of which also have routine Sanger hsp65 reads.

- **hsp65 alone gave the genome species for 125 of 126 isolates** when the amplicon was extracted from the assembly (the ceiling for a perfect read). The one miss is a mixed culture (kansasii + ostraviense, CheckM contamination flagged), where hsp65 correctly reports the ostraviense component.
- **Real Sanger reads: 82 of 83 calls were correct (99%)**, aligned to the hsp65 amplicons of 72 reference genomes at a 98% identity threshold. 9 of 92 isolates (10%) got no call, all because of the read, not the marker.
- **16S cannot do this job**: it does not separate kansasii from gastri, persicum, pseudokansasii or ostraviense (12/21 species pairs resolved, against 17/21 for hsp65 and 21/21 for gyrA among 72 reference genomes).
- The calls are now built into `sanger-microsynth`, which fetches Microsynth results by mail, pairs forward and reverse reads per sample, applies fixed action rules and sends a PDF report.

## Methods in brief

| Program | What it does |
|---|---|
| main1 | Extracts 9 candidate loci from 72 GTDB r232 genomes of the 7 species and tests which pairs of species each locus separates (DNA-barcoding gap). |
| main2 | MLSA-style: aligns each isolate's hsp65/16S read with main1's reference alignment, estimates within-species diversity from all references, calls `unambiguous` only if the second-closest species is further than that diversity. |
| main3 | Plain alignment of each read to references; call = fewest differences; `ok` needs ≥98% identity and ≥2 differences fewer than the closest other species. References: one genome per species, or (hsp65) the amplicons of all 72 genomes. |
| main4 | Benchmark of the above against whole-genome species; also extracts hsp65 in silico from each isolate's own assembly. |

Truth set: 126 isolates with a GTDB-Tk species in the complex (85 kansasii, 24 persicum, 14 pseudokansasii, 2 attenuatum, 1 innocens; no gastri or ostraviense). 92 have a representative Sanger hsp65 read.

## Results

### 1. At reference level, hsp65 is almost enough, and the exception is one lineage

| Locus | Species pairs separated (of 21) |
|---|---|
| gyrA, gyrB, rpoB, recA, tuf | 21 |
| hsp65 | 17 (21 without three Korean genomes) |
| ITS | 16 |
| 16S | 12 |
| secA1 | 6 |

The four failing hsp65 pairs (kansasii vs gastri, ostraviense, persicum, pseudokansasii) all come from one cause: three near-clonal Korean isolates (GCF_002705785.1/825.1/865.1) that are *M. kansasii* by ANI (98.5–99.3% to other kansasii, ≤94% to any other species) and identical to the type strain at gyrA, rpoB, 16S and others, but carry a persicum-like hsp65 (11 SNPs; the "atypical hsp65 types" of Jagielski et al. 2019). Without them, hsp65 separates all 21 pairs. They are real, so the practical rule is: **hsp65 can misplace a kansasii that carries an atypical hsp65 as persicum.** None of our 85 kansasii isolates behaves like this (all within 0–1 differences of the ATCC 12478 hsp65), and none of the 252 hsp65 reads in the lab archive matches those genomes.

### 2. A perfect hsp65 read identifies 125/126 isolates

For every isolate the in-silico amplicon had status `ok` (≥98% identity, margin ≥4 differences to the closest other species; median margin 8). Distance to the correct reference: 0 differences in 119/126, 1 in 3, 2 in 4.

### 3. Real Sanger reads: accurate when they work

| Method (92 isolates with a hsp65 read) | Called | Correct | No call |
|---|---|---|---|
| main3, 7 representative genomes, ≥99% identity | 79 | 78 | 13 |
| main3, hsp65 amplicons of 72 genomes, ≥98% identity | 83 | 82 | 9 |
| main2, nearest species | 82 | 81 | 10 |

The single error in each is the mixed culture (Mkan329-035). The 4 extra calls at 98% are noisy reads at 98.3–99.0% identity (mixed-peak fractions up to 10%) whose closest reference was always the correct species. The 9 remaining non-calls are 8 failed reads (78–82% identity to everything, one of them only 173 bp, one at 93.5%) and 1 read with no alignment. So about 10% of routine reads need to be repeated, and nothing in the remainder is wrong.

Which references are used hardly matters for the call: the 7 reference strains alone give the same result (83 ok, 82 correct) as all 72 genomes. Adding one Korean genome to the 7 turns 2 persicum isolates ambiguous (81 ok), because the Korean amplicon lies within 1 difference of them. The 72-genome set is kept because it also flags the atypical lineage.

### 4. Why main2 looks much worse than it is

Counting main2's `unambiguous` flag, hsp65 resolves only 8 of 127 isolates (3 of 42 with 16S added), because the three Korean genomes inflate kansasii's within-species tolerance to 2.5%; removing them gives 110/110 unambiguous. The *nearest species* was right in both runs (81/82). The flag is conservative by design: it measures what hsp65 cannot guarantee in the worst reference lineage, not how it performs on the isolates we see. Plain alignment (main3) agrees with main2 on all 78 shared calls. It has no such tolerance, which is why the atypical-hsp65 reference group matters.

## What this means for the lab

1. **Keep hsp65 as the routine test and stop adding 16S** for this question. 16S ties between several complex species and adds no information.
2. **Add gyrA only for defined cases:** the read matches the atypical-hsp65 kansasii genomes, two species are not clearly separated (margin <2), or an isolate needs certainty for reporting. gyrA separates all 21 species pairs and was the recommended extra locus for 30 of 31 isolates still ambiguous under hsp65+16S.
3. **Treat mixed cultures as a QC flag, not a species:** the one discordant isolate was a contaminated culture. A high share of secondary peaks should trigger a repeat from a pure colony.

## Automation with sanger-microsynth

`sanger-microsynth` reads Microsynth result mails (IMAP, read-only), unzips `.ab1`/FASTA, Mott-trims, BLASTs against NCBI, aligns each read to the references, and emails a PDF with a sample table and one alignment PDF per read. Implemented rules (`sanger_ms/call.py`):

| Situation | Action |
|---|---|
| ≥98% identity, ≥2 differences ahead of the closest other species, mixed peaks <10%, read ≥200 bp | report species |
| closest = atypical-hsp65 kansasii genome, or species not separated | confirm with gyrA |
| no reference ≥98% identical, or read <200 bp | repeat sequencing |
| ≥10% of base calls with a secondary peak | possible mixed culture, repeat from pure colony |
| 16S read | not species-informative for this complex |

**Forward and reverse reads** (TB-11 / TB-12w) are paired by sample number and locus. Each read is called on its own and the sample call requires all usable reads to agree; otherwise the sample is "repeat, reads disagree". In the lab archive (183 forward, 64 reverse reads) 51 samples have both reads, 37 have a call on both, and none disagree, so no consensus sequence is built. Each read covers about 85% of the amplicon.

hsp65 reads are aligned to the amplicons of all 72 genomes (within-species diversity, atypical lineage explicit); the alignment PDF shows only the 7 reference strains (ATCC 12478, AFPC-000227, MK142, MK13, MK41, 241/15, DSM 43505), plus the Korean genome if it is the closest hit. On the 2026-09-30 batch the pipeline reported kansasii (2025341808) and innocens (2025500758), asked to repeat two hsp65 reads of 84 and 94 bp, and showed that two of five 16S reads were not mycobacteria of the complex.

## Limits

- 126 isolates from one lab, skewed to kansasii: no gastri or ostraviense, one innocens, two attenuatum. For the rare species the claim rests on reference genomes only (1–3 each).
- Whole-genome species is the truth; the 3 contaminated isolates are the only discordance.
- Atypical-hsp65 kansasii would be flagged only if it resembles the three Korean genomes; none was seen in the lab data.
- 34 of the 126 genome-typed isolates have no hsp65 Sanger read yet.
- The rules were validated retrospectively on the archive; the first prospective batches should be checked by hand.

## Reproduction

mlsa-kansasii: `main3_reference_alignment/submit_hsp65_amplicon_alignment.sbatch` (references from `scripts/make_hsp65_references.py`) and `main4_hsp65_vs_wgs/`; outputs in `output/mlsa/`. sanger-microsynth: `sanger_ms/call.py`, `main2_analyse_report/run_analyse_report.py`.
