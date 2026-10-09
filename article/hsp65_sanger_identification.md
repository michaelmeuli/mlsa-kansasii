# One Sanger read is enough: hsp65 identifies the *Mycobacterium kansasii* complex, and the lab can automate it

*Draft, 2026-10-09. Numbers come from `output/mlsa/main1`–`main4` and `output/screening_map_results.csv`.*

## Summary

The *M. kansasii* complex contains at least seven species (kansasii, persicum, pseudokansasii, innocens, attenuatum, ostraviense, gastri) with different clinical and resistance profiles. Our routine identifies them by Sanger sequencing of a 441 bp hsp65 (groEL1) amplicon (primers TB-11/TB-12w) and, sometimes, 16S. We asked whether this is enough, using 126 clinical isolates whose species is known from whole-genome sequencing (GTDB-Tk), 92 of which also have routine Sanger hsp65 reads.

- **hsp65 alone gave the genome species for 125 of 126 isolates** when the amplicon was extracted from the assembly (the ceiling for a perfect read). The one miss is a mixed culture (kansasii + ostraviense, CheckM contamination flagged), where hsp65 correctly reports the ostraviense component.
- **Real Sanger reads: 78 of 79 calls were correct (99%)**, with plain alignment to one reference per species. 13 of 92 isolates (14%) got no call, almost all because of the read, not the marker.
- **16S cannot do this job**: it does not separate kansasii from gastri, persicum, pseudokansasii or ostraviense (11/21 species pairs resolved, against 17/21 for hsp65 and 21/21 for gyrA among 72 reference genomes).
- A plain alignment of each read to seven reference genomes (main3) takes about 3 s and is already built into `sanger-microsynth`, which fetches Microsynth results by mail and sends a PDF back. With a few changes (below) the lab can drop manual interpretation for routine reads.

## Methods in brief

| Program | What it does |
|---|---|
| main1 | Extracts 9 candidate loci from 72 GTDB r232 genomes of the 7 species and tests which pairs of species each locus separates (DNA-barcoding gap). |
| main2 | MLSA-style: aligns each isolate's hsp65/16S read with main1's reference alignment, estimates within-species diversity from all references, calls `unambiguous` only if the second-closest species is further than that diversity. |
| main3 | Plain alignment of each read to one genome per species; call = fewest differences; `ok` needs ≥99% identity and ≥2 differences fewer than the runner-up. No locus extraction needed. |
| main4 | Benchmark of the above against WGS species; also extracts hsp65 in silico from each isolate's own assembly. |

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

The four failing hsp65 pairs (kansasii vs gastri, ostraviense, persicum, pseudokansasii) all come from one cause: three near-clonal Korean isolates (GCF_002705785.1/825.1/865.1) that are *M. kansasii* by ANI (98.5–99.3% to other kansasii, ≤94% to any other species) and identical to the type strain at gyrA, rpoB, 16S and others, but carry a persicum-like hsp65 (11 SNPs; the "atypical hsp65 types" of Jagielski et al. 2019). Without them, hsp65 separates all 21 pairs. They are real, so the practical rule is: **hsp65 can misplace a kansasii that carries an atypical hsp65 as persicum.** None of our 85 kansasii isolates behaves like this (all within 0–1 differences of the ATCC 12478 hsp65).

### 2. A perfect hsp65 read identifies 125/126 isolates

For every isolate the in-silico amplicon had status `ok` (≥99% identity, margin ≥4 differences to the runner-up; median margin 8). Distance to the correct reference: 0 differences in 119/126, 1 in 3, 2 in 4. Species-level margins were 4–8 differences for persicum and kansasii, the closest pair.

### 3. Real Sanger reads: accurate when they work

| Method (92 isolates with a hsp65 read) | Called | Correct | No call |
|---|---|---|---|
| main2, nearest species | 82 | 81 | 10 |
| main3, status `ok` | 79 | 78 | 13 |

The single error in each is the mixed culture (Mkan329-035). The 13 main3 non-calls are:

- 8 reads at 78–82% identity to everything (and one 173 bp read, one at 93.5%): failed or non-hsp65 sequences, not a species-limit.
- 4 reads at 98.3–99.0% identity (noisy reads, mixed-peak fractions up to 10%) whose closest reference was nevertheless the correct species. A 98% identity threshold would have called all four correctly (82/83).
- 1 read with no alignment.

So roughly **10% of routine reads need to be repeated**, and nothing in the remainder is wrong.

### 4. Why main2 looks much worse than it is

Counting main2's `unambiguous` flag, hsp65 resolves only 8 of 127 isolates (3 of 42 with 16S added), because the three Korean genomes inflate kansasii's within-species tolerance to 2.5%; removing them gives 110/110 unambiguous. The *nearest species* was right in both runs (81/82). The flag is conservative by design: it measures what hsp65 cannot guarantee in the worst reference lineage, not how it performs on the isolates we see. Plain alignment (main3) has no such tolerance and agrees with main2 on all 78 shared calls, but it is also blind to this limitation, which is why the atypical-hsp65 caveat above matters.

## What this means for the lab

1. **Keep hsp65 as the routine test and stop adding 16S** for this question. 16S ties between several complex species and adds no information.
2. **Add gyrA only for defined cases:** hsp65 calls persicum (or a kansasii/persicum tie), a mixed-peak read, or an isolate needing certainty for reporting. gyrA separates all 21 pairs (it was the recommended extra locus for 30 of 31 isolates still ambiguous under hsp65+16S).
3. **Treat mixed cultures as a QC flag, not a species:** the one discordance was a contaminated culture. A high share of secondary peaks should trigger a pure-culture repeat.

## Automating the analysis with sanger-microsynth

`sanger-microsynth` already: reads Microsynth result mails (IMAP, read-only), unzips `.ab1`/FASTA, Mott-trims, BLASTs against NCBI, aligns each read to the seven reference genomes (the same code as main3), and emails a PDF report with a per-read alignment. Proposed rules for hands-off use:

| Situation | Automatic output |
|---|---|
| hsp65 read ≥98% identity, ≥2 differences fewer than runner-up, mixed fraction <10% | Species reported, no human review |
| Closest reference persicum, or margin <2 | "Confirm with gyrA", queued for review |
| Identity <98% to every reference, or read <200 bp | "Repeat sequencing" (about 10% of reads) |
| Mixed fraction ≥10% | "Possible mixed culture", repeat from pure colony |
| 16S read | "Not species-informative for this complex" |

Changes needed: lower the `ok` threshold from 99% to 98% (validated here on 4 reads), add forward/reverse consensus (not done yet), add several references per species instead of one (intra-species diversity is not covered), and include the lab's hsp65 species list in the report.

## Limits

- 126 isolates from one lab, species-skewed: no gastri or ostraviense isolates, one innocens and two attenuatum. The claim for the rare species rests on reference genomes only.
- Reference sets are small for the rare species (1–3 genomes).
- WGS species is the truth; 3 isolates are flagged as contaminated and are the only discordance.
- Atypical-hsp65 kansasii (persicum-like) would be misidentified; none was seen here, but nothing in hsp65 alone would reveal it.
- 34 of the 126 WGS isolates have no hsp65 Sanger read yet.

## Reproduction

`main4_hsp65_vs_wgs/` (sbatch job `submit_hsp65_vs_wgs.sbatch`), outputs in `output/mlsa/main4_hsp65_vs_wgs/` (`sanger_vs_wgs.tsv`, `insilico_vs_wgs.tsv`, `SUMMARY.txt`).
