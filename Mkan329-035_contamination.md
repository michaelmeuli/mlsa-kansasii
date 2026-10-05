# Mkan329-035 (TNR 2024500472): species disagreement and contamination

Checked 2026-10-05. TNR taken from `output/screening_map_results.csv`.

## Summary

**Mkan329-035 is a mixed culture of *M. kansasii* (about 30%) and *M. ostraviense* (about 70%). It can be excluded
from the kansasii analyses.**

- **Both species calls are correct.** The two methods read different parts of the same mixed DNA.
  - WGS/GTDB says *kansasii* because the genome contains a complete *M. kansasii* genome (58x depth, 97.6% covered).
  - Sanger hsp65 / MLSA says *ostraviense* because most of the template in the PCR was *M. ostraviense* (139x depth, 96.7% covered).
    The best hsp65 read matches ostraviense at 99.47%.
  - Neither method made an error. The sample is not a single species, so there is no single correct answer.
- **The mixture is independently confirmed** by CheckM (113.37% contamination), a 13.3 Mb assembly (about twice one genome),
  mixed Sanger hsp65 traces, and the read mapping below.
- **Why exclude it:** a mixed sample has no valid species label, genome, resistance profile or MIC. The assembly joins
  contigs of two species, so SNPs, gene content (including resistance genes) and the tree position of this isolate are not
  those of one organism. Keeping it would add false variants to the GWAS and a wrongly placed tip to the MLSA tree,
  and its MIC reflects both populations. Excluding it costs one isolate, and the reason for excluding it
  (contamination) does not depend on the phenotype.

## The disagreement

| Source | Call |
|---|---|
| WGS / GTDB (`species`) | kansasii (ANI 99.66, AF 0.953, ref GCF_000157895.3) |
| WGS CheckM contamination | **113.37%**, `contamination_flag = contaminated` |
| main2 hsp65 (Sanger) | ostraviense (dist 0.0053), gastri second (0.0214), unambiguous |
| main2 16S (Sanger) | kansasii and ostraviense **tied** (0.00113 each), not unambiguous |
| main2 hsp65+16S | ostraviense |
| main3 reference alignment (hsp65, 65kDa 5ul read) | ostraviense 99.47% (2 diffs / 374 bp), gastri 97.86% |

`species_mlsa1`, `species_mlsa2` and `species_ref` in screening_map_results.csv all say ostraviense.

**Conclusion:** the contamination flag, the Sanger traces and the read mapping together show that the culture
contains both species. The two calls therefore do not conflict (see Summary). The 16S read cannot break the
tie between the two species, so the Sanger ostraviense call rests on hsp65, which is from the dominant component.

## Sanger trace check (secondary-peak fraction)

Method: for each base in the Q>=20 window of the six `.ab1` files with TNR 2024500472, compare the
second-highest and highest channel at the base-call position (`PLOC2`). "Mixed" means second/primary >= 0.25.
Script: `mix.py` (scratchpad, Biopython, env `biopython`).

| Read | HQ bp | mean Q | mixed positions | median sec/pri |
|---|---|---|---|---|
| 2024-09-16 LC9_Mbak14 (16S) | 893 | 57 | 6 (0.7%) | 0.002 |
| 2024-09-23 hsp65_11 | 386 | 35 | 31 (8.0%) | 0.069 |
| 2025-01-20 65kDa 5ul | 383 | 41 | 13 (3.4%) | 0.044 |
| 2025-01-20 65kDa 10ul | none | | | |
| 2025-01-21 hsp65 FastStart | 381 | 21 | 112 (29.4%) | 0.167 |
| 2025-01-29 65kDal | 370 | 21 | 148 (40.0%) | 0.195 |

Reading this:
- The 16S read is clean, so the 16S template was a single sequence (or the mixture is invisible at 16S).
  That locus does not separate kansasii from ostraviense here, so a clean trace says nothing about species.
- hsp65 is mixed in every read. The fraction runs from 3.4% to 40% between reads of the same culture,
  which suggests either a variable ratio between PCRs or poor reads. The 20.01 5ul read has the best Q (41)
  and the fewest mixed positions, and it is the representative read picked by main2 and main3. The
  ostraviense call depends on that read. The two later reads are too noisy to call.
- The 10ul read has no Q>=20 bases at all.

## Can mixed bases alone exclude an isolate from Sanger identification?

Not by a fixed threshold. `mlsa/sanger_io.py` now has `mixed_peak_fraction` (share of base calls in the
Mott-trimmed read whose second peak is at least 25% of the primary one), and main2 writes `mixed_fraction` (chosen read)
and `max_mixed_fraction` (worst usable read) to `representative_reads.tsv`. No flag is applied. For Mkan329-035 the
best hsp65 read has 2.1% (second read 5.3%, 16S 0.1%), because the two noisiest reads are trimmed below 100 bp and never reach main2.
Across all 252 hsp65 reads in `data/sanger/seq_kansasii` the median is 1.6%, 25% of reads are at or above 5% and 10% are at or above 18%,
so the isolate looks ordinary by this measure. The CheckM `contamination_flag` (and the read mapping above) are
the exclusion signal. The GWAS table already excludes NR 35 (`kansasii-gwas/scripts/01_build_table.py`).

## Ways to detect contamination in Sanger data

1. **Secondary peaks in the trace.** A contaminant at 10-50% gives a second peak at the same positions
   along the read. Count mixed bases in the HQ region (as above) or use Tracy `decompose`
   or `sangerseqR::makeBaseCalls`. Noise at the read ends is not contamination.
2. **Forward vs reverse reads.** If both directions are clean but give different species, the template is mixed.
3. **Indel-driven mixing.** A clean trace that turns into a double signal that never recovers means
   two templates of different length. Tracy `decompose` can separate the two sequences.
4. **Several loci.** Compare hsp65, rpoB, 16S, gyrA on the same DNA. Disagreement across the
   species border (kansasii vs ostraviense) with clean traces is a mixed culture or a real mosaic.
5. **BLAST margin.** A small gap between the top two species means the locus does not resolve them
   (as the 16S tie here).
6. **Neighbours and controls.** Look for an identical sequence in a nearby well or run, and check the no-template control.
7. **WGS.** CheckM contamination, Kraken2 read classification or mapping to both reference species
   (the route that already caught this isolate).


## Read mapping to both references (2026-10-05)

Kraken2 could not be used: the shared `kraken_standard_db/hash.k2d` is truncated (68.64 GB on disk,
71.38 GB required by its header), so both attempts (jobs 7031513, 7032204) aborted. Instead the
Illumina reads (`Mkan329-035_r1/r2`) were mapped with bwa-mem to one combined reference of
*M. kansasii* GCF_000157895.3 and *M. ostraviense* GCF_002705925.1 (job 7033025, files in
`output/mapping_Mkan329-035/`). Only primary reads with MAPQ >= 20 are counted, i.e. reads that tell the two genomes apart.

| | kansasii | ostraviense |
|---|---|---|
| unique reads | 2,581,067 (31%) | 5,705,895 (69%) |
| mean depth | 58x | 139x |
| breadth of coverage | 97.6% | 96.7% |
| mean mismatches per read (NM) | 0.88 | 0.80 |

96.1% of reads mapped, 93.8% properly paired, so little of the library comes from other species.

**Interpretation:** the sample contains two near-complete genomes, each covered over about 97% of its
length with almost no mismatches. A single isolate would give few reads on the wrong reference, or reads with many more
mismatches, because these species differ by several percent in ANI. The depths suggest roughly 30% kansasii and
70% ostraviense, which matches the 13.3 Mb assembly (about twice the size of one genome), the
CheckM contamination of 113% and the mixed hsp65 traces. The Sanger hsp65 call (ostraviense) matches the
dominant component, and GTDB's kansasii call picks up the minor one. So neither call is wrong, and this
is a mixed culture of *M. kansasii* and *M. ostraviense*.

Caveats: this uses one reference per species. The mix ratio is an estimate from depth, not a measured fraction of cells.

## Recommendation

1. **Exclude Mkan329-035** from the kansasii GWAS, the MLSA tree and the MIC association analyses. Reason to give if
   an exclusion list is kept: "mixed culture, *M. kansasii* + *M. ostraviense* (CheckM 113%, read mapping 31%/69%)".
   Do not use either species call as the isolate's label. In `isolate_classification`, mark it as contaminated and do not use the ostraviense call.
2. **Optional rescue:** re-streak, pick single colonies, and re-sequence (WGS and hsp65 in both directions).
   A split-by-mapping reassembly is possible, but a single-colony genome is more reliable.
3. **Check neighbours:** other isolates from the same plate or run (check `contamination_flag` in
   `screening_map_results.csv`) may have the same issue.
