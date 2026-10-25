# Mkan329-035 (TNR 2024500472): species disagreement and contamination

Checked 2026-10-05. TNR taken from `output/screening_map_results.csv`.

## Summary

**Mkan329-035 is most likely a mixed culture of *M. kansasii* (about 30%) and *M. ostraviense* (about 70%). It can be
excluded from the kansasii analyses.** This rests on sequence data only; single-colony re-isolation has not been done
(see "Alternative explanations" below).

- **Both species calls are consistent with the mixture.** The two methods read different parts of the same mixed DNA.
  - WGS/GTDB says *kansasii* because the genome contains a complete *M. kansasii* genome (58x depth, 97.6% covered).
  - Sanger hsp65 / MLSA says *ostraviense* because most of the template in the PCR was *M. ostraviense* (139x depth, 96.7% covered).
    The best hsp65 read matches ostraviense at 99.47%.
  - Neither method needs to have made an error. If the sample is not a single species, there is no single correct answer.
- **Four independent lines of evidence agree with a mixture:** CheckM (113.37% contamination), a 13.3 Mb assembly
  (about twice one genome), mixed Sanger hsp65 traces, and the read mapping below.
- **Why exclude it:** a mixed sample has no valid species label, genome, resistance profile or MIC. The assembly joins
  contigs of two species, so SNPs, gene content (including resistance genes) and the tree position of this isolate are not
  those of one organism. Keeping it would add false variants to the association analyses and a wrongly placed tip to the MLSA tree,
  and its MIC reflects both populations. Excluding it costs one isolate, and the reason for excluding it
  (contamination) does not depend on the phenotype.

## The disagreement

| Source | Call |
|---|---|
| WGS / GTDB (`species`) | kansasii (ANI 99.66, AF 0.953, ref GCF_000157895.3) |
| WGS CheckM contamination | **113.37%**, `contamination_flag = contaminated` |
| sanger_mlsa hsp65 (Sanger) | ostraviense (dist 0.0053), gastri second (0.0214), unambiguous |
| sanger_mlsa 16S (Sanger) | kansasii and ostraviense **tied** (0.00113 each), not unambiguous |
| sanger_mlsa hsp65+16S | ostraviense |
| sanger_refalign reference alignment (hsp65, 65kDa 5ul read) | ostraviense 99.47% (2 diffs / 374 bp), gastri 97.86% |

`species_mlsa1`, `species_mlsa2` and `species_ref` in screening_map_results.csv all say ostraviense.

**Conclusion:** the contamination flag, the Sanger traces and the read mapping together indicate that the sample
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
  and the fewest mixed positions, and it is the representative read picked by sanger_mlsa and sanger_refalign. The
  ostraviense call depends on that read. The two later reads are too noisy to call.
- The 10ul read has no Q>=20 bases at all.

## Can mixed bases alone exclude an isolate from Sanger identification?

Not by a fixed threshold. `mlsa/sanger_io.py` now has `mixed_peak_fraction` (share of base calls in the
Mott-trimmed read whose second peak is at least 25% of the primary one), and sanger_mlsa writes `mixed_fraction` (chosen read)
and `max_mixed_fraction` (worst usable read) to `representative_reads.tsv`. No flag is applied. For Mkan329-035 the
best hsp65 read has 2.1% (second read 5.3%, 16S 0.1%), because the two noisiest reads are trimmed below 100 bp and never reach sanger_mlsa.
Across all 252 hsp65 reads in `data/sanger/seq_kansasii` the median is 1.6%, 25% of reads are at or above 5% and 10% are at or above 18%,
so the isolate looks ordinary by this measure. The CheckM `contamination_flag` (and the read mapping above) are
the exclusion signal.

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
7. **WGS.** CheckM contamination or mapping to both reference species
   (the route that already caught this isolate).


## Read mapping to both references (2026-10-05)

The Illumina reads (`Mkan329-035_r1/r2`) were mapped with bwa-mem to one combined reference of
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
dominant component, and GTDB's kansasii call picks up the minor one. So neither call needs to be wrong, and the
data are consistent with a mixed culture of *M. kansasii* and *M. ostraviense*.

Caveats: this uses one reference per species. The mix ratio is an estimate from depth, not a measured fraction of cells.

### Alternative explanations

- **A new species or hybrid** is unlikely but is not excluded by an explicit test. A single novel genome would have one
  depth across its length and would not cover two reference genomes at about 97% breadth with <1 mismatch per read, since the
  species differ by several percent in ANI. This is an argument from the numbers above, not a measurement (no
  ANI or core-genome analysis of separated contigs was done).
- **Cross-contamination** (library prep, sequencing run, or the Sanger PCR) would also give two organisms in the data
  without a mixed culture in the strict sense. Nothing here separates the two. Checking neighbouring isolates and
  the no-template control (see Recommendation 3) is the only handle.
- **Decisive test:** single colonies from a re-streak, each sequenced (WGS or hsp65 in both directions). Two species recovered
  from separate colonies would settle it. Until then, "mixed culture" is the best-supported interpretation, not proven.

## Recommendation

1. **Exclude Mkan329-035** from the MLSA tree and the MIC association analyses. Reason to give if
   an exclusion list is kept: "probable mixed culture, *M. kansasii* + *M. ostraviense* (CheckM 113%, read mapping 31%/69%)".
   Do not use either species call as the isolate's label. In `isolate_classification`, mark it as contaminated and do not use the ostraviense call.
2. **Optional rescue and confirmation:** re-streak, pick single colonies, and re-sequence (WGS and hsp65 in both directions).
   This is also the decisive test for mixed culture vs. novel species. A split-by-mapping reassembly is possible,
   but a single-colony genome is more reliable.
3. **Check neighbours:** other isolates from the same plate or run (check `contamination_flag` in
   `screening_map_results.csv`) may have the same issue.
