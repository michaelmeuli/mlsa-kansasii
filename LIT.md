# Literature and genome evidence notes

## Are GCF_002705785.1, GCF_002705825.1 and GCF_002705865.1 real *M. kansasii*?

**Conclusion: yes.** They are *M. kansasii* sensu stricto (the former subtype I).
What is atypical is their hsp65 gene, not the genomes. Excluding them is not
justified as a misidentification fix, and it inflates hsp65's apparent
performance (see "Consequences for the analyses" below).

### The genomes

| Accession | Strain | Origin | Submitter | BioSample |
|---|---|---|---|---|
| GCF_002705785.1 | K4 | Ulsan, South Korea, sputum (SMC, 2008) | Forschungszentrum Borstel | SAMN07193056 |
| GCF_002705865.1 | K14 | Seoul, South Korea, sputum (SMC, 2011) | Forschungszentrum Borstel | SAMN07193077 |
| GCF_002705825.1 | K19 | Seoul, South Korea, sputum (SMC, 2012) | Forschungszentrum Borstel | SAMN07193078 |

All three are BioProject PRJNA224116, draft scaffolds (308-375 contigs, checkM2
completeness >99.9%, contamination <0.03%). GTDB r232 places all three in
*M. kansasii*, with ATCC 12478 (GCF_000157895.3) as the species representative.
They have identical scaffold statistics (a 6,518,727 bp scaffold, 268 spanned
gaps), which suggests reference-guided scaffolding.

### Genome-wide identity (skani 0.3.2, all 72 genomes of the 7 species)

- ANI to the other *M. kansasii* genomes: 98.5-99.3%.
- ANI to the type strain ATCC 12478: 98.96-99.19% (alignment fraction 93-94%).
- ANI to every other species: at most 93.95% (innocens 93.5-93.95%, persicum
  93.0-93.9%, pseudokansasii 92.1-93.0%, ostraviense 91.5-92.7%, gastri
  91.2-91.8%, attenuatum 89.7-90.8%).
- The other *M. kansasii* genomes span 98.6-100% ANI among themselves, so the
  three sit at the low end of the kansasii range, not outside it. The usual
  species cutoff is 95-96%.
- The three are near-clonal: K14/K19 ANI 100.0%, K4 vs K14/K19 99.9%. For
  purposes of independent evidence they count as one lineage.

### Marker genes (from `main1_locus_discovery` alignments)

| Locus | Distance to type strain | Nearest species |
|---|---|---|
| gyrA, gyrB, rpoB, recA, tuf, secA1 | 0.0000-0.0004 | kansasii (other species 2-12% away) |
| 16S, ITS | 0.0000 | kansasii |
| **hsp65** | **0.0249** | **persicum (0.0023)** |

At hsp65 the amplicon (441 bp, TB-11/TB-12w primers) lies in the same GroEL
gene as in the type strain, at the same relative position, with no N bases.
At 11 sites inside it the three genomes carry the *M. persicum* allele, where
39 other kansasii genomes carry the kansasii allele. hsp65 alone therefore
calls them persicum. Only the 441 bp amplicon was examined, so the extent of
the persicum-like block is unknown. Gene transfer between complex members is
a plausible explanation (see Tagini 2021 below), but this was not tested.

### Literature (via PubMed)

- Jagielski T, et al. *Genomic insights into the Mycobacterium kansasii
  complex: an update.* Front Microbiol 2019;10:2918.
  [doi:10.3389/fmicb.2019.02918](https://doi.org/10.3389/fmicb.2019.02918)
  (PMC6974680). Sequenced K4, K14 and K19. Its whole-genome grouping puts all
  three in subtype I. Its PCR-RFLP/sequencing typing had called K14 and K19
  "atypical type II (IIb)" and K4 "atypical type I (Ib)", based on a unique
  hsp65 sequence most similar to type II (classification of Iwamoto and Saito
  2005, as cited by Jagielski; not retrieved separately). The paper also
  concludes that none of the common taxonomic markers (16S, hsp65, rpoB, tuf,
  ITS) distinguishes all complex species, and that whole-genome approaches
  should be used.
- Tagini F, et al. *Phylogenomics reveal that Mycobacterium kansasii subtypes
  are species-level lineages.* Int J Syst Evol Microbiol 2019;69:1696-1704.
  [doi:10.1099/ijsem.0.003378](https://doi.org/10.1099/ijsem.0.003378).
  ANI between subtypes 88.4-94.2%, below the 95-96% species cutoff; described
  *M. pseudokansasii*, *M. innocens* and *M. attenuatum*; *M. persicum* shares
  99.77% ANI with subtype 2.
- Tagini F, et al. *Pathogenic determinants of the Mycobacterium kansasii
  complex: an unsuspected role for distributive conjugal transfer.*
  Microorganisms 2021;9:348.
  [doi:10.3390/microorganisms9020348](https://doi.org/10.3390/microorganisms9020348).
  Reports gene transfer between complex members.

### Discrepancy not resolved

Jagielski et al. describe K14 and K19 as having a type II ITS sequence. The
ITS extracted here (the 16S-23S gap from the GFF rRNA features) is identical
to the kansasii type strain for all three, and 7.6% from the nearest persicum
genome. This may reflect a different region or sequevar definition. It does
not affect the ANI conclusion.

### Consequences for the analyses

- gyrA is unaffected and still separates all 7 species with or without the
  three genomes.
- The excluded run (hsp65 alone: 123/127 isolates unambiguous, versus 9/127
  with all genomes) is an over-estimate. The three genomes are the only ones
  that inflate kansasii's within-species tolerance at hsp65, so dropping them
  assumes atypical-hsp65 kansasii do not occur.
- Treat the full-reference runs as primary and the `*_excluded` runs as a
  sensitivity analysis.
- hsp65 calls 11 of the 127 isolates persicum, each with kansasii as the
  second-closest species. With all genomes, all 11 are ambiguous: main2's
  tolerance is the larger of the two species' diversity, which here is
  kansasii's 0.025 from these three genomes. Some of the 11 could be
  atypical kansasii like these; only gyrA would tell them apart. Not checked.

## Extraction artifact: contaminant 16S copies (fixed)

The spurious within-kansasii 16S tolerance (0.42 to 0.45) came from
GCF_900565995.1 (strain MK7). Its 16S was not truncated. It came from a
different organism. The assembly has 13 contigs below 55% GC (about 42 kb in
a 66% GC genome), and NCBI annotated two extra 16S copies on them:
NZ_UPHJ01000162.1 (37% GC, partial, 1,128 bp) and NZ_UPHJ01000138.1 (52% GC,
1,810 bp). The old `extract_16s` took the first annotated copy, which was the
37% GC contaminant, and the 161 bp ITS came from the same wrong operon. The
genome's own 16S (NZ_UPHJ01000094.1, 1,531 bp) is typical *M. kansasii*.
persicum GCF_900566005.1 has the same contamination pattern (likely the same
sequencing batch). It was unaffected only because its real copy is listed
first.

Two more genomes only have partial 16S copies at contig edges:
pseudokansasii GCF_019751035.1 (95 bp, which inflated the pseudokansasii 16S
tolerance to 0.41) and persicum GCF_902825395.1 (495 bp).

`mlsa.loci.extract_16s`/`extract_its` now skip partial copies and keep the
full-length copy sharing the most 8-mers with the type-strain 16S
(GCF_000157895.3). A copy sharing under 50% is rejected. The real copies
share 0.98 to 1.0, the contaminants 0.05. GCF_019751035.1 and GCF_902825395.1
therefore have no 16S or ITS. The other 68 genomes are unchanged.

GCF_002086895.1 also has an ITS 8% from the type strain and a persicum-type
16S. See "What 16S can and cannot resolve" below.

## What 16S can and cannot resolve

### Reference genomes

GTDB r232 gives 72 genomes for the 7 species. Not every locus could be
extracted from every genome:

| Species | Genomes | with hsp65 | with 16S (= with both) |
|---|---|---|---|
| kansasii | 42 | 42 | 42 |
| persicum | 12 | 12 | 11 (GCF_902825395.1: only a partial 16S) |
| pseudokansasii | 6 | 6 | 5 (GCF_019751035.1: only a partial 16S) |
| attenuatum | 4 | 4 | 4 |
| innocens | 3 | 3 | 3 |
| ostraviense | 3 | 3 | 3 |
| gastri | 2 | 2 | 2 |
| **total** | **72** | **72** | **70** |

The `*_excluded` runs drop 3 kansasii genomes (39 kansasii, 69 genomes with
hsp65, 67 with 16S).

Over the 1,511 16S positions all 70 genomes share, there are only 7 distinct
sequences (H1 to H7):

| Variant | Genomes | Differences from H1 |
|---|---|---|
| H1 | kansasii (41), ostraviense (3), gastri (2) | - |
| H2 | persicum (11), kansasii GCF_002086895.1 | 8 |
| H3 | pseudokansasii (5) | 8 (1 from H2) |
| H4 | innocens (3) | 2 |
| H5 to H7 | attenuatum (2 + 1 + 1) | 9 to 10 (1 to 3 from H2/H3) |

- kansasii, ostraviense and gastri have identical 16S sequences. 16S cannot
  separate them at all.
- innocens differs from that group by 2 positions. persicum, pseudokansasii
  and attenuatum differ from each other by only 1 to 3.
- 16S mainly splits the complex into two groups:
  kansasii/ostraviense/gastri/innocens and persicum/pseudokansasii/attenuatum,
  8 to 12 positions apart.
- kansasii GCF_002086895.1 carries the persicum 16S. It is the only source of
  within-kansasii 16S variation (tolerance 0.0053; 0.0 without it). Its ITS is
  also 8% from the type strain. Genome-wide it is a real *M. kansasii*: ANI
  99.15% to the type strain (AF 93%) and 98.6-99.7% to the other kansasii,
  versus at most 93.3% to persicum and 93.5% to any other species
  (`genome_identity_check`). So it is a kansasii with a persicum-like rRNA
  operon, like the persicum-like hsp65 of K4/K14/K19. Gene transfer or an
  assembly artifact are both possible; neither was tested.
- By the barcoding-gap criterion (`single_locus_pair_separation.tsv`), 16S
  separates 12 of 21 species pairs. None of the pairs involving kansasii are
  separated.

### Isolates (main2, Sanger reads)

The Sanger reads cover only part of the gene. 59 of 65 isolates are exactly
the same distance from two species: kansasii/ostraviense (41),
kansasii/persicum (10) or pseudokansasii/attenuatum (8). Three are
unambiguous: Mkan329-124 (kansasii), Mkan329-171 (pseudokansasii) and
Mkan329-056 (attenuatum). These calls are fragile. kansasii and ostraviense
have identical reference 16S over the shared region, so Mkan329-124's
separation must come from positions at the ends where not all reference
sequences have data (gaps are skipped). Mkan329-056 is 2.7% from its
nearest reference, which points to a noisy read.

### Is hsp65+16S better than hsp65 alone? No, apart from one isolate.

Species pairs separated in the reference genomes (barcoding gap), and
unambiguous isolates among the 42 that have both reads:

| | hsp65 | 16S | hsp65+16S |
|---|---|---|---|
| Reference pairs separated, all genomes | 17/21 | 12/21 | **15/21** |
| Reference pairs separated, `_excluded` | 21/21 | 12/21 | **20/21** |
| Unambiguous isolates, all genomes | 3/42 | 0/42 | **4/42** |
| Unambiguous isolates, `_excluded` | 42/42 | 0/42 | **34/42** |

- The combined locus calls the same nearest species as hsp65 for all 42
  isolates.
- With all genomes, adding 16S resolves one more isolate: Mkan329-092
  (persicum, second kansasii). With hsp65 alone its margin (0.0087) is below
  kansasii's 0.025 tolerance. In the combined locus the 16S positions shrink
  kansasii's tolerance more than the margin. That is a side effect of the
  dilution described below, not extra information from 16S: its 16S alone
  is ambiguous. Relative to hsp65, adding 16S loses the attenuatum-kansasii
  and innocens-kansasii reference pairs.
- In the `_excluded` run, it loses 8 isolates that hsp65 alone resolves
  (Mkan329-006, -020, -028, -039, -044, -049, -050, -146), and the
  gastri-kansasii reference pair.

Why adding 16S makes things worse:

1. 16S is about 1,500 of the 1,993 aligned positions and has almost no
   variation, so it shrinks every distance and margin by roughly 4x.
2. The kansasii tolerance is set by the most divergent pair of kansasii
   genomes, and GCF_002086895.1's persicum 16S makes that pair more
   divergent. Combined kansasii tolerance: 0.0096 (0.0056 without
   GCF_002086895.1). In the `_excluded` run: 0.0046 (0.0005 without it).

So 16S adds nothing to hsp65 for species identification in this complex.
hsp65 alone (and gyrA, main1's winning locus) is the better basis. A
combined analysis only makes sense after GCF_002086895.1 is resolved, and
even then 16S can at best confirm the broad two-group split.

## Reproducing

The ANI results come from `genome_identity_check/`:

```bash
cd genome_identity_check
sbatch submit_genome_identity_check.sbatch
```

`run_genome_identity_check.py` lists the same 72 GTDB r232 genomes main1
uses (`mlsa.loci.discover_genomes`) and runs `skani triangle -E -s 70`
(skani 0.3.2 from the GTDB-Tk container, `mlsa.align.run_skani_triangle`).
It writes to `output/mlsa/genome_identity_check/`:

- `genomes.txt`: the `.fna` paths.
- `skani_ani.tsv`: raw skani output, one row per genome pair.
- `skani_long.tsv`: each pair in both directions with accession, species,
  ANI and AF (the smaller of the two aligned fractions).
- `species_ani_summary.tsv`: ANI range per species pair.
- `skani_version.txt`

The sbatch file excludes node `u24-cva0000-129`, which has no
`/etc/resolv.conf`, so Singularity fails there.

The original run (job 6226012, 2026-09-21) is in
`output-old/mlsa/genome_identity_check/`. The folder was moved there
together with the rest of `output/`. The rerun on 2026-09-26 gives an
identical `skani_long.tsv`.
