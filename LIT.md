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
- The excluded run (hsp65 alone: 117/119 isolates unambiguous, versus 22/119
  with all genomes) is an over-estimate. The three genomes are the only ones
  that inflate kansasii's within-species tolerance at hsp65, so dropping them
  assumes atypical-hsp65 kansasii do not occur.
- Treat the full-reference runs as primary and the `*_excluded` runs as a
  sensitivity analysis.
- hsp65 called 11 of the 119 isolates persicum. Some could be atypical
  kansasii like these; only gyrA would tell them apart. Not checked.

## Extraction artifact: GCF_900565995.1

The spurious 0.45 within-kansasii tolerance at 16S comes from
GCF_900565995.1, whose extracted 16S is truncated (1,128 bp instead of about
1,520 bp) and whose ITS is truncated (161 bp instead of about 285 bp). This is
an extraction artifact, not biology. GCF_002086895.1 also has an ITS 8% from
the type strain, which has not been investigated.

## Reproducing

ANI outputs are in `output/mlsa/genome_identity_check/` (`skani_ani.tsv`,
`skani_long.tsv`, `genomes.txt`). skani ships in the GTDB-Tk container:

```bash
module load apptainer
singularity exec --bind /shares \
  /shares/sander.imm.uzh/software/pipelines/IMMense/IMMense_dependencies/containers/quay.io-biocontainers-gtdbtk-2.7.2--pyhdfd78af_1.img \
  skani triangle -l genomes.txt -E -t 8 -s 70 -o skani_ani.tsv
```

Run it as an sbatch job rather than on the login node. Node `u24-cva0000-129`
has no `/etc/resolv.conf` and Singularity fails there, so use
`--exclude=u24-cva0000-129`.
