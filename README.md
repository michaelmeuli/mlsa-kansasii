# mlsa-kansasii

Can the *Mycobacterium kansasii* species complex (kansasii, persicum,
pseudokansasii, innocens, attenuatum, ostraviense, gastri) be differentiated
by Sanger sequencing / MLSA with as few loci as possible, and how far does
the lab's existing hsp65 + 16S Sanger data already get?

## Layout

- `scripts/link_sanger_kansasii.sh` — symlinks the `.ab1` Sanger traces of
  the TNR isolates in `data/imm/screening_map.csv` into
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
```

Outputs land in `results/main1_locus_discovery/` and
`results/main2_sanger_differentiation/` (tables + `figures/`, not tracked in
git — see `.gitignore`).

The TNR isolates will also be Illumina-sequenced and speciated with the
`kansasii` branch of
[immensekansasii](https://gitlab.uzh.ch/appliedmicrobiologyresearch/immense.git)
(`../immensekansasii`); this repo answers the same species-complex question
computationally ahead of that NGS data being available.
