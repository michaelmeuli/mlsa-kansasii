# CLAUDE.md

See README.md for the question, layout, setup, run order and how main2
classifies isolates. This file only lists what is easy to get wrong.

- **Outputs are on the shares, not in the repo:**
  `/shares/sander.imm.uzh/MM/kansasii/output/mlsa/<program>/`. A repo-local
  `results/` folder is obsolete. Don't read from it or write to it.
- **Environment:** the system `python3` has no pandas or Biopython. Use
  `conda activate env_mlsa`. Real runs go through the `submit_*.sbatch` SLURM
  jobs (`module load apptainer`, mafft/iqtree via Singularity). Don't run the
  full pipelines on the login node. Job logs (`*-<jobid>.log`) are written to the
  directory `sbatch` was run from.
- **Paths are hardcoded:** data paths in `mlsa/__init__.py`, which resolve to
  `/shares/sander.imm.uzh/MM/kansasii/data`; `RESULTS`/`MAIN1_RESULTS` in
  each `run_*.py`; `REPO_ROOT` in the sbatch files. Update all of them if
  anything moves.
- **Inputs are symlinks.** Run `scripts/link_sanger_kansasii.sh` (Sanger
  `.ab1`) and `scripts/link_gtdb_kansasii_complex.sh` (GTDB genomes) again
  after the source data changes. The isolate list is
  `data/imm/screening_map_link.csv`.
- **Isolates are keyed by `PROBENNUMMER`.** Reads under all TNR columns
  (`SCREENING_MAP_TNR_COLUMNS`) are pooled per isolate. Output tables with a
  `tnr` column come from older code and are out of date.
- **`*_excluded` variants** rerun main1/main2 without the 3 kansasii genomes
  GCF_002705785.1, GCF_002705825.1 and GCF_002705865.1. They call the base
  scripts with different parameters, so change the base script and not a
  copy. LIT.md explains why those genomes are real *M. kansasii*, so
  excluding them is a sensitivity check, not a correction.
- **Order:** main1 before main2, and `main1_*_excluded` before
  `main2_*_excluded`. main2 reuses main1's reference alignments and
  pair-separation table.
- **Species list and colors** are defined once in `mlsa/__init__.py`
  (`SPECIES`, `SPECIES_COLORS`). Reuse them; don't redefine them.
