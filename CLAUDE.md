# CLAUDE.md

See README.md for the question, layout, setup, run order and how sanger_mlsa
classifies isolates. This file only lists what is easy to get wrong.

- **Outputs are on the shares, not in the repo:**
  `/shares/sander.imm.uzh/MM/kansasii/output/mlsa/<program>/`. A repo-local
  `results/` folder is obsolete. Don't read from it or write to it.
- **Environment:** the system `python3` has no pandas or Biopython. Use
  `conda activate env_mlsa`. Real runs go through the `submit_*.sbatch` SLURM
  jobs (`module load apptainer`, mafft/iqtree via Singularity). Don't run the
  full pipelines on the login node. Job logs (`*-<jobid>.log`) are written to
  `output/mlsa/logs/` (the directory must exist before `sbatch`).
- **Paths are hardcoded:** data paths in `mlsa/__init__.py`, which resolve to
  `/shares/sander.imm.uzh/MM/kansasii/data`; `RESULTS`/`LOCUS_DISCOVERY_RESULTS` in
  each `run_*.py`; `REPO_ROOT` in the sbatch files. Update all of them if
  anything moves.
- **Inputs are symlinks.** Run `scripts/03_link_sanger_kansasii.sh` (Sanger
  `.ab1`) and `scripts/02_link_gtdb_kansasii_complex.sh` (GTDB genomes) again
  after the source data changes. The isolate list is
  `data/imm/screening_map_link.csv`.
- **Isolates are keyed by `PROBENNUMMER`.** Reads under all TNR columns
  (`SCREENING_MAP_TNR_COLUMNS`) are pooled per isolate. Output tables with a
  lowercase `tnr` column come from older code and are out of date. sanger_mlsa's
  uppercase `TNR` column is only the TNR of the read that was picked.
  `screening_map_results.csv` (immensekansasii) flags reads whose TNR differs
  from the isolate's `TNR` in `TNR_MLSA`.
- **`*_excluded` variants** rerun locus_discovery/sanger_mlsa without the 3 kansasii genomes
  GCF_002705785.1, GCF_002705825.1 and GCF_002705865.1. They call the base
  scripts with different parameters, so change the base script and not a
  copy. LIT.md explains why those genomes are real *M. kansasii*, so
  excluding them is a sensitivity check, not a correction.
- **Order:** locus_discovery before sanger_mlsa, and `locus_discovery_excluded` before
  `sanger_mlsa_excluded`. sanger_mlsa reuses locus_discovery's reference alignments and
  pair-separation table.
- **Species list and colors** are defined once in `mlsa/__init__.py`
  (`SPECIES`, `SPECIES_COLORS`). Reuse them; don't redefine them.
- **immensekansasii pipeline runs** go in their own subdirectory of
  `/shares/sander.imm.uzh/MM/kansasii/runs/<run_name>/`, never in `output/`. The Nextflow `work/`
  directory is large and temporary; only the end results are copied to
  `/shares/sander.imm.uzh/MM/kansasii/output/<run_name>/` (see README.md).
- **`mlsa/refalign.py`** (sanger_refalign) is copied to sanger-microsynth as `sanger_ms/refalign.py`; keep them in sync.
