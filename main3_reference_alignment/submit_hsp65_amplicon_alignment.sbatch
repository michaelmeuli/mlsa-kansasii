#!/usr/bin/env bash
#SBATCH --job-name=mlsa-refalign-hsp65
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
ROOT="/shares/sander.imm.uzh/MM/kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

# hsp65 reads against the 72 per-genome hsp65 amplicons (make with scripts/make_hsp65_references.py)
python "$REPO_ROOT/main3_reference_alignment/run_reference_alignment.py" \
    --ref-dir "$ROOT/data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_hsp65_amplicons" \
    --out "$ROOT/output/mlsa/main3_reference_alignment_hsp65_amplicons" \
    "$ROOT/data/sanger/seq_kansasii"
