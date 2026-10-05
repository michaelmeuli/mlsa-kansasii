#!/usr/bin/env bash
#SBATCH --job-name=mlsa-reference-alignment
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

# every linked Sanger read (about 3 s per read)
python "$REPO_ROOT/main3_reference_alignment/run_reference_alignment.py" \
    /shares/sander.imm.uzh/MM/kansasii/data/sanger/seq_kansasii
