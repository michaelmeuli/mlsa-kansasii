#!/usr/bin/env bash
#SBATCH --job-name=mlsa-genome-identity-check
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.log
# Singularity fails on this node (no /etc/resolv.conf).
#SBATCH --exclude=u24-cva0000-129

set -euo pipefail

module load apptainer

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO_ROOT/genome_identity_check/run_genome_identity_check.py"
