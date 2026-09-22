#!/usr/bin/env bash
#SBATCH --job-name=mlsa-locus-discovery
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail

module load apptainer

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO_ROOT/main1_locus_discovery/run_locus_discovery.py"
