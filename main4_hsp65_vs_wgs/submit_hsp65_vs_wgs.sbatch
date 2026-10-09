#!/usr/bin/env bash
#SBATCH --job-name=mlsa-hsp65-vs-wgs
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO_ROOT/main4_hsp65_vs_wgs/run_hsp65_vs_wgs.py"
