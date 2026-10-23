#!/usr/bin/env bash
#SBATCH --job-name=mlsa-hsp65-vs-wgs
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO_ROOT/sanger_vs_wgs/run_sanger_vs_wgs.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-docs/bin/write_manifest.sh \
  "/shares/sander.imm.uzh/MM/kansasii/output/mlsa/sanger_vs_wgs" \
  /shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii submit_sanger_vs_wgs.sbatch /shares/sander.imm.uzh/MM/kansasii/data/imm/screening_map_link.csv
