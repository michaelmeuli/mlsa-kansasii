#!/usr/bin/env bash
#SBATCH --job-name=mlsa-locus-discovery-excluded
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail

module load apptainer

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO_ROOT/main1_locus_discovery_excluded/run_locus_discovery_excluded.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-docs/bin/write_manifest.sh \
  "/shares/sander.imm.uzh/MM/kansasii/output/mlsa/main1_locus_discovery_excluded" \
  /shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii submit_locus_discovery_excluded.sbatch /shares/sander.imm.uzh/MM/kansasii/data/imm/screening_map_link.csv
