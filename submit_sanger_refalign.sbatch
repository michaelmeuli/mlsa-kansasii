#!/usr/bin/env bash
# 09: plain alignment of every linked Sanger read to one genome per species -> output/mlsa-kansasii/results/sanger_refalign/
#SBATCH --job-name=mlsa-reference-alignment
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa-kansasii/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/mlsa-kansasii
OUTPUT=$OUT_MLSA/results
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

# every linked Sanger read (about 3 s per read)
python "$REPO/scripts/09_sanger_refalign.py" "$SANGER_LINK_DIR"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/sanger_refalign" "$REPO" submit_sanger_refalign.sbatch "$SANGER_LINK_DIR" "$GTDB_REPRESENTATIVES"
