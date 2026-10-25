#!/usr/bin/env bash
# 09 with the 72 per-genome hsp65 amplicons as references -> output/mlsa/sanger_refalign_hsp65_amplicons/ (needs 06)
#SBATCH --job-name=mlsa-refalign-hsp65
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/mlsa-kansasii
OUTPUT=$OUTPUT_ROOT/mlsa
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO/scripts/09_sanger_refalign.py" --ref-dir "$HSP65_AMPLICONS" --out "$OUTPUT/sanger_refalign_hsp65_amplicons" "$SANGER_LINK_DIR"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/sanger_refalign_hsp65_amplicons" "$REPO" submit_hsp65_amplicon_alignment.sbatch "$SANGER_LINK_DIR" "$HSP65_AMPLICONS"
