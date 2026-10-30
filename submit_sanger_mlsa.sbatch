#!/usr/bin/env bash
# 08: how far do the Sanger hsp65 + 16S reads resolve species -> output/mlsa-kansasii/results/sanger_mlsa/ (after 04)
#SBATCH --job-name=mlsa-sanger-differentiation
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa-kansasii/logs/%x-%j.log

set -euo pipefail
module load apptainer
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/mlsa-kansasii
OUTPUT=$OUT_MLSA/results
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO/scripts/08_sanger_mlsa.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/sanger_mlsa" "$REPO" submit_sanger_mlsa.sbatch "$SCREENING_MAP_LINK" "$SANGER_LINK_DIR" "$OUTPUT/locus_discovery"
