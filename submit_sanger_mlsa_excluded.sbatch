#!/usr/bin/env bash
# 08 against the reference set without the atypical genomes -> output/mlsa/sanger_mlsa_excluded/ (after 04 --exclude-atypical)
#SBATCH --job-name=mlsa-sanger-differentiation-excluded
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail
module load apptainer
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/mlsa-kansasii
OUTPUT=$OUTPUT_ROOT/mlsa
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO/scripts/08_sanger_mlsa.py" --exclude-atypical

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/sanger_mlsa_excluded" "$REPO" submit_sanger_mlsa_excluded.sbatch "$SCREENING_MAP_LINK" "$SANGER_LINK_DIR" "$OUTPUT/locus_discovery_excluded"
