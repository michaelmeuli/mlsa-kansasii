#!/usr/bin/env bash
# 10: Sanger hsp65 calls vs whole-genome species -> output/mlsa/sanger_vs_wgs/ (after 08, 09 and the results table)
#SBATCH --job-name=mlsa-hsp65-vs-wgs
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/mlsa-kansasii
OUTPUT=$OUTPUT_ROOT/mlsa
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO/scripts/10_sanger_vs_wgs.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/sanger_vs_wgs" "$REPO" submit_sanger_vs_wgs.sbatch "$SCREENING_MAP_RESULTS" "$GTDB_REPRESENTATIVES" "$ASSEMBLIES" "$OUTPUT/sanger_mlsa" "$OUTPUT/sanger_mlsa_excluded" "$OUTPUT/sanger_refalign"
