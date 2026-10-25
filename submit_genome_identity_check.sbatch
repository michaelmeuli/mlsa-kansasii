#!/usr/bin/env bash
# 05: all-vs-all ANI (skani) of the 72 genomes -> output/mlsa/genome_identity_check/
#SBATCH --job-name=mlsa-genome-identity-check
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --exclude=u24-cva0000-129
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail
module load apptainer
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/mlsa-kansasii
OUTPUT=$OUTPUT_ROOT/mlsa
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO/scripts/05_genome_identity_check.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/genome_identity_check" "$REPO" submit_genome_identity_check.sbatch "$GTDB_COMPLEX_GENOMES"
