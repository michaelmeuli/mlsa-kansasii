#!/usr/bin/env bash
# 04: smallest locus set that separates the 7 species -> output/mlsa-kansasii/results/locus_discovery/
#SBATCH --job-name=mlsa-locus-discovery
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

python "$REPO/scripts/04_locus_discovery.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/locus_discovery" "$REPO" submit_locus_discovery.sbatch "$GTDB_COMPLEX_GENOMES"
