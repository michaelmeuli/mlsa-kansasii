#!/usr/bin/env bash
#SBATCH --job-name=mlsa-hsp65-public
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mlsa/logs/%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

python "$REPO_ROOT/sanger_vs_wgs/run_insilico_public.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-docs/bin/write_manifest.sh \
  "/shares/sander.imm.uzh/MM/kansasii/output/mlsa/sanger_vs_wgs_public" \
  /shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii submit_insilico_public.sbatch \
  /shares/sander.imm.uzh/MM/kansasii/data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_hsp65_amplicons
