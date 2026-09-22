#!/bin/bash
#
# Pull the Singularity containers used by mlsa-kansasii into the shared
# container path also used by immensekansasii. blast and iqtree are already
# provided there by immensekansasii's own pull script and are reused as-is;
# this script only adds mafft, and skips any image that already exists so it
# never re-downloads multi-GB images unnecessarily.
#
# Usage:
#   conda activate env_mlsa
#   cd /shares/sander.imm.uzh/MM/kansasii/repos/mlsa-kansasii
#   bash Singularity/pull_singularity_img.sh /shares/sander.imm.uzh/software/pipelines/IMMense/IMMense_dependencies/containers

set -euo pipefail

module load apptainer

INSTALL_PATH=${1:?"Usage: $0 <install_path>"}
mkdir -p "$INSTALL_PATH"
echo ""
echo "Installing containers at:            ${INSTALL_PATH}"
echo ""

pull_if_missing() {
  local img_name="$1"
  local docker_uri="$2"
  local target="${INSTALL_PATH}/${img_name}"
  if [ -d "$target" ]; then
    echo "Already present, skipping: ${img_name}"
  else
    echo "Building: ${img_name} <- ${docker_uri}"
    singularity build --sandbox "$target" "docker://${docker_uri}"
  fi
}

pull_if_missing "quay.io-biocontainers-mafft-7.525--h031d066_1.img" "quay.io/biocontainers/mafft:7.525--h031d066_1"

# Reused from immensekansasii's shared container set (built there already):
#   quay.io-biocontainers-blast-2.17.0--h66d330f_0.img
#   quay.io-biocontainers-iqtree-3.1.3--h8471819_0.img
for reused in \
  "quay.io-biocontainers-blast-2.17.0--h66d330f_0.img" \
  "quay.io-biocontainers-iqtree-3.1.3--h8471819_0.img"
do
  if [ -d "${INSTALL_PATH}/${reused}" ]; then
    echo "Confirmed shared container present: ${reused}"
  else
    echo "WARNING: expected shared container not found: ${reused} (run immensekansasii's Singularity/pull_singularity_img.sh first)" >&2
  fi
done
