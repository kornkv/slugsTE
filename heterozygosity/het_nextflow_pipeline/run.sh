#!/usr/bin/env bash
# Convenience launcher for the Nextflow pipeline.
# Usage:
#   ./run.sh                  # one sample at a time (default profile)
#   ./run.sh -profile parallel    # multiple samples in parallel (big instance)
#   ./run.sh -resume          # resume after a failure / parameter tweak
set -euo pipefail

cd "$(dirname "$0")"

# Java is required by Nextflow (and by the Nextflow self-installer below).
# Install OpenJDK 17 if missing — must happen BEFORE the nextflow install step.
if ! command -v java >/dev/null 2>&1; then
  echo "[run] Installing OpenJDK 17 ..."
  sudo apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y openjdk-17-jre-headless
fi

# Install nextflow locally if not on PATH (writes to ./nextflow binary)
if ! command -v nextflow >/dev/null 2>&1 && [[ ! -x ./nextflow ]]; then
  echo "[run] Installing Nextflow into $(pwd)/nextflow ..."
  curl -fsSL get.nextflow.io | bash
fi
NF=$(command -v nextflow || echo ./nextflow)

# Pass --samplesheet / --window / etc. through to nextflow
"$NF" run main.nf "$@"
