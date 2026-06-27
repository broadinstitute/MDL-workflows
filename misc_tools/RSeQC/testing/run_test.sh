#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repo_root}"

awk -f RSeQC/testing/gtf_to_bed12.awk \
  RSeQC/testing/data/SIRVs1-7.annot.reduced.gtf \
  | sort -k1,1 -k2,2n > RSeQC/testing/data/SIRVs1-7.annot.reduced.bed12

mkdir -p RSeQC/testing/runs
miniwdl run RSeQC/GeneBodyCoverage.wdl \
  -i RSeQC/testing/inputs.json \
  --dir RSeQC/testing/runs
