#!/usr/bin/env bash
set -euo pipefail

SUBSET_DIR="$1"   # directory containing subset_*.vcf
THREADS="${2:-8}"

for vcf in "${SUBSET_DIR}"/subset_*.vcf; do
  base="$(basename "$vcf" .vcf)"
  phy="${SUBSET_DIR}/${base}.phy"
  out="${SUBSET_DIR}/${base}"

  # Convert VCF -> PHYLIP
  python3 scripts/snp_subsampling/vcf_to_phylip.py "$vcf" "$phy"

  # Run IQ-TREE2 (GTR+F+ASC, UFBoot and SH-aLRT)
  iqtree2 -s "$phy" -m GTR+F+ASC -T "$THREADS" -B 1000 --alrt 1000 -pre "$out"

  echo "DONE: $base"
done
