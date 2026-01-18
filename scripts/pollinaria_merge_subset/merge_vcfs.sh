#!/usr/bin/env bash
set -euo pipefail

SPECIES_VCF="$1"     # Filtered species VCF.gz
POLLINARIA_VCF="$2"  # Filtered pollinaria VCF.gz
OUT="$3"             # Output merged VCF.gz

# Merge to a combined multi-sample VCF
bcftools merge "$SPECIES_VCF" "$POLLINARIA_VCF" -Oz -o "$OUT"
bcftools index "$OUT"

echo "DONE: merged -> $OUT"
