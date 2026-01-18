#!/usr/bin/env bash
set -euo pipefail


source make_subsets_idempotent.sh

# Set to FORCE=true 
FORCE="${FORCE:-false}"

command -v python3 >/dev/null || { echo "[ERR] python3 missing"; exit 1; }
command -v iqtree2  >/dev/null || { echo "[ERR] iqtree2 missing"; exit 1; }

mkdir -p "$ALIGNDIR" "$TREEDIR" "$LOGDIR"

# Trova un VCF anche se la percentuale è formattata 0.5 vs 0.50
find_vcf() {
  local strat="$1" pct="$2" rep="$3"
  local cand="${OUTDIR}/sub.${strat}.p${pct}.${rep}.vcf.gz"
  [[ -s "$cand" ]] && { echo "$cand"; return 0; }
  local base="${OUTDIR}/sub.${strat}.p"
  local try1="${base}${pct}.${rep}.vcf.gz"
  local try2="${base}${pct%0}.${rep}.vcf.gz"  # 0.50 -> 0.5
  local try3="${base}${pct}"*"."${rep}".vcf.gz"
  for g in "$try1" "$try2" $try3; do
    for f in $g 2>/dev/null; do
      [[ -s "${f:-}" ]] && { echo "$f"; return 0; }
    done
  done
  return 1
}

iqtree_done() {
  local pre="$1"

  [[ -s "${pre}.treefile" || -s "${pre}.contree" || -s "${pre}.iqtree" ]]
}

run_one () {
  local STRAT="$1" PCT="$2" R="$3"
  local VCF PHY PRE LOGC

  if ! VCF="$(find_vcf "$STRAT" "$PCT" "$R")"; then
    echo "[WARN] missing VCF: strat=${STRAT} pct=${PCT} rep=${R}"
    return 0
  fi

  PHY="${ALIGNDIR}/sub.${STRAT}.p${PCT}.${R}.phy"
  PRE="${TREEDIR}/iq.${STRAT}.p${PCT}.${R}"
  LOGC="${LOGDIR}/convert.${STRAT}.p${PCT}.${R}.log"

  # --- VCF -> PHYLIP 
  if [[ "$FORCE" != "true" && -s "$PHY" ]]; then
    echo "[SKIP] PHYLIP exists: $PHY"
  else
    echo "[DO]   Convert VCF→PHYLIP: $(basename "$VCF") -> $(basename "$PHY")"
    python3 ./vcf_to_phylip.py "$VCF" > "$PHY" 2> "$LOGC"
  fi

  # --- IQ-TREE 
  if [[ "$FORCE" != "true" && $(iqtree_done "$PRE") ]]; then
    echo "[SKIP] IQ-TREE already done: $PRE.(treefile/contree/iqtree)"
  else
    echo "[DO]   IQ-TREE: $(basename "$PHY")  →  prefix=$(basename "$PRE")"
    iqtree2 -s "$PHY" -m GTR+ASC -B 1000 --alrt 1000 -T AUTO -pre "$PRE"
  fi
}

# Loop 
for P in $PCT_LIST; do
  for STRAT in random bias; do
    for R in $(seq 1 "$REPS"); do
      run_one "$STRAT" "$P" "$R"
    done
  done
done

echo "[OK] IQ-TREE finished (idempotent)."
