#!/usr/bin/env bash
set -euo pipefail
#source make_subsets.sh 
run_one () {
  STRAT="$1"  # random | bias
  PCT="$2"
  R="$3"
  VCF="${OUTDIR}/sub.${STRAT}.p${PCT}.${R}.vcf.gz"
  PHY="${ALIGNDIR}/sub.${STRAT}.p${PCT}.${R}.phy"
  PRE="${TREEDIR}/iq.${STRAT}.p${PCT}.${R}"

  if [ ! -s "$VCF" ]; then echo "Missing $VCF"; return; fi
  if [ ! -s "$PHY" ]; then
    python3 vcf_to_phylip.py "$VCF" > "$PHY" 2> "${LOGDIR}/convert.${STRAT}.p${PCT}.${R}.log"
  fi

  iqtree2 -s "$PHY" -m GTR+ASC -B 1000 --alrt 1000 -T AUTO -pre "$PRE"
}

export -f run_one
export OUTDIR ALIGNDIR TREEDIR LOGDIR

# Con GNU parallel
for P in $PCT_LIST; do
  parallel --will-cite -j 8 run_one {1} ${P} {2} ::: random bias ::: $(seq 1 $REPS)
done

echo "[OK] IQ-TREE finished for all subsets."

