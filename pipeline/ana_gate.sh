#!/bin/bash
# ============================================================================
#  ana_gate.sh -- first node of the analysis DAG: wait (up to 10 h, 5-min
#  polls) until every sample has its 500 tree files on dCache, then exit 0.
# ============================================================================
TREES=/dcache/alice/sbhatta/JewelV0pT/trees_v2
SAMPLES="pp_vac OO_vac OO_NR OO_RC py_pp"
NFILES=500
for poll in $(seq 1 120); do
  complete=1
  for s in $SAMPLES; do
    n=$(ls $TREES/$s/*.root 2>/dev/null | wc -l); [ $n -ge $NFILES ] || complete=0
  done
  [ $complete = 1 ] && { echo "all samples complete"; exit 0; }
  sleep 300
done
echo "timeout: samples incomplete"
for s in $SAMPLES; do echo "$s $(ls $TREES/$s/*.root | wc -l)"; done
exit 1
