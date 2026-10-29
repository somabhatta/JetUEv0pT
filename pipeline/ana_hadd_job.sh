#!/bin/bash
# ============================================================================
#  ana_hadd_job.sh <stage> <sample> <tree> <iacc> <cfg> <nslice>
#
#  Merge the <nslice> Phase1 parts of one stage and remove the parts. If any
#  part is missing nothing is merged and the parts are kept, so a failed slice
#  can be rerun.
#
#  Input : $PHASE1_OUT/<cfg>/acc<iacc>/parts/<sample>_<tree>_<stage>_<i>.root
#  Output: $PHASE1_OUT/<cfg>/acc<iacc>/<sample>_<tree>_<stage>.root
# ============================================================================
STAGE=$1; SAMPLE=$2; TREE=$3; IACC=$4; CFG=$5; NSLICE=$6
source ${PIPELINE_DIR:-/data/alice/sbhatta/JewelV0pT/pipeline}/ana_config.sh
TAG=${SAMPLE}_${TREE}
OUTDIR=$PHASE1_OUT/$CFG/acc$IACC
MERGED=${TAG}_${STAGE}.root

source $ROOT_SETUP >/dev/null 2>&1
work_dir

parts=()
for ((i=0; i<NSLICE; i++)); do
  f=$OUTDIR/parts/${TAG}_${STAGE}_$i.root
  [ -s $f ] || { echo "missing part $f -- not merging"; exit 1; }
  parts+=($f)
done
hadd -f $MERGED "${parts[@]}" > hadd.log 2>&1 || { echo "hadd failed"; tail -20 hadd.log; exit 1; }
rm -f $OUTDIR/$MERGED
dcache_put $MERGED $OUTDIR/$MERGED || { echo "copy to dCache failed -- parts kept"; exit 1; }
rm -f "${parts[@]}"
echo "merged ${#parts[@]} parts -> $OUTDIR/$MERGED, parts removed"
