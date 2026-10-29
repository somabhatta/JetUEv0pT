#!/bin/bash
# ============================================================================
#  ana_phase1_job.sh <stage> <sample> <tree> <iacc> <cfg> <islice> <nslice>
#
#  One Phase1 slice: stage ini | mean | corr of the framework in analysis/
#  (Phase1.C) for one sample, one tree variant (t_<tree>), one acceptance
#  (0: |eta|<0.9, 1: |eta|<2.5) and one analysis configuration <cfg>
#  (ana_config.sh). The slice takes every <nslice>-th tree file of the sample,
#  starting at file <islice>.
#
#  Input : $TREES/<sample>/*.root
#          corr only: merged mean $PHASE1_OUT/<cfg>/acc<iacc>/<sample>_<tree>_mean.root
#  Output: $PHASE1_OUT/<cfg>/acc<iacc>/parts/<sample>_<tree>_<stage>_<islice>.root
#          (merged and removed by ana_hadd_job.sh)
#  Env   : PIPELINE_DIR (location of ana_config.sh; set by the DAG submit files)
# ============================================================================
STAGE=$1; SAMPLE=$2; TREE=$3; IACC=$4; CFG=$5; ISLICE=$6; NSLICE=$7
source ${PIPELINE_DIR:-/data/alice/sbhatta/JewelV0pT/pipeline}/ana_config.sh
cfg_env $CFG || exit 1
TAG=${SAMPLE}_${TREE}
OUTDIR=$PHASE1_OUT/$CFG/acc$IACC
mkdir -p $OUTDIR/parts

source $ROOT_SETUP >/dev/null 2>&1
work_dir

# input files of this slice
ls -1 $TREES/$SAMPLE/*.root | awk -v i=$ISLICE -v n=$NSLICE '(NR-1)%n==i' > files.list
[ -s files.list ] || { echo "no input files for slice $ISLICE"; exit 1; }

# run options (Rule.h); the event offset keeps bootstrap seeds distinct per slice
export V0_TREE=$TREE V0_IACC=$IACC
export V0_EVTOFFSET=$((ISLICE * 1000000000))
if [ $STAGE = corr ]; then
  cp $OUTDIR/${TAG}_mean.root mean.root || { echo "merged mean missing"; exit 1; }
  export V0_MEANFILE=$W/mean.root
fi

mkdir -p out
root -l -b > p1.log 2>&1 <<ROOT
gSystem->Load("$PHASE1_LIB");
Phase1 p("$W/files.list","$TAG","$W/out",$ISLICE,-1);
p.Exec_$STAGE();
.q
ROOT
tail -5 p1.log
PART=${TAG}_${STAGE}_${ISLICE}.root
[ -s out/$PART ] || { echo "Phase1 $STAGE failed"; tail -30 p1.log; exit 1; }
dcache_put out/$PART $OUTDIR/parts/$PART || { echo "copy to dCache failed"; exit 1; }
echo "done -> $OUTDIR/parts/$PART"
