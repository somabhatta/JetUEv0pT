#!/bin/bash
# ============================================================================
#  ana_phase2_job.sh <iacc> <cfg> <tag0,tag1,tag2> <name>
#
#  Phase2 (analysis/Phase2.C) on three merged Phase1 samples; tag = <sample>_<tree>.
#  The first tag is the reference of the stored differences dv0*. Output pT
#  grid: coarse unless V0_PTGRID is set in the environment.
#
#  Input : $PHASE1_OUT/<cfg>/acc<iacc>/<tag>_{mean,corr}.root
#  Output: $RESULTS/<cfg>/acc<iacc>/v0pt_<name>.root (+ v0pt_<name>_*.pdf)
# ============================================================================
IACC=$1; CFG=$2; TRIPLE=$3; NAME=$4
source ${PIPELINE_DIR:-/data/alice/sbhatta/JewelV0pT/pipeline}/ana_config.sh
cfg_env $CFG || exit 1
INDIR=$PHASE1_OUT/$CFG/acc$IACC
RESDIR=$RESULTS/$CFG/acc$IACC
mkdir -p $RESDIR

source $ROOT_SETUP >/dev/null 2>&1
work_dir
mkdir -p in
for tag in ${TRIPLE//,/ }; do
  for st in mean corr; do cp $INDIR/${tag}_$st.root in/ || { echo "missing ${tag}_$st"; exit 1; }; done
done

export V0_IACC=$IACC V0_PTGRID=${V0_PTGRID:-coarse}
root -l -b > p2.log 2>&1 <<ROOT
gSystem->Load("$PHASE2_LIB");
Phase2("$W/in","$TRIPLE","v0pt_$NAME.root");
.q
ROOT
tail -5 p2.log
[ -s in/v0pt_$NAME.root ] || { echo "Phase2 failed"; tail -30 p2.log; exit 1; }
cp in/v0pt_$NAME.root in/v0pt_${NAME}_*.pdf $RESDIR/ 2>/dev/null
[ -s $RESDIR/v0pt_$NAME.root ] && echo "done -> $RESDIR/v0pt_$NAME.root"
