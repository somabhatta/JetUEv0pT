# ============================================================================
#  ana_config.sh -- shared settings of the analysis job scripts
#  (ana_phase1_job.sh, ana_hadd_job.sh, ana_phase2_job.sh); sourced, not run.
#
#  Analysis configuration name CFG = <ref>_<cent>_jt<n>, e.g.
#     refnom_call_jt0     nominal reference, no centrality cut, jetPtCorr classes
#     refnom_c0-40_jt0    same, centrality 0-40 %
#     refaw_call_jt2      away-side-excluded reference, jetPtTruth classes
#   <ref>  refnom refgap refaw refnear reftrans refaway  -> V0_IREF = 0..5 (Rule.h REFN)
#   <cent> call = no cut, c<lo>-<hi> = centrality lo <= centr < hi (percent)
#                                     -> V0_CENTLO, V0_CENTHI
#   jt<n>  jet-pT definition of the jet classes -> V0_JETPT
#          0 jetPtCorr, 1 csJetPt, 2 jetPtTruth
# ============================================================================
B=/data/alice/sbhatta/JewelV0pT        # work area (results, logs)
D=/dcache/alice/sbhatta/JewelV0pT      # dCache: input trees and Phase1 outputs
TREES=${TREES:-$D/trees_v2}            # input trees: $TREES/<sample>/*.root
PHASE1_OUT=${PHASE1_OUT:-$D/phase1_v3} # Phase1 outputs: $PHASE1_OUT/<cfg>/acc<i>/
RESULTS=${RESROOT:-$B/results_v3}      # Phase2 outputs: $RESULTS/<cfg>/acc<i>/

# Compiled framework (analysis/, ACLiC). The production used Phase1 built with
# -DV0_NBOOT=40 and Phase2 built with the default NBOOT=200 (Phase2 skips
# replicas missing from its input; the v0boot* TH2D then have 201 y bins, of
# which bins 1..41 are filled). Override to use another build.
PHASE1_LIB=${PHASE1_LIB:-$B/v0pT_v3_b40/Phase1_C.so}
PHASE2_LIB=${PHASE2_LIB:-$B/v0pT_v3/Phase2_C.so}

ROOT_SETUP="/cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107 x86_64-el9-gcc14-opt"

REFN=(refnom refgap refaw refnear reftrans refaway)

# cfg_env <cfg> : export V0_IREF, V0_CENTLO/V0_CENTHI (or unset), V0_JETPT
cfg_env(){
  local cfg=$1 ref cent jt i range
  IFS=_ read -r ref cent jt <<< "$cfg"
  export V0_IREF=-1
  for i in "${!REFN[@]}"; do [ "${REFN[$i]}" = "$ref" ] && export V0_IREF=$i; done
  [ $V0_IREF -ge 0 ] || { echo "bad cfg $cfg"; return 1; }
  unset V0_CENTLO V0_CENTHI
  if [ "$cent" != call ]; then
    range=${cent#c}
    export V0_CENTLO=${range%-*} V0_CENTHI=${range#*-}
  fi
  export V0_JETPT=${jt#jt}
}

# dcache_put <src> <dst> : copy and verify (cmp) with retries; dCache can show
# a file before it is complete, so the comparison is retried for up to 2 min.
dcache_put(){
  local src=$1 dst=$2 try k
  for try in 1 2 3; do
    rm -f $dst; cp $src $dst || { sleep 20; continue; }
    for k in 1 2 3 4 5 6 7 8 9 10 11 12; do cmp -s $src $dst && return 0; sleep 10; done
  done
  return 1
}

# work_dir : cd into the condor scratch dir, or a temporary dir removed at exit
work_dir(){
  if [ -n "$_CONDOR_SCRATCH_DIR" ]; then W=$_CONDOR_SCRATCH_DIR
  else W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; fi
  cd $W
}
