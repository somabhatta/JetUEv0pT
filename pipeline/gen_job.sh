#!/bin/bash
# ============================================================================
#  gen_job.sh <sample> <jobid> <nevents>
#
#  One generation job of the jet-UE v0(pT) production (v2 generator setup,
#  jetv0ana v3 options): generate <nevents> events of <sample>, run jetv0ana on
#  them inside the job, copy the compact ROOT trees to dCache and delete the
#  HepMC file. Readable version of pipeline/gen_job_v2.sh, same behaviour.
#
#  Samples:
#    pp_vac   JEWEL vacuum, pp, CT18ANLO (LHAPDF 14600), ISOCHANNEL PP
#    OO_vac   JEWEL vacuum, O+O nPDF EPPS21nlo_CT18Anlo_O16 (902750), A = 16, Z = 8
#    OO_NR    JEWEL simple medium, O+O minimum bias, recoils off (energy loss only)
#    OO_RC    JEWEL simple medium, O+O minimum bias, recoils + scattering centres
#             + dummies kept (medium response; analysed with -recoil)
#    py_pp    PYTHIA 8 pp inside jetv0ana (Monash, MPI on, pT-hat^4.5 bias)
#  O+O medium samples are minimum bias (CENTRMIN 0, CENTRMAX 100, see
#  config/medium_OO_MB.dat); the sampled centrality is written into the HepMC H
#  record by the patched JEWEL (generation/jewel) and cut offline.
#
#  Output: /dcache/alice/sbhatta/JewelV0pT/trees_v2/<sample>/<sample>_<jobid 5 digits>.root
# ============================================================================
SAMPLE=$1; JOB=$2; NEV=$3

# ---- locations (production paths)
B=/data/alice/sbhatta/JewelV0pT                     # study area on stbc
JEWEL_DIR=$B/jewel-2.4.0-gcc14-cent                 # JEWEL 2.4.0, gcc 14, centrality patch
ANALYZER=$B/ana/jetv0ana                            # jetv0ana v3 (generation/jetv0ana.cc)
TABLES=$B/tables                                    # g14_* JEWEL tables (built by the gcc 14 binaries)
MEDIUM=$B/config/medium_OO_MB.dat                   # simple medium, O+O MB (also read by vacuum runs)
OUT=/dcache/alice/sbhatta/JewelV0pT/trees_v2/$SAMPLE
LOGDIR=$B/prod_v2/logs/$SAMPLE                      # JEWEL logs are kept here
ATTEMPTS=$B/prod_v2/attempts                        # per-job attempt counter (retry reseeding)
mkdir -p $OUT $LOGDIR $ATTEMPTS

# ---- software: the analyzer is built with gcc 16 (LCG_110); JEWEL must run
#      with the gcc 14 view it was built with (gcc 16 miscompiles JEWEL)
ANAENV="source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_110 x86_64-el9-gcc16-opt"
GENENV="source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107 x86_64-el9-gcc14-opt"

# Composition test (variants embMBmatch / embCmatch): mean pT of the reference
# window (0.2-3 GeV/c) in PYTHIA pp events with 40-60 GeV/c jets, R = 0.4,
# |eta| < 0.9 (measured 2026-10-28)
REFPT_MATCH=0.822

# ---- copy to dCache with verification: up to 3 copies, each verified by cmp
#      for up to 12 x 10 s (dCache may expose the file before it is complete)
dcache_put(){ local src=$1 dst=$2
  for t in 1 2 3; do
    rm -f $dst; cp $src $dst || { sleep 20; continue; }
    for k in 1 2 3 4 5 6 7 8 9 10 11 12; do cmp -s $src $dst && return 0; sleep 10; done
  done
  return 1; }

# ---- working directory: condor scratch, or a temporary directory removed on exit
if [ -n "$_CONDOR_SCRATCH_DIR" ]; then W=$_CONDOR_SCRATCH_DIR; else W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; fi
cd $W

# ---- seed: 100000 * (first 3 digits of cksum("<sample>_v2")) + jobid, shifted by
#      7919 per previous attempt so a DAGMan retry (e.g. after a JEWEL hang)
#      does not repeat the same events. The same seed drives JEWEL (NJOB) and
#      jetv0ana (-seed).
SEED=$((100000 * $(echo ${SAMPLE}_v2 | cksum | cut -c1-3) + JOB))
NATT=$(cat $ATTEMPTS/${SAMPLE}_$JOB 2>/dev/null || echo 0); echo $((NATT + 1)) > $ATTEMPTS/${SAMPLE}_$JOB
SEED=$((SEED + 7919 * NATT))
TAG=${SAMPLE}_$(printf %05d $JOB)
echo "job $TAG seed $SEED on $(hostname) in $W"

# ---- UE source lists for -ueembed: 20 consecutive files per job, starting at
#      (jobid * 20) mod nfiles (each file holds 2000 Pb+Pb events; the analyzer
#      cycles through those inside the multiplicity window)
pick(){ ls -1 $1/*.root | awk -v j=$JOB 'BEGIN{} {f[NR]=$0} END{n=NR; s=(j*20)%n; for(i=0;i<20;i++) print f[(s+i)%n+1]}'; }
pick /dcache/alice/sbhatta/Regen/Angantyr/shoving > ue_ang.list
pick /dcache/alice/sbhatta/Regen/HIJING/nominal  > ue_hij.list

# ---- analyzer options, identical for every sample (see generation/jetv0ana.cc):
#   -thermcent                      embNp   thermal UE from N_part(centrality)
#   -ueembed Ang  ... 30 60         embAng  Angantyr event, 30-60 ch. in |eta|<0.5 (O+O MB-like)
#   -ueembed AngC ... 110 160       embAngC Angantyr event, 110-160 (O+O central-like)
#   -ueembed Hij  ... 30 60         embHij  HIJING event, 30-60
#   -corr                           corr/ directory (EEC, jet-UE flow, Eq. 4 moments)
#   -thermpt 0.822                  embMBmatch, embCmatch
#   -thermflow                      embFlow (multiplicity and <pT> co-fluctuate with centrality)
#   -thermptscan 0.60 / 0.80        embMBpt060, embMBpt080
#   -spectra                        spectra/ directory (truth spectra for R_AA)
UEOPT="-thermcent -ueembed Ang ue_ang.list 30 60 -ueembed AngC ue_ang.list 110 160 -ueembed Hij ue_hij.list 30 60 -corr -thermpt $REFPT_MATCH -thermflow -thermptscan 0.60 -thermptscan 0.80 -spectra"

case $SAMPLE in
  py_pp)
    # PYTHIA runs inside the analyzer
    ( eval $ANAENV >/dev/null 2>&1; $ANALYZER -pythia $NEV -seed $SEED $UEOPT -out $TAG.root > ana.log 2>&1 ) || { tail -20 ana.log; exit 1; } ;;

  pp_vac|OO_vac|OO_NR|OO_RC)
    #            binary       LHAPDF set    A         Z         isospin              recoils
    case $SAMPLE in
      pp_vac) EXE=vac;    PDF=14600;  MASS=1;  NPROT=1; ISO="ISOCHANNEL PP"; REC=F ;;
      OO_vac) EXE=vac;    PDF=902750; MASS=16; NPROT=8; ISO="";              REC=F ;;
      OO_NR)  EXE=simple; PDF=902750; MASS=16; NPROT=8; ISO="";              REC=F ;;
      OO_RC)  EXE=simple; PDF=902750; MASS=16; NPROT=8; ISO="";              REC=T ;;
    esac
    # JEWEL parameters:
    #   NJOB = seed; sqrt(s_NN) = 5.36 TeV; hard scattering pT-hat > 5 GeV/c,
    #   |eta| < 3 for the matrix element, weighted generation with pT-hat^4.5
    #   (event weights in the HepMC E record); hadronisation on; short HepMC
    #   output; for OO_RC recoils, scattering centres (status 3) and dummies are
    #   written so that jetv0ana -recoil can subtract the thermal momenta.
    cat > params.dat <<EOF
NEVENT $NEV
NJOB $SEED
LOGFILE $W/$TAG.log
HEPMCFILE $W/$TAG.hepmc
SPLITINTFILE $TABLES/g14_splitint.dat
PDFFILE $TABLES/g14_pdfs_$PDF.dat
XSECFILE $TABLES/g14_xsecs_$PDF.dat
MEDIUMPARAMS $MEDIUM
SQRTS 5360
PDFSET $PDF
MASS $MASS
NPROTON $NPROT
$ISO
PTMIN 5.
PTMAX -1.
ETAMAX 3.
WEIGHTED T
WEXPO 4.5
HADRO T
SHORTHEPMC T
KEEPRECOILS $REC
WRITESCATCEN $REC
WRITEDUMMIES $REC
EOF
    # JEWEL with a 3 h timeout (a hang fails the job; the retry gets a new seed)
    ( eval $GENENV >/dev/null 2>&1; timeout 10800 $JEWEL_DIR/jewel-2.4.0-$EXE params.dat > jewel.out 2>&1 ) || { echo "JEWEL failed or timed out (attempt $((NATT+1)), seed $SEED)"; tail -20 jewel.out; exit 1; }
    RFLAG=""; [ $REC = T ] && RFLAG="-recoil"
    ( eval $ANAENV >/dev/null 2>&1; $ANALYZER -hepmc $TAG.hepmc $RFLAG -seed $SEED $UEOPT -out $TAG.root > ana.log 2>&1 ) || { echo "jetv0ana failed"; tail -20 ana.log; exit 1; }
    cp $TAG.log $LOGDIR/
    rm -f $TAG.hepmc                                   # only the trees are kept
    ;;

  *) echo "unknown sample $SAMPLE"; exit 1 ;;
esac
# model summary lines of the analyzer (Glauber fit, slopes, UE source usage)
grep -E "Glauber|slope|UE source|flow UE|pT scan" ana.log

dcache_put $TAG.root $OUT/$TAG.root || { echo "copy to dCache failed"; exit 1; }
echo "done -> $OUT/$TAG.root"
