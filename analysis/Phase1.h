// ============================================================================
//  Phase1.h -- per-event record (Event) and the staged event loop (Phase1).
//
//  Phase1 runs one of three stages over one slice of one tree variant, one
//  acceptance and one reference choice (Rule.h run options). Each stage is
//  Init -> Run (event loop: Fill_Evt + Fill_<stage>) -> Final:
//    Exec_ini()  -> <OutDir>/<Tag>_ini_<IJob>.root   QA distributions
//    Exec_mean() -> <OutDir>/<Tag>_mean_<IJob>.root  <f_b>, <[pT]_ref> vs Nch
//    Exec_corr() -> <OutDir>/<Tag>_corr_<IJob>.root  <f_b [pT]_ref>, C2 vs Nch
//  corr needs the MERGED mean of the same sample/config (mu(Nch) for C2).
//  Phase2.C combines the merged mean + corr files into v0(pT).
//
//  Profile arrays are indexed [replica k][probe variant][R index][jet class];
//  k = 0 is the nominal set, k = 1..NBOOT the Poisson bootstrap replicas.
// ============================================================================
#ifndef Phase1_h
#define Phase1_h
#include "Rule.h"
#include "TProfile2D.h"
#include "Tool.C"

// One tree entry unpacked for acceptance IACC (filled by Phase1::Fill_Evt).
class Event {
public:
  double weight;                              // event weight (tree branch "weight")
  double nch;                                 // Nch estimator = nchMid
  // reference particles of reference IREF, per R:
  int    refCount [RMAX];                     // N_ref
  double refSumPt [RMAX];                     // S1 = sum pT
  double refSumPt2[RMAX];                     // S2 = sum pT^2 (< 0 => absent)
  double refMeanPt[RMAX];                     // [pT]_ref = S1/N_ref (-1 if N_ref = 0)
  // leading-jet spectra; non-zero only for the jet class the event belongs to:
  bool   inJetClass[RMAX][NSEL];              // event is in jet class c and has a raw in-cone yield
  double nInCone   [RMAX][NSEL][BMAX];        // n_b : raw in-cone yield per probe pT bin
  double nPerpCone [RMAX][NSEL][BMAX];        // p_b : perpendicular-cone yield (UE estimate)
  double nCsJet    [RMAX][NSEL][BMAX];        // cs_b: yield of the constituent-subtracted jet
  Event(){ memset(this, 0, sizeof(Event)); weight = kDefaultEventWeight; }
};

class Phase1 {
private:
  string InInput, Tag, OutDir, Phase;   // file list, sample tag, output dir, stage name
  int    IJob, NEv;                     // slice index (output suffix) and events per slice (<=0: all)
  int    nR, nB;                        // number of R values, number of probe pT bins
  double Rval[RMAX];
  double ptLo[BMAX], ptHi[BMAX], ptCenter[BMAX], ptEdges[BMAX+1];   // probe pT grid

  TChain *fChain;
  TFile  *fMean, *fOut;                 // merged mean (corr stage only), stage output

  // ---- ini: QA distributions (event weighted) ------------------------------
  TH1D *hNch;                           // "hNch"            Nch estimator
  TH1D *hCentrality;                    // "hLeadHi"         tree centr (percent)
  TH1D *hRefCount   [RMAX];             // "hNref_R<i>"      N_ref
  TH1D *hRefMeanPt  [RMAX];             // "hRef_R<i>"       [pT]_ref
  TH1D *hLeadJetPt  [RMAX];             // "hNcone_R<i>"     leading-jet jetPtCorr
  TH1D *hInConeYield[RMAX];             // "hNin_R<i>"       sum_b n_b (all events)
  TH1D *hJetClass   [RMAX];             // "hJetFlag_R<i>"   jet class index, -1 = none
  TH2D *hRefMeanPtVsNch[RMAX];          // "hRefVsNch_R<i>"  [pT]_ref vs Nch

  // ---- mean stage: first moments vs fine Nch ---------------------------------
  TProfile2D *pFrac     [NBOOT+1][NVAR][RMAX][NSEL];   // <f_b>(Nch,pT)    "pN_NchPt<var>[_boot<k>]_R<i>_<sel>"
  TProfile   *pRefPt    [NBOOT+1][RMAX][NSEL];         // <[pT]_ref>(Nch), events of f      "pRef_Nch..."
  TProfile   *pRefPtFsub[NBOOT+1][RMAX][NSEL];         // <[pT]_ref>(Nch), events with sum(n-p)>0 "pRefFsub_Nch..."
  TProfile   *pRefPtFcs [NBOOT+1][RMAX][NSEL];         // <[pT]_ref>(Nch), events with sum cs>0   "pRefFcs_Nch..."

  // ---- corr stage: second moments vs fine Nch ---------------------------------
  TProfile2D *pFracRefPt[NBOOT+1][NVAR][RMAX][NSEL];   // <f_b [pT]_ref>(Nch,pT)  "pNR_NchPt<var>[_boot<k>]_R<i>_<sel>"
  TProfile   *pC2       [NBOOT+1][RMAX][NSEL];         // <c2>(Nch), mu = pRef_Nch      "pC2_Nch..."
  TProfile   *pC2Fsub   [NBOOT+1][RMAX][NSEL];         // <c2>(Nch), mu = pRefFsub_Nch  "pC2Sub_Nch..."
  TProfile   *muRefPt    [RMAX][NSEL];                 // nominal pRef_Nch from the merged mean
  TProfile   *muRefPtFsub[RMAX][NSEL];                 // nominal pRefFsub_Nch from the merged mean

  // ---- tree buffers (all acceptance slots; IACC selects the ones used) -------
  float br_weight, br_centr;
  int   br_nchMid;
  float br_jetPt[TREE_NK], br_jetPtCorr[TREE_NK], br_csJetPt[TREE_NK], br_jetPtTruth[TREE_NK];
  float br_nIn[TREE_NK*TREE_NB], br_nPerp[TREE_NK*TREE_NB], br_nCS[TREE_NK*TREE_NB];
  float br_refN[TREE_NK], br_refS1[TREE_NK], br_refS2[TREE_NK];   // reference IREF (Rule.h REFBR)

public:
  Phase1(){}
  Phase1(string input, string tag, string outdir, int ijob = 0, int nevts = -1);
  ~Phase1(){}

  void Exec_ini();
  void InitHists_ini(); void Run_ini(); void Fill_ini(Event*); void Final_ini();

  void Exec_mean();
  void InitHists_mean(); void Run_mean(); void Fill_mean(Event*, Long64_t globalEvt); void Final_mean();

  void Exec_corr();
  void ReadHists_corr();
  void InitHists_corr(); void Run_corr(); void Fill_corr(Event*, Long64_t globalEvt); void Final_corr();

private:
  void InitTree(); void ReadAxis();
  void Fill_Evt(Event*);
  int  slot(int iR) const { return IACC*TREE_NR + iR; }   // tree array slot of (IACC, iR)
};
#endif
