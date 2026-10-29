// ============================================================================
//  Rule.h -- constants, tree layout and run options of the jet--UE v0(pT)
//  analysis (Phase1 = staged ini/mean/corr event loop, Phase2 = combination).
//
//  Observable (per leading jet, per event):
//      f_b      = n_b / N            in-cone charged-particle fraction in probe
//                                    pT bin b (N = sum over b of n_b)
//      [pT]_ref = S1 / N_ref          mean pT of the reference particles
//      v0(pT_b) = Cov(f_b, [pT]_ref) / ( <f_b> sqrt(C2,ref) )
//  Everything is accumulated in fine Nch bins (mid-rapidity multiplicity) and
//  combined over Nch in Phase2 (average of per-Nch-bin ratios).
//
//  Input : jetv0ana trees "t_<V0_TREE>" (one tree per particle-list variant:
//          nom, cs10, embMB, embFlow, cs10_embFlow, ...). Per event, arrays over
//          (acceptance, R) slots and, for the in-cone spectra, the probe pT grid.
//
//  Run options (environment variables, read once at library load):
//    V0_TREE      tree variant, tree name = "t_<V0_TREE>"           (default nom)
//    V0_IACC      acceptance slot: 0 = |eta|<0.9, 1 = |eta|<2.5      (default 0)
//    V0_IREF      reference region index into REFN[] below          (default 0)
//    V0_CENTLO/HI centrality class V0_CENTLO <= centr < V0_CENTHI on the tree's
//                 centr (percent); events with centr < 0 always pass (default:
//                 no cut)
//    V0_JETPT     jet-pT definition of the jet classes: 0 jetPtCorr,
//                 1 csJetPt, 2 jetPtTruth                           (default 0)
//    V0_NCHBINS   number of fine Nch bins                          (default 100)
//    V0_NCHMAX    upper edge of the Nch axis                       (default 500)
//    V0_EVTOFFSET global-event offset of this slice (bootstrap seed) (default 0)
//    V0_MEANFILE  merged mean file read by the corr stage (Phase1_corr.C)
//    V0_PTGRID    "coarse" = coarse output pT grid in Phase2 (Phase2.C)
//  Compile-time option:
//    -DV0_NBOOT=n number of Poisson bootstrap replicas (default 200;
//                 the v3 production uses 40)
// ============================================================================
#ifndef Rule_H_
#define Rule_H_
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TNamed.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TRandom3.h"
using namespace std;

enum { RMAX = 16, BMAX = 64 };   // array caps: number of R values, number of probe pT bins

inline int    EnvInt   (const char* var, int    def){ const char* e = getenv(var); return e ? atoi(e) : def; }
inline double EnvDouble(const char* var, double def){ const char* e = getenv(var); return e ? atof(e) : def; }

// ---- jetv0ana tree layout (keep in sync with ana/jetv0ana.cc) ---------------
// Per-event arrays of length TREE_NK are indexed by the slot
//     k = iAcc*TREE_NR + iR ;
// arrays of length TREE_NK*TREE_NB (in-cone spectra) by k*TREE_NB + iPt.
enum { TREE_NR = 4, TREE_NA = 2, TREE_NB = 25, TREE_NK = TREE_NR*TREE_NA };
static const double TREE_RVALS [TREE_NR]   = {0.1, 0.2, 0.4, 0.8};   // anti-kT R
static const double TREE_ACCS  [TREE_NA]   = {0.9, 2.5};             // |eta| acceptance
static const double TREE_PTEDGE[TREE_NB+1] = {0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0,
                                              5.0, 6.0, 8.0, 10., 12., 15., 20., 25., 30., 40.,
                                              50., 60., 80., 100., 150., 300.};   // probe pT grid (GeV/c)

// acceptance analysed by this run (V0_IACC) and the tree variant (V0_TREE)
int IACC = []{ int v = EnvInt("V0_IACC", 0); return (v >= 0 && v < TREE_NA) ? v : 0; }();
inline string treeName(){ const char* e = getenv("V0_TREE"); return string("t_") + (e ? e : "nom"); }

// ---- wide centrality classes in Phase2 (Nch cuts) ---------------------------
// Only the inclusive class is used: NCENT = 1 wide class + the inclusive class
// (index ICENT_INCL = NCENT, written as "_c1_"). The wide class (index 0) is
// formed only if Phase2 finds INPUT/cutNch_<sample>.txt; the production has no
// such file. Event-level centrality selection is done in Phase1 instead
// (V0_CENTLO/HI, see below).
enum { NCENT = 1, NCENTX = NCENT + 1 };
static const int ICENT_INCL = NCENT;

// ---- fine Nch axis ----------------------------------------------------------
// Nch = tree branch nchMid (charged particles |eta|<0.8, 0.2<pT<5 GeV, full
// event). All moments are accumulated per fine Nch bin and combined in Phase2
// as an event-weighted average of per-bin ratios, which removes the bias from
// multiplicity (volume) fluctuations. 100 bins of width 5 up to 500 keep a
// 40-replica mean/corr job at ~0.3 GB.
int    NCHBINS = EnvInt   ("V0_NCHBINS", 100);
double NCHMAX  = EnvDouble("V0_NCHMAX",  500.);

// ---- jet classes ("selections") ---------------------------------------------
// An event enters jet class c for radius R if its leading jet (this R, this
// acceptance) has jet pT (definition V0_JETPT) in [SELLO[c], SELHI[c]).
// SELN[c] is the class name in every output object. ISEL_JET (20-40 GeV/c) is
// used for the Phase2 quick-look plots; 10-20 is exploratory.
enum { NSEL = 5 };
static const char*  SELN [NSEL] = {"jet10_20", "jet20_40", "jet40_60", "jet60_80", "jet80_120"};
static const double SELLO[NSEL] = { 10.,        20.,        40.,        60.,        80.       };
static const double SELHI[NSEL] = { 20.,        40.,        60.,        80.,        120.      };
static const int ISEL_JET = 1;

// Phase2 output pT grid (fine) = the probe grid.
enum { NCEDGE = TREE_NB + 1 };
static const double* CEDGES = TREE_PTEDGE;

// ---- reference particles ([pT]_ref) -----------------------------------------
// All references use 0.2 < pT < 3 GeV/c inside the acceptance and differ by the
// excluded region (tree branches <REFBR>N, <REFBR>S1, <REFBR>S2 = count, sum pT,
// sum pT^2 per slot):
//   0 refnom   : outside the jet cone and both perpendicular cones
//   1 refgap   : |eta - eta_jet| > R + 0.5
//   2 refaw    : refnom without the away side |dphi - pi| < 1
//   3 refnear  : refnom, |dphi| < pi/3
//   4 reftrans : refnom, pi/3 < |dphi| < 2pi/3
//   5 refaway  : refnom, |dphi| > 2pi/3        (near + trans + away = refnom)
// One run analyses ONE reference (V0_IREF); REFN[IREF] names its output dir.
enum { NREF = 6 };
static const double PTREF_LO[NREF] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
static const double PTREF_HI[NREF] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
static const char*  REFN    [NREF] = {"refnom", "refgap", "refaw", "refnear", "reftrans", "refaway"};
static const char*  REFBR   [NREF] = {"ref",    "gap",    "refAw", "near",    "trans",    "away"};
int IREF = []{ int v = EnvInt("V0_IREF", 0); return (v >= 0 && v < NREF) ? v : 0; }();

// ---- reference-window stamp (TNamed "REFWIN") -------------------------------
// IREF comes from the environment of the shell that ran a stage, not from the
// data, so the mean stage stamps the reference it used into its output:
// title = "<ptlo> <pthi> <iref> <name>". TNamed has no Merge(), so hadd keeps
// one copy (it is a per-job constant).
inline void WriteRefWin(double ptLo, double ptHi, int iref, const char* refName){
  TNamed stamp("REFWIN", Form("%.2f %.2f %d %s", ptLo, ptHi, iref, refName));
  stamp.Write();   // writes into gDirectory: cd() to the output file first
}
inline void WriteRefWin(){ WriteRefWin(PTREF_LO[IREF], PTREF_HI[IREF], IREF, REFN[IREF]); }
inline bool ReadRefWin(TFile* f, double& ptLo, double& ptHi, int& iref, TString& refName){
  if (!f || f->IsZombie()) return false;
  TNamed* stamp = (TNamed*)f->Get("REFWIN"); if (!stamp) return false;
  char buf[64];
  if (sscanf(stamp->GetTitle(), "%lf %lf %d %63s", &ptLo, &ptHi, &iref, buf) != 4) return false;
  refName = buf; return true;
}

// ---- probe definitions --------------------------------------------------------
// Per event and jet class, with n_b the raw in-cone yield in probe bin b, p_b
// the perpendicular-cone yield (average of the cones at phi_jet +- pi/2) and
// cs_b the yield of the leading jet after event-wise constituent subtraction:
//   IVAR_F    "f"    : f_b     = n_b / sum n                    (raw)
//   IVAR_FSUB "f_sub": f_sub,b = (n_b - p_b) / sum (n - p)       (perp-cone subtracted)
//   IVAR_FCS  "f_cs" : f_cs,b  = cs_b / sum cs                   (constituent subtracted)
// VARSUF is the name suffix of the variant in every output object (f keeps the
// unsuffixed names). f_sub and f_cs exist only for events with a positive
// total; their reference mean <[pT]_ref> is averaged over exactly that event
// subset (pRefFsub_Nch, pRefFcs_Nch) so that sum_b Cov(f_b,[pT]_ref) = 0 holds
// exactly for every variant (sum_b f_b = 1).
enum { IVAR_F = 0, IVAR_FSUB = 1, IVAR_FCS = 2, NVAR = 3 };
static const char* VARSUF[NVAR] = {"",  "_fsub", "_fcs"};
static const char* VARLAB[NVAR] = {"f", "fsub",  "fcs"};

static const double kDefaultEventWeight = 1.0;   // the tree branch "weight" overrides it

// ---- Poisson bootstrap ---------------------------------------------------------
// The mean and corr stages fill NBOOT+1 parallel copies of every profile in the
// same event loop: index 0 = nominal (event weight w), index k = 1..NBOOT uses
// w * Poisson(1) with the draw seeded by (global event index, k). Consequences:
//  * mean replica k and corr replica k reweight each event identically, so
//    <f>_k, mu_k, <f [pT]_ref>_k and C2_k fluctuate together;
//  * the seed does not depend on the job: rerunning a job reproduces it, and
//    hadd of per-slice outputs equals one full-sample run;
//  * tree variants of the SAME generated events (e.g. nom vs embFlow) get the
//    same weights in replica k, so replica-wise differences give correlated
//    errors (Phase2 v0boot*, dv0*).
#ifndef V0_NBOOT
#define V0_NBOOT 200
#endif
static const int NBOOT = V0_NBOOT;

inline double PoissonBootWeight(Long64_t globalEvt, int replica){
  ULong64_t seed = (ULong64_t)globalEvt*2654435761ULL + (ULong64_t)replica*40503ULL + 1ULL;
  static TRandom3 rng;
  rng.SetSeed((UInt_t)(seed ^ (seed >> 32)));
  return (double)rng.Poisson(1.0);
}

// Global-event offset of this slice (V0_EVTOFFSET). File-split slices start
// their chain at entry 0; each slice adds its own offset (pipeline: slice*1e9)
// so bootstrap weights differ between slices and agree between the mean and
// corr job of the same slice.
inline Long64_t evtOffset(){ const char* e = getenv("V0_EVTOFFSET"); return e ? atoll(e) : 0; }

// ---- event selection -------------------------------------------------------------
// Centrality: V0_CENTLO <= centr < V0_CENTHI (percent; JEWEL-sampled for O+O
// medium, drawn per event in the embedded vacuum/PYTHIA variants). Events with
// centr < 0 (no added UE) pass every class: nothing in them depends on
// centrality, so they serve as the baseline of each class.
double CENTLO = EnvDouble("V0_CENTLO", -1e9);
double CENTHI = EnvDouble("V0_CENTHI",  1e9);
// Jet-pT definition used for the jet classes (V0_JETPT).
int    JETPTDEF = []{ int v = EnvInt("V0_JETPT", 0); return (v >= 0 && v <= 2) ? v : 0; }();
static const char* JETPTN[3] = {"jetPtCorr",    // raw pT - rho_perp * A (perp-cone UE correction)
                                "csJetPt",      // jet of the constituent-subtracted event
                                "jetPtTruth"};  // leading jet of the same event without added UE

// Event range of job slice ijob of nev events (nev <= 0: all events).
inline void jobRange(Long64_t nEntries, int ijob, int nev, Long64_t& first, Long64_t& last){
  if (nev <= 0){ first = 0; last = nEntries; return; }
  first = (Long64_t)ijob*nev; last = first + nev;
  if (last > nEntries) last = nEntries;
  if (first > nEntries) first = nEntries;
}
#endif
