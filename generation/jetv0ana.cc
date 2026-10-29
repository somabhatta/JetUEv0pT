// ============================================================================
//  jetv0ana.cc -- per-event building blocks for the jet--underlying-event
//                 v0(pT) study (JEWEL 2.4.0 / PYTHIA 8, pp and O+O at 5.36 TeV)
//
//  The program reads generator events, builds several particle lists per event
//  ("variants": the event as generated, the event with recoil thermal momenta
//  subtracted, the event with an added underlying event, ...) and writes, for
//  every variant, one TTree entry per event with the leading-jet, background and
//  reference-region quantities that the offline v0(pT) analysis needs.
//
//  This file is the readable version of the production analyzer ("v3",
//  2026-10-28). All tree, branch, directory and histogram names, the order of
//  every random-number call and every arithmetic expression are unchanged, so
//  the output is bit-identical to the production binary.
//
// ----------------------------------------------------------------------------
//  COMMAND LINE
// ----------------------------------------------------------------------------
//  Input (exactly one of):
//    -hepmc FILE        JEWEL HepMC2 ASCII (IO_GenEvent) file
//    -pythia N          generate N PYTHIA 8 events in process (needs -DWITH_PYTHIA8):
//                       HardQCD:all, pT-hat > pthatmin, pT-hat^4.5 bias (same
//                       weighting as JEWEL WEXPO 4.5), Monash default tune, MPI on
//  Input modifiers:
//    -recoil            JEWEL sample with KEEPRECOILS/WRITESCATCEN/WRITEDUMMIES T:
//                       write the recoil-subtraction variants nosub, cs10, cs05
//    -ecm E             PYTHIA centre-of-mass energy in GeV       (default 5360)
//    -pthatmin P        PYTHIA minimum pT-hat in GeV/c            (default 5)
//    -nompi             PYTHIA with PartonLevel:MPI = off
//    -nev N             stop after N events (default: all)
//  Output and random numbers:
//    -out FILE          output ROOT file                          (default jetv0.root)
//    -seed S            seed of the random generators (thermal UE, drawn
//                       centrality, PYTHIA, start of the UE-file cycle; default 1)
//  Added underlying-event (UE) variants (all optional; production uses all):
//    -thermcent         embNp: thermal UE with dN/deta from N_part(centrality)
//    -thermflow         embFlow: thermal UE with co-fluctuating dN/deta and <pT>
//                       (requires -thermcent)
//    -thermpt X         embMBmatch, embCmatch: thermal UE with mean pT X (GeV/c)
//                       inside the reference window 0.2-3 GeV/c
//    -thermptscan X     embMBptNNN (NNN = round(100 X)): thermal UE at dN/deta 42
//                       with window mean pT X; repeatable
//    -ueembed TAG LIST LO HI
//                       emb<TAG>: embed one event from the CompactEvent files in
//                       LIST (Angantyr, HIJING) whose charged multiplicity in
//                       |eta| < 0.5 lies in [LO, HI]; repeatable
//  Extra output directories:
//    -corr              directory "corr"    (EEC, jet-UE flow, Eq. 4 moments)
//    -spectra           directory "spectra" (truth spectra for R_AA / R_CP)
//  Event display (off by default; never changes the main output):
//    -display FILE      write the intermediate objects of the first events to the
//                       separate ROOT file FILE (see EVENT DISPLAY below)
//    -displayEvents N   number of events written to FILE          (default 5)
//    -displayVariants A,B,...
//                       variants written to FILE (default: the no-UE variant, nom
//                       or cs10, and the thermal MB variant, embMB or cs10_embMB)
//
//  Production call (pipeline/gen_job.sh), with REFPT_MATCH = 0.822:
//    jetv0ana (-pythia N | -hepmc F [-recoil]) -seed SEED
//      -thermcent -ueembed Ang ue_ang.list 30 60 -ueembed AngC ue_ang.list 110 160
//      -ueembed Hij ue_hij.list 30 60 -corr -thermpt 0.822 -thermflow
//      -thermptscan 0.60 -thermptscan 0.80 -spectra -out TAG.root
//
// ----------------------------------------------------------------------------
//  VARIANTS (one TTree "t_<variant>" each, in this order)
// ----------------------------------------------------------------------------
//  Samples without recoils (pp_vac, OO_vac, OO_NR, py_pp): base list "nom".
//  Recoil sample (OO_RC, -recoil): base lists "nosub", "cs10", "cs05" and every
//  added-UE variant is built on cs10 and carries the prefix "cs10_".
//
//    nom          status-1 final state as generated (JEWEL dummies e < 1e-5 dropped)
//    nosub        recoil sample, status-1 particles, no thermal subtraction
//    cs10, cs05   recoil sample, thermal momenta of the scattering centres removed
//                 by constituent subtraction, dR_max = 1.0 / 0.5 (port of the
//                 JEWEL 2.6 subtraction.C)
//    embMB, embC  + thermal charged pions, dN/deta = 42 (O+O MB) / 135 (O+O 0-5%),
//                 Gamma(2,T) spectrum with <pT> = 0.7 GeV/c (window <pT> 0.767)
//    embNp        + thermal pions, dN/deta = A * N_part(c)^beta (optical Glauber
//                 16O+16O), c = JEWEL medium centrality, or drawn uniformly in
//                 0-100% for vacuum / PYTHIA events
//    embMBmatch, embCmatch
//                 as embMB / embC but with window <pT> = -thermpt value (0.822,
//                 the PYTHIA jet-event reference window <pT>)
//    embAng, embAngC, embHij
//                 + one Angantyr (shoving) event with 30-60 / 110-160 charged
//                 particles in |eta| < 0.5, or one HIJING event with 30-60
//    embFlow      + thermal pions, dN/deta and window <pT> co-fluctuating with c
//                 (see FlowingUEModel)
//    embMBpt060, embMBpt080
//                 + thermal pions at dN/deta = 42 with window <pT> 0.60 / 0.80
//
// ----------------------------------------------------------------------------
//  TREE BRANCHES (identical in every t_<variant>; one entry per event)
// ----------------------------------------------------------------------------
//  Array index k = iacc * 4 + iR, with iacc 0: |eta_track| < 0.9 (ALICE),
//  1: |eta_track| < 2.5 (wide), and iR 0..3: R = 0.1, 0.2, 0.4, 0.8.
//  Probe-pT histograms are flattened as [k * 25 + bin], bin edges kProbePtEdges.
//  Tracks: charged, status 1, pT > 0.15 GeV/c, |eta| < acceptance; massless.
//  Quantities of a missing leading jet (or of an acceptance without tracks) are 0.
//
//    weight     F       event weight (JEWEL HepMC weight / PYTHIA info.weight),
//                       carries the pT-hat^4.5 bias
//    centr      F       centrality in percent: JEWEL medium sampled centrality;
//                       -1 for vacuum JEWEL and PYTHIA; for embNp / embFlow the
//                       centrality that set the UE (drawn if the event had none)
//    nchMid     I       charged particles |eta| < 0.8, 0.2 <= pT < 5 GeV/c
//    nchFwdA    I       charged particles 3.5 < eta < 4.9,  0.2 < pT < 10 GeV/c (FT0A-like)
//    nchFwdC    I       charged particles -3.3 < eta < -2.1, 0.2 < pT < 10 GeV/c (FT0C-like)
//                       (trees of the 2026-10-28 production have no pT cut on these two)
//    rhoEv[2]   F       per acceptance: rho (GeV/c per unit area) of the R = 0.4
//                       perpendicular cones; background of the constituent
//                       subtraction (0 if there is no R = 0.4 jet)
//    ueMult     F       added UE: thermal dN/deta used, or the number of charged
//                       particles in |eta| < 0.5 of the embedded event; 0 if none
//    jetPt      F[8]    leading anti-kT jet pT (GeV/c), |eta_jet| < acc - R
//    jetPtCorr  F[8]    jetPt - rho * jetArea (GeV/c)
//    jetEta, jetPhi F[8] leading jet eta, phi in (-pi, pi]
//    jetArea    F[8]    active jet area (explicit ghosts, ghost area 0.05)
//    leadTrk    F[8]    highest constituent pT of the leading jet (GeV/c)
//    rho        F[8]    perpendicular-cone pT density for this R (GeV/c / area)
//    csJetPt, csJetEta, csJetPhi, csLeadTrk F[8]
//                       leading jet after event-wise constituent subtraction
//    nIn        F[200]  leading-jet constituents per probe-pT bin
//    nPerp      F[200]  tracks per probe-pT bin in the two perpendicular cones,
//                       averaged over the two cones (0.5 per track)
//    nCS        F[200]  CS-jet constituents per probe-pT bin
//    refN, refS1, refS2 F[8]
//                       reference region: tracks with 0.2 <= pT < 3 GeV/c outside
//                       the jet cone (constituents and dR < R to the axis) and
//                       outside both perpendicular cones: count, sum pT (GeV/c),
//                       sum pT^2 (GeV^2/c^2)
//    refNa..refS2a, refNb..refS2b F[8]
//                       same, split into sub-events a: eta < 0, b: eta >= 0
//    gapN, gapS1, gapS2 F[8]
//                       eta-gap reference: 0.2 <= pT < 3, outside the jet cone,
//                       |eta - eta_jet| > R + 0.5 (perpendicular cones NOT removed)
//    refAwN..refAwS2b F[8]
//                       reference region with the away side |dphi - pi| < 1
//                       removed (kept: |dphi| <= pi - 1), full / a / b
//    nearN/S1/S2   F[8] reference region with |dphi| < pi/3
//    transN/S1/S2  F[8] reference region with pi/3 <= |dphi| < 2pi/3
//    awayN/S1/S2   F[8] reference region with |dphi| >= 2pi/3
//                       (near + trans + away = ref, per event)
//    jetPtTruth F[8]    jetPt of the same event without added UE (variant nom,
//                       or cs10 for the recoil sample)
//    dRTruth    F[8]    distance (eta, phi) between this leading jet and the
//                       no-UE leading jet; -1 if either is missing
//
// ----------------------------------------------------------------------------
//  DIRECTORIES
// ----------------------------------------------------------------------------
//  corr/    (-corr) per variant, acceptance a<iacc>, R index R2 (0.4) / R3 (0.8),
//           jet class j20_40 / j40_60 GeV/c in jetPtCorr; all event weighted,
//           sums only (hadd-able). Names "<variant>_a<i>_R<ir>_j<lo>_<hi>_<h>":
//             eec       sum_{i<j in jet} pT,i pT,j / pT,jet^2 vs dR_ij (log bins)
//             flowXS    X * pT_j, flowXN X, flowS pT_j, flowN 1, vs angle of
//                       reference track j to the jet axis
//                       (X = sum of constituent pT / pT,jet)
//             eq4_nS, eq4_nNr, eq4_n, eq4_fr, eq4_f
//                       per probe bin b: n_b S, n_b N_r, n_b, f_b r, f_b
//                       (n_b raw in-cone counts, f_b = n_b / N, S / N_r sum pT /
//                       count of the reference tracks, r = S / N_r)
//             scal      bins: sum w, N S, N N_r, N, S, N_r, r, X, n_jets
//  spectra/ (-spectra) truth level (no added UE; recoil sample: cs10), per
//           centrality class c0_10, c0_40, c40_80, c0_100 (events without a
//           centrality only enter c0_100), event weighted:
//             trk_<c>        charged-hadron pT, |eta| < 0.9, pT >= 0.2 GeV/c
//             jet_R02_<c>, jet_R04_<c>
//                            inclusive charged jets pT > 5, |eta_jet| < 0.9 - R
//             norm_<c>       bins: sum w, number of events, sum w^2
//
// ----------------------------------------------------------------------------
//  EVENT DISPLAY (-display FILE; read by plotting/PlotEventDisplay.C)
// ----------------------------------------------------------------------------
//  The display is recorded inside the production functions (selectChargedTracks,
//  findLeadingJet, fillPerpendicularCones, fillReferenceSums, the constituent
//  subtraction, subtractRecoilThermalMomenta) by passing them an optional record
//  pointer; it only reads objects that exist anyway (the jet ghosts are read from
//  the ClusterSequenceArea, the CS ghosts from the ConstituentSubtractor) and
//  draws no random numbers, so the main output is unchanged.
//  Slots written: k = 2 (|eta_trk| < 0.9, R = 0.4) and k = 7 (|eta_trk| < 2.5,
//  R = 0.8).
//    tree "display", one entry per (event, displayed variant, slot):
//      event, variant, iacc, iR, R, trackEtaMax, jetEtaMax (= trackEtaMax - R),
//      weight, centr, ueMult, the tree quantities of slot k (leadPt = jetPt,
//      leadPtCorr, leadEta, leadPhi, leadArea, leadTrk, rho, rhoCS = rhoEv,
//      csJetPt/Eta/Phi, jetPtTruth, dRTruth, refN, refS1), perpPhi[2],
//      probeEdges[26], per probe bin b: nIn (n_b), nPerp (b_b), nSub (n_b - b_b),
//      nCS, fRaw = n_b / sum n, fSub = (n_b - b_b) / sum (n - b) (0 if
//      sum (n - b) <= 0), fCS; jetGhostArea, csGhostArea, csMaxDR;
//      tracks entering the clustering (index = PseudoJet user_index):
//        trkPt, trkEta, trkPhi, trkOrigin (0 generator, 1 added UE, 2 generator
//        particle whose pT was reduced by the recoil thermal subtraction),
//        trkJet (index into jet*, -1: in no displayed jet), trkFlags (bits
//        kDisplayFlag*: 1 leading-jet constituent, 2 within R of the jet axis,
//        4 perpendicular cone, 8 reference, 16 reference without away side,
//        32 eta-gap reference, 64 CS-jet constituent), trkPtCS (pT after the
//        event-wise CS, 0 if removed);
//      anti-kT jets with pT > 3 GeV/c (and the leading jet): jetPt, jetEta,
//        jetPhi, jetArea, jetNConst, jetFiducial (|eta_jet| < jetEtaMax);
//        leadJet = index of the selected leading jet (-1: none);
//      active-area ghosts of those jets: ghostEta, ghostPhi, ghostJet;
//      CS ghost grid: csGhostEta, csGhostPhi.
//    tree "recoil" (-recoil only), one entry per (event, dRmax = 1.0, 0.5):
//      finPt, finY, finPhi, finPtAfter, finCharged (final-state particles,
//      dummies excluded), thPt, thY, thPhi, thPtAfter (scattering centres,
//      status 3), nPairs (pairs with dR <= dRmax).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVector2.h"

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"

#ifdef WITH_PYTHIA8
#include "Pythia8/Pythia.h"
#endif

using namespace std;
using fastjet::PseudoJet;

// ===========================================================================
//  Configuration
// ===========================================================================
static constexpr int kNumJetRadii = 4, kNumAcceptances = 2, kNumProbeBins = 25;
static const double kJetRadii[kNumJetRadii]            = {0.1, 0.2, 0.4, 0.8};
static const double kTrackEtaMax[kNumAcceptances]      = {0.9, 2.5};   // ALICE, wide
static const double kProbePtEdges[kNumProbeBins + 1]   = {0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0,
                                                          5.0, 6.0, 8.0, 10., 12., 15., 20., 25., 30., 40.,
                                                          50., 60., 80., 100., 150., 300.};
static constexpr int kNumAccR = kNumJetRadii * kNumAcceptances;        // k = iacc * kNumJetRadii + iR

static const double kTrackPtMin = 0.15;                    // GeV/c, all tracks
static const double kRefPtMin = 0.2, kRefPtMax = 3.0;      // GeV/c, soft reference window [min, max)
static const double kEtaGap = 0.5;                         // eta-gap reference: |eta - eta_jet| > R + kEtaGap
static const double kAwayHalfWidth = 1.0;                  // refAw: drop |dphi - pi| < 1
static const double kJetGhostArea = 0.05;                  // active area of the raw jets (rho * A)
static const double kRhoConeR = 0.4;                       // perp-cone R whose rho drives the CS
// event-wise constituent subtraction (fjcontrib ConstituentSubtractor)
static const double kCSMaxDeltaR = 0.25, kCSAlpha = 0.0, kCSGhostArea = 0.01;
// recoil thermal-momentum subtraction (variants cs10, cs05)
static const double kRecoilSubDRmaxLarge = 1.0, kRecoilSubDRmaxSmall = 0.5;
static const double kDummyEnergyMax = 1e-5;                // JEWEL dummy particles carry e < 1e-5 GeV

// thermal background: uncorrelated charged pions with Gamma(2,T) pT spectrum
struct ThermalConfig { string tag; double dndeta; double meanPt; };
static const vector<ThermalConfig> kThermalConfigs = {
  {"embMB", 42.0,  0.70},   // O+O 5.36 TeV minimum bias (CMS: 41.8)
  {"embC",  135.0, 0.70}};  // O+O 5.36 TeV 0-5%        (CMS: 135)
// meanPt is the mean of the full Gamma(2, T) spectrum, i.e. T = meanPt / 2 = 0.35 GeV
static const double kThermalEtaMax = 3.5;
static const double kFwdPtMin = 0.2, kFwdPtMax = 10.0;   // forward-estimator pT window (GeV/c)
static const double kThermalDefaultSlope = 0.35;           // GeV, used by embNp (same as embMB/embC)
static const double kPionMass = 0.13957;                   // GeV/c^2

// ===========================================================================
//  Particles and kinematics
// ===========================================================================
struct Particle { double px, py, pz, e; int pid, status; };

// three times the electric charge from the PDG code (HepPID convention)
static int threeTimesCharge(int pid)
{
  static const int quarkCharge3[7] = {0, -1, 2, -1, 2, -1, 2};   // -, d u s c b t
  int a = abs(pid), sign = pid > 0 ? 1 : -1;
  if (a <= 6) return sign * quarkCharge3[a];
  if (a == 11 || a == 13 || a == 15) return -3 * sign;
  if (a == 12 || a == 14 || a == 16 || a == 21 || a == 22 || a == 23) return 0;
  if (a == 24) return 3 * sign;
  if (a < 100 || a > 9999999) return 0;
  int q1 = (a / 1000) % 10, q2 = (a / 100) % 10, q3 = (a / 10) % 10;
  if (q2 > 6 || q3 > 6 || q1 > 6) return 0;
  int c;
  if (q1 == 0) c = (q2 == 3 || q2 == 5) ? quarkCharge3[q3] - quarkCharge3[q2] : quarkCharge3[q2] - quarkCharge3[q3];  // meson
  else         c = quarkCharge3[q1] + quarkCharge3[q2] + quarkCharge3[q3];                                           // baryon
  return sign * c;
}

static double transverseMomentum(const Particle &p) { return hypot(p.px, p.py); }
static double pseudorapidity(double px, double py, double pz)
{
  double pt = hypot(px, py), p = sqrt(pt * pt + pz * pz);
  if (pt < 1e-12) return pz > 0 ? 1e9 : -1e9;
  return 0.5 * log((p + pz) / (p - pz));
}
// |phi_a - phi_b| folded into [0, pi]
static double deltaPhi(double a, double b)
{
  double d = fabs(a - b);
  while (d > 2 * M_PI) d -= 2 * M_PI;
  return d > M_PI ? 2 * M_PI - d : d;
}

// charged tracks (status 1, charge != 0, pT > 0.15 GeV/c, |eta| < etaMax) as
// massless PseudoJets; user_index = position in the returned vector.
// sourceIndex (event display only): position of each track in event.
static vector<PseudoJet> selectChargedTracks(const vector<Particle> &event, double etaMax,
                                             vector<int> *sourceIndex = nullptr)
{
  vector<PseudoJet> tracks;
  if (sourceIndex) sourceIndex->clear();
  for (size_t iSource = 0; iSource < event.size(); iSource++) {
    const Particle &p = event[iSource];
    if (p.status != 1 || threeTimesCharge(p.pid) == 0) continue;
    double pt = transverseMomentum(p);
    if (pt < kTrackPtMin) continue;
    double eta = pseudorapidity(p.px, p.py, p.pz);
    if (fabs(eta) >= etaMax) continue;
    PseudoJet track(p.px, p.py, p.pz, sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz));
    track.set_user_index((int)tracks.size());
    tracks.push_back(track);
    if (sourceIndex) sourceIndex->push_back((int)iSource);
  }
  return tracks;
}

// ===========================================================================
//  Particle lists of the recoil sample
// ===========================================================================
// Status-1 particles without the JEWEL dummies (variant nom / nosub).
static vector<Particle> finalStateParticles(const vector<Particle> &event)
{
  vector<Particle> out;
  for (auto &p : event) if (p.status == 1 && p.e > kDummyEnergyMax) out.push_back(p);
  return out;
}

// Recoil thermal-momentum subtraction (variants cs10, cs05). Port of
// jewel-2.6.0/subtraction.C: every pair (final-state particle, thermal scattering
// centre = status 3) with dR = sqrt(dphi^2 + dy^2) <= dRmax is visited in order
// of increasing dR; the smaller of the two pT (and of m_delta = mT - pT) is
// removed from both. Particles left with pT > 0 are kept (massive, rebuilt from
// pT, m_delta, phi, y); unsubtracted thermal remainders are dropped.
// Event display only: the inputs and the pT left after the subtraction.
struct RecoilSubtractionRecord {
  double dRmax = 0; long nPairs = 0;
  vector<float> finalPt, finalY, finalPhi, finalPtAfter; vector<int> finalCharged;
  vector<float> thermalPt, thermalY, thermalPhi, thermalPtAfter;
  vector<int> keptFinal;              // per output particle: index into final*
  // origin of every output particle: 2 if its pT was reduced, else 0 (generator)
  vector<int> outputOrigin() const
  {
    vector<int> origin;
    for (int i : keptFinal) origin.push_back(finalPtAfter[i] != finalPt[i] ? 2 : 0);
    return origin;
  }
};
static vector<Particle> subtractRecoilThermalMomenta(const vector<Particle> &event, double dRmax,
                                                     RecoilSubtractionRecord *record = nullptr)
{
  struct Kin { double pt, mDelta, phi, y; int index; };
  vector<Kin> finals, thermals;
  for (size_t i = 0; i < event.size(); i++) {
    const Particle &p = event[i];
    if (p.e - fabs(p.pz) <= 0.) continue;                    // rapidity undefined
    double pt = transverseMomentum(p), m2 = p.e * p.e - pt * pt - p.pz * p.pz;
    Kin q{pt, sqrt(max(0., m2) + pt * pt) - pt, atan2(p.py, p.px),
          0.5 * log((p.e + p.pz) / (p.e - p.pz)), (int)i};
    if (p.status == 1 && p.e > kDummyEnergyMax) finals.push_back(q);   // skip JEWEL dummies
    else if (p.status == 3)                      thermals.push_back(q);
  }
  struct Pair { double dR; int iFinal, iThermal; };
  vector<Pair> pairs;
  pairs.reserve(finals.size() * thermals.size() / 4 + 16);
  for (int i = 0; i < (int)finals.size(); i++)
    for (int j = 0; j < (int)thermals.size(); j++) {
      double dR = hypot(deltaPhi(finals[i].phi, thermals[j].phi), finals[i].y - thermals[j].y);
      if (dR <= dRmax) pairs.push_back({dR, i, j});
    }
  sort(pairs.begin(), pairs.end(), [](const Pair &a, const Pair &b) { return a.dR < b.dR; });
  if (record) {
    *record = RecoilSubtractionRecord(); record->dRmax = dRmax; record->nPairs = (long)pairs.size();
    for (auto &f : finals) {
      record->finalPt.push_back(f.pt); record->finalY.push_back(f.y); record->finalPhi.push_back(f.phi);
      record->finalCharged.push_back(threeTimesCharge(event[f.index].pid) != 0);
    }
    for (auto &t : thermals) { record->thermalPt.push_back(t.pt); record->thermalY.push_back(t.y); record->thermalPhi.push_back(t.phi); }
  }
  for (auto &pr : pairs) {
    Kin &f = finals[pr.iFinal], &t = thermals[pr.iThermal];
    if (f.pt > t.pt)         { f.pt -= t.pt; t.pt = 0.; }         else { t.pt -= f.pt; f.pt = 0.; }
    if (f.mDelta > t.mDelta) { f.mDelta -= t.mDelta; t.mDelta = 0.; } else { t.mDelta -= f.mDelta; f.mDelta = 0.; }
  }
  if (record) {
    for (auto &f : finals) record->finalPtAfter.push_back(f.pt);
    for (auto &t : thermals) record->thermalPtAfter.push_back(t.pt);
  }
  vector<Particle> out;
  for (size_t iFinal = 0; iFinal < finals.size(); iFinal++) {
    const Kin &f = finals[iFinal];
    if (f.pt <= 0.) continue;
    if (record) record->keptFinal.push_back((int)iFinal);
    double mt = f.pt + f.mDelta;
    out.push_back({f.pt * cos(f.phi), f.pt * sin(f.phi), mt * sinh(f.y), mt * cosh(f.y),
                   event[f.index].pid, 1});
  }
  return out;
}

// ===========================================================================
//  Thermal underlying event
// ===========================================================================
// Appends uncorrelated charged pions to a copy of the event:
//   N ~ Poisson(dndeta * 2 * 3.5), |eta| < 3.5 uniform, phi uniform,
//   pT ~ Gamma(2, T) = pT exp(-pT/T) (<pT> = 2T), charge +-1 with equal probability.
// Random-number order per pion: 2 x pT, eta, phi, charge.
static vector<Particle> addThermalPions(const vector<Particle> &event, double dndeta, double T, TRandom3 &rng)
{
  vector<Particle> out = event;
  int n = rng.Poisson(dndeta * 2 * kThermalEtaMax);
  for (int i = 0; i < n; i++) {
    double pt  = -T * log(rng.Rndm() * rng.Rndm());                   // Gamma(2,T)
    double eta = (2 * rng.Rndm() - 1) * kThermalEtaMax, phi = 2 * M_PI * rng.Rndm();
    double px = pt * cos(phi), py = pt * sin(phi), pz = pt * sinh(eta), m = kPionMass;
    out.push_back({px, py, pz, sqrt(px * px + py * py + pz * pz + m * m),
                   rng.Rndm() < 0.5 ? 211 : -211, 1});
  }
  return out;
}

// Mean pT of Gamma(2,T) restricted to the window [lo, hi) (midpoint rule).
static double gammaWindowMeanPt(double T, double lo, double hi)
{
  double s0 = 0, s1 = 0; int n = 4000;
  for (int i = 0; i < n; i++) { double p = lo + (hi - lo) * (i + 0.5) / n, f = p * exp(-p / T); s0 += f; s1 += f * p; }
  return s1 / s0;
}
// Slope T with gammaWindowMeanPt(T, lo, hi) = target (bisection in [0.05, 3] GeV).
static double slopeForWindowMeanPt(double target, double lo, double hi)
{
  double a = 0.05, b = 3.0;
  for (int it = 0; it < 100; it++) { double m = 0.5 * (a + b); (gammaWindowMeanPt(m, lo, hi) < target ? a : b) = m; }
  return 0.5 * (a + b);
}
// Tabulated inverse of gammaWindowMeanPt on the reference window, T in
// [0.10, 1.20] GeV in steps of 1 MeV, linear interpolation (per-event use).
struct WindowMeanPtInverse {
  vector<double> slope, windowMean;
  void build()
  {
    for (double t = 0.10; t <= 1.2001; t += 0.001) { slope.push_back(t); windowMean.push_back(gammaWindowMeanPt(t, kRefPtMin, kRefPtMax)); }
  }
  double slopeFor(double target) const
  {
    if (target <= windowMean.front()) return slope.front();
    if (target >= windowMean.back()) return slope.back();
    size_t i = upper_bound(windowMean.begin(), windowMean.end(), target) - windowMean.begin();
    double f = (target - windowMean[i - 1]) / (windowMean[i] - windowMean[i - 1]);
    return slope[i - 1] * (1 - f) + slope[i] * f;
  }
};

// ---------------------------------------------------------------------------
//  Optical Glauber for 16O+16O: N_part vs centrality percentile.
//  Harmonic-oscillator density rho(r) ~ (1 + alpha r^2/a^2) exp(-r^2/a^2),
//  a = 1.833 fm, alpha = 1.544; sigma_NN = 6.8 fm^2 (68 mb).
//  Multiplicity model dN/deta = A * N_part^beta, with A and beta fixed by the CMS
//  values <dN/deta>(0-5%) = 135 and <dN/deta>(0-100%) = 41.8 (centrality-class
//  averages = averages over uniformly distributed centrality).
// ---------------------------------------------------------------------------
struct OpticalGlauberOxygen {
  vector<double> centrality, npart;   // monotonic table, centrality in percent
  double A = 0, beta = 1;
  static double density(double r) { double a = 1.833, al = 1.544; return (1 + al * r * r / (a * a)) * exp(-r * r / (a * a)); }
  void build()
  {
    // nuclear thickness T_A(s) on a grid in s, normalised to A = 16
    const int NS = 400; const double SMAX = 8.0, sigma = 6.8, Anuc = 16;
    vector<double> thick(NS + 1);
    for (int i = 0; i <= NS; i++) {
      double s = SMAX * i / NS, sum = 0; int nz = 800; double zmax = 10;
      for (int j = 0; j < nz; j++) { double z = -zmax + 2 * zmax * (j + 0.5) / nz; sum += density(sqrt(s * s + z * z)); }
      thick[i] = sum * 2 * zmax / nz;
    }
    double norm = 0; for (int i = 0; i < NS; i++) { double s = SMAX * (i + 0.5) / NS; norm += 2 * M_PI * s * (thick[i] + thick[i + 1]) / 2 * SMAX / NS; }
    for (auto &t : thick) t *= Anuc / norm;
    auto TA = [&](double s) { double x = s / SMAX * NS; if (x >= NS) return 0.; int i = (int)x; double f = x - i; return thick[i] * (1 - f) + thick[i + 1] * f; };
    // N_part(b) and dsigma/db on an impact-parameter grid, transverse grid h = 0.1 fm
    const int NBIMP = 160; const double BMAX = 10.0; const double h = 0.1; const int NG = 140;
    vector<double> bb(NBIMP + 1), np(NBIMP + 1), dsig(NBIMP + 1);
    for (int ib = 0; ib <= NBIMP; ib++) {
      double b = BMAX * ib / NBIMP, Np = 0, TAB = 0;
      for (int ix = 0; ix < NG; ix++) for (int iy = 0; iy < NG; iy++) {
        double x = -7 + h * (ix + 0.5), y = -7 + h * (iy + 0.5);
        double ta = TA(hypot(x + b / 2, y)), tb = TA(hypot(x - b / 2, y));
        Np += ta * (1 - pow(1 - sigma * tb / Anuc, Anuc)) + tb * (1 - pow(1 - sigma * ta / Anuc, Anuc));
        TAB += ta * tb;
      }
      Np *= h * h; TAB *= h * h;
      bb[ib] = b; np[ib] = Np; dsig[ib] = 2 * M_PI * b * (1 - pow(1 - sigma * TAB / (Anuc * Anuc), Anuc * Anuc));
    }
    // centrality percentile = cumulative inelastic cross section
    vector<double> cum(NBIMP + 1, 0);
    for (int ib = 1; ib <= NBIMP; ib++) cum[ib] = cum[ib - 1] + 0.5 * (dsig[ib] + dsig[ib - 1]) * BMAX / NBIMP;
    for (int ib = 0; ib <= NBIMP; ib++) { centrality.push_back(100. * cum[ib] / cum[NBIMP]); npart.push_back(np[ib]); }
    // fit A, beta to the two CMS constraints (bisection in beta on the 0-5% / 0-100% ratio)
    auto avg = [&](double c0, double c1, double be) {
      double s = 0; int n = 2000;
      for (int i = 0; i < n; i++) s += pow(max(npartAt(c0 + (c1 - c0) * (i + 0.5) / n), 1e-6), be);
      return s / n;
    };
    double lo = 0.3, hi = 2.5, target = 135. / 41.8;
    for (int it = 0; it < 80; it++) { double m = 0.5 * (lo + hi); (avg(0, 5, m) / avg(0, 100, m) < target ? lo : hi) = m; }
    beta = 0.5 * (lo + hi); A = 41.8 / avg(0, 100, beta);
  }
  double npartAt(double c) const
  {
    if (c <= centrality.front()) return npart.front();
    if (c >= centrality.back()) return npart.back();
    size_t i = upper_bound(centrality.begin(), centrality.end(), c) - centrality.begin();
    double f = (c - centrality[i - 1]) / (centrality[i] - centrality[i - 1] + 1e-12);
    return npart[i - 1] * (1 - f) + npart[i] * f;
  }
  double dndetaAt(double c) const { return A * pow(max(npartAt(c), 1e-6), beta); }
  double classMean(double c0, double c1, bool multiplicity) const
  { double s = 0; int n = 2000; for (int i = 0; i < n; i++) { double c = c0 + (c1 - c0) * (i + 0.5) / n; s += multiplicity ? dndetaAt(c) : npartAt(c); } return s / n; }
  double classMeanLogDndeta(double c0, double c1) const
  { double s = 0; int n = 2000; for (int i = 0; i < n; i++) s += log(dndetaAt(c0 + (c1 - c0) * (i + 0.5) / n)); return s / n; }
};

// ---------------------------------------------------------------------------
//  Flowing thermal UE (variant embFlow). For centrality c:
//    m(c) = M_MB + K * (ln dN/deta(c) - <ln dN/deta>_0-100)
//           window <pT> (0.2-3 GeV/c); M_MB = 0.767 GeV/c (= embMB);
//           K fixed by m(0-5%) - m(70-80%) = 0.11 GeV/c
//  Per event one z ~ N(0,1), shared by multiplicity and slope:
//    dN/deta = max(0, dN/deta(c) * (1 + z / sqrt(max(N_part(c), 2))))
//    T       = max(0.05, T_window^-1(m(c)) * (1 + 0.05 z))
//  c is the JEWEL medium centrality (O+O medium), so quenching and UE
//  co-fluctuate; vacuum/PYTHIA use the same drawn c as embNp.
// ---------------------------------------------------------------------------
struct FlowingUEModel {
  static constexpr double M_MB = 0.767, SPREAD = 0.11, TFLUC = 0.05;
  double K = 0, L0 = 0;
  WindowMeanPtInverse inverse;
  void build(const OpticalGlauberOxygen &gl)
  {
    inverse.build();
    L0 = gl.classMeanLogDndeta(0, 100);
    K = SPREAD / (gl.classMeanLogDndeta(0, 5) - gl.classMeanLogDndeta(70, 80));
  }
  double windowMeanPt(const OpticalGlauberOxygen &gl, double c) const { return M_MB + K * (log(gl.dndetaAt(c)) - L0); }
  double classWindowMeanPt(const OpticalGlauberOxygen &gl, double c0, double c1) const
  { double s = 0; int n = 2000; for (int i = 0; i < n; i++) s += windowMeanPt(gl, c0 + (c1 - c0) * (i + 0.5) / n); return s / n; }
};

// ---------------------------------------------------------------------------
//  Generator UE events (Angantyr / HIJING) from CompactEvent v1h files (tree
//  "events"; pt, eta stored x6000, phi x 65535/2pi). Events are taken in file
//  order, starting at a seed-dependent entry and cycling through the chain;
//  only events with LO <= N_ch(|eta| < 0.5) <= HI are used. All stored particles
//  (charged and neutral) are embedded with the pion mass.
// ---------------------------------------------------------------------------
struct UEEventSource {
  string tag; double lo, hi;
  TChain *chain = nullptr; Long64_t pos = 0, nEntries = 0;
  static const int NMAX = 40000;
  Int_t n; UShort_t pt[NMAX], phi[NMAX]; Short_t eta[NMAX]; Int_t pdg[NMAX]; Char_t charge[NMAX];
  long scanned = 0, used = 0;
  bool open(const string &listFile, unsigned seed)
  {
    chain = new TChain("events");
    ifstream in(listFile); string l; while (getline(in, l)) if (!l.empty() && l[0] != '#') chain->Add(l.c_str());
    nEntries = chain->GetEntries(); if (nEntries <= 0) return false;
    chain->SetBranchStatus("*", 0);
    for (auto b : {"n", "pt", "eta", "phi", "pdg", "charge"}) chain->SetBranchStatus(b, 1);
    chain->SetBranchAddress("n", &n); chain->SetBranchAddress("pt", pt); chain->SetBranchAddress("eta", eta);
    chain->SetBranchAddress("phi", phi); chain->SetBranchAddress("pdg", pdg); chain->SetBranchAddress("charge", charge);
    pos = (Long64_t)((seed * 2654435761ULL) % (unsigned long long)nEntries);   // Knuth hash of the seed
    return true;
  }
  // next UE event inside the multiplicity window; false if a full cycle finds none
  bool next(vector<Particle> &out)
  {
    for (long tries = 0; tries < nEntries; tries++) {
      chain->GetEntry(pos); pos = (pos + 1) % nEntries; scanned++;
      if (n > NMAX) continue;
      int nChargedMid = 0;
      for (int i = 0; i < n; i++) if (charge[i] != 0 && fabs(eta[i] / 6000.) < 0.5) nChargedMid++;
      if (nChargedMid < lo || nChargedMid > hi) continue;
      out.clear();
      for (int i = 0; i < n; i++) {
        double p = pt[i] / 6000., e = eta[i] / 6000., f = phi[i] / (65535. / (2 * M_PI)), m = kPionMass;
        double px = p * cos(f), py = p * sin(f), pz = p * sinh(e);
        out.push_back({px, py, pz, sqrt(px * px + py * py + pz * pz + m * m), pdg[i], 1});
      }
      used++;
      return true;
    }
    return false;
  }
};

// ===========================================================================
//  HepMC2 ASCII reader (IO_GenEvent), enough for JEWEL output
// ===========================================================================
class HepMC2Reader {
public:
  explicit HepMC2Reader(const string &f) : in(f) {}
  bool ok() const { return in.good(); }
  // Reads the next event; returns false at end of file.
  //   weight     : first entry of the E-record weight list (1 if none)
  //   centrality : 10th number of the H (heavy-ion) record, where the patched
  //                JEWEL writes the sampled centrality as a fraction; converted
  //                to percent. Vacuum JEWEL writes -1; -1 if no H record.
  bool next(vector<Particle> &event, double &weight, double &centrality)
  {
    event.clear(); weight = 1.; centrality = -1.;
    string line;
    if (pendingEventLine.empty()) {
      while (getline(in, line)) if (!line.empty() && line[0] == 'E') { pendingEventLine = line; break; }
      if (pendingEventLine.empty()) return false;
    }
    parseEventLine(pendingEventLine, weight);
    pendingEventLine.clear();
    while (getline(in, line)) {
      if (line.empty()) continue;
      char type = line[0];
      if (type == 'E') { pendingEventLine = line; return true; }
      if (type == 'P') {
        // P barcode pid px py pz e m status ...
        istringstream s(line); string tag; int barcode, pid, status; double px, py, pz, e, m;
        s >> tag >> barcode >> pid >> px >> py >> pz >> e >> m >> status;
        event.push_back({px, py, pz, e, pid, status});
      } else if (type == 'H') {
        istringstream s(line); string tag; double x;
        vector<double> v; s >> tag; while (s >> x) v.push_back(x);
        if (v.size() > 9) centrality = v[9];                         // impact_parameter slot
        if (centrality > 0 && centrality <= 1.0) centrality *= 100.; // fraction -> percent
      }
    }
    return true;                                                     // last event in file
  }
private:
  static void parseEventLine(const string &l, double &weight)
  {
    istringstream s(l); string tag; vector<double> v; double x;
    s >> tag; while (s >> x) v.push_back(x);
    // E evnum nmpi scale aqcd aqed sigid bcvtx nvtx bp1 bp2 nrand [r..] nw [w..]
    if (v.size() < 11) return;
    size_t i = 10, nRandom = (size_t)v[i]; i += 1 + nRandom;
    if (i < v.size()) { size_t nWeights = (size_t)v[i]; if (nWeights > 0 && i + 1 < v.size()) weight = v[i + 1]; }
  }
  ifstream in;
  string pendingEventLine;
};

// ===========================================================================
//  Per-variant output tree (branches documented in the header)
// ===========================================================================
struct VariantTree {
  string name;
  TTree *tree = nullptr;
  float weight, centr, ueMult = 0; int nchMid, nchFwdA, nchFwdC;
  float rhoEv[kNumAcceptances];
  // per k = iacc * kNumJetRadii + iR
  float jetPt[kNumAccR], jetPtCorr[kNumAccR], jetEta[kNumAccR], jetPhi[kNumAccR], jetArea[kNumAccR], leadTrk[kNumAccR], rho[kNumAccR];
  float csJetPt[kNumAccR], csJetEta[kNumAccR], csJetPhi[kNumAccR], csLeadTrk[kNumAccR];
  float nIn[kNumAccR * kNumProbeBins], nPerp[kNumAccR * kNumProbeBins], nCS[kNumAccR * kNumProbeBins];
  float refN[kNumAccR], refS1[kNumAccR], refS2[kNumAccR];
  float refNa[kNumAccR], refS1a[kNumAccR], refS2a[kNumAccR], refNb[kNumAccR], refS1b[kNumAccR], refS2b[kNumAccR];
  float gapN[kNumAccR], gapS1[kNumAccR], gapS2[kNumAccR];
  float refAwN[kNumAccR], refAwS1[kNumAccR], refAwS2[kNumAccR];
  float refAwNa[kNumAccR], refAwS1a[kNumAccR], refAwS2a[kNumAccR], refAwNb[kNumAccR], refAwS1b[kNumAccR], refAwS2b[kNumAccR];
  float nearN[kNumAccR], nearS1[kNumAccR], nearS2[kNumAccR];
  float transN[kNumAccR], transS1[kNumAccR], transS2[kNumAccR];
  float awayN[kNumAccR], awayS1[kNumAccR], awayS2[kNumAccR];
  float jetPtTruth[kNumAccR], dRTruth[kNumAccR];

  // Branch order and leaf lists are part of the file format: do not reorder.
  void book(TFile *f)
  {
    f->cd();
    tree = new TTree(("t_" + name).c_str(), ("jet-UE v0 building blocks, variant " + name).c_str());
    string perK = "[" + to_string(kNumAccR) + "]/F", perKBin = "[" + to_string(kNumAccR * kNumProbeBins) + "]/F";
    tree->Branch("weight", &weight, "weight/F");
    tree->Branch("centr", &centr, "centr/F");
    tree->Branch("nchMid", &nchMid, "nchMid/I");
    tree->Branch("nchFwdA", &nchFwdA, "nchFwdA/I");
    tree->Branch("nchFwdC", &nchFwdC, "nchFwdC/I");
    tree->Branch("rhoEv", rhoEv, "rhoEv[2]/F");
    tree->Branch("ueMult", &ueMult, "ueMult/F");
    auto br = [&](const char *n, float *a, const string &s) { tree->Branch(n, a, (string(n) + s).c_str()); };
    br("jetPt", jetPt, perK); br("jetPtCorr", jetPtCorr, perK); br("jetEta", jetEta, perK); br("jetPhi", jetPhi, perK);
    br("jetArea", jetArea, perK); br("leadTrk", leadTrk, perK); br("rho", rho, perK);
    br("csJetPt", csJetPt, perK); br("csJetEta", csJetEta, perK); br("csJetPhi", csJetPhi, perK); br("csLeadTrk", csLeadTrk, perK);
    br("nIn", nIn, perKBin); br("nPerp", nPerp, perKBin); br("nCS", nCS, perKBin);
    br("refN", refN, perK); br("refS1", refS1, perK); br("refS2", refS2, perK);
    br("refNa", refNa, perK); br("refS1a", refS1a, perK); br("refS2a", refS2a, perK);
    br("refNb", refNb, perK); br("refS1b", refS1b, perK); br("refS2b", refS2b, perK);
    br("gapN", gapN, perK); br("gapS1", gapS1, perK); br("gapS2", gapS2, perK);
    br("refAwN", refAwN, perK); br("refAwS1", refAwS1, perK); br("refAwS2", refAwS2, perK);
    br("refAwNa", refAwNa, perK); br("refAwS1a", refAwS1a, perK); br("refAwS2a", refAwS2a, perK);
    br("refAwNb", refAwNb, perK); br("refAwS1b", refAwS1b, perK); br("refAwS2b", refAwS2b, perK);
    br("nearN", nearN, perK); br("nearS1", nearS1, perK); br("nearS2", nearS2, perK);
    br("transN", transN, perK); br("transS1", transS1, perK); br("transS2", transS2, perK);
    br("awayN", awayN, perK); br("awayS1", awayS1, perK); br("awayS2", awayS2, perK);
    br("jetPtTruth", jetPtTruth, perK); br("dRTruth", dRTruth, perK);
  }
  // zero all per-event jet quantities (weight, centr, ueMult and the
  // multiplicities are set explicitly for every event)
  void reset()
  {
    rhoEv[0] = rhoEv[1] = 0;
    for (int i = 0; i < kNumAccR; i++) {
      jetPt[i] = jetPtCorr[i] = jetEta[i] = jetPhi[i] = jetArea[i] = leadTrk[i] = rho[i] = 0;
      csJetPt[i] = csJetEta[i] = csJetPhi[i] = csLeadTrk[i] = 0;
      refN[i] = refS1[i] = refS2[i] = refNa[i] = refS1a[i] = refS2a[i] = refNb[i] = refS1b[i] = refS2b[i] = 0;
      gapN[i] = gapS1[i] = gapS2[i] = 0;
      refAwN[i] = refAwS1[i] = refAwS2[i] = refAwNa[i] = refAwS1a[i] = refAwS2a[i] = refAwNb[i] = refAwS1b[i] = refAwS2b[i] = 0;
      nearN[i] = nearS1[i] = nearS2[i] = transN[i] = transS1[i] = transS2[i] = awayN[i] = awayS1[i] = awayS2[i] = 0;
      jetPtTruth[i] = 0; dRTruth[i] = -1;
    }
    for (int i = 0; i < kNumAccR * kNumProbeBins; i++) nIn[i] = nPerp[i] = nCS[i] = 0;
  }
};

// count, sum pT and sum pT^2 of one track into the three reference sums at k
static void addToMoments(float *N, float *S1, float *S2, int k, double pt)
{
  N[k] += 1; S1[k] += pt; S2[k] += pt * pt;
}

// probe-pT bin index in kProbePtEdges, -1 outside [0.15, 300)
static int probePtBin(double pt)
{
  if (pt < kProbePtEdges[0] || pt >= kProbePtEdges[kNumProbeBins]) return -1;
  return int(upper_bound(kProbePtEdges, kProbePtEdges + kNumProbeBins + 1, pt) - kProbePtEdges) - 1;
}

// (eta, phi) distance to a cone axis strictly below R
static bool insideCone(const PseudoJet &p, double eta, double phi, double R)
{
  double de = p.eta() - eta, dp = deltaPhi(p.phi(), phi);
  return de * de + dp * dp < R * R;
}

// ===========================================================================
//  Event display records (-display; see header). Filled by the production
//  functions below through optional pointers; they only read existing objects.
// ===========================================================================
enum DisplayTrackFlag {
  kDisplayFlagLeadJet = 1, kDisplayFlagJetCone = 2, kDisplayFlagPerpCone = 4, kDisplayFlagRef = 8,
  kDisplayFlagRefAway = 16, kDisplayFlagGap = 32, kDisplayFlagCSJet = 64
};
static const double kDisplayJetPtMin = 3.0;                           // GeV/c, jets written to the display
static const int kDisplaySlots[] = {0 * kNumJetRadii + 2, 1 * kNumJetRadii + 3};   // (0.9, R 0.4), (2.5, R 0.8)

// one (variant, acceptance, R) slot of one event
struct DisplaySlotRecord {
  vector<float> trkPt, trkEta, trkPhi, trkPtCS; vector<int> trkOrigin, trkJet, trkFlags;
  vector<float> jetPt, jetEta, jetPhi, jetArea; vector<int> jetNConst, jetFiducial;
  vector<float> ghostEta, ghostPhi; vector<int> ghostJet;
  vector<float> csGhostEta, csGhostPhi;
  int leadJet = -1;
  void clear() { *this = DisplaySlotRecord(); }
  // tracks entering the clustering; origin of each source particle (empty: 0)
  void recordTracks(const vector<PseudoJet> &tracks, const vector<int> &sourceIndex, const vector<int> &sourceOrigin)
  {
    for (size_t i = 0; i < tracks.size(); i++) {
      trkPt.push_back(tracks[i].pt()); trkEta.push_back(tracks[i].eta()); trkPhi.push_back(tracks[i].phi_std());
      int s = sourceIndex[i];
      trkOrigin.push_back(s < (int)sourceOrigin.size() ? sourceOrigin[s] : 0);
      trkJet.push_back(-1); trkFlags.push_back(0); trkPtCS.push_back(0.f);
    }
  }
  // all jets of the cluster sequence above kDisplayJetPtMin (and the selected
  // leading jet), their constituents and their active-area ghosts
  void recordJets(const fastjet::ClusterSequence &clusterSeq, const PseudoJet *leading, double jetEtaMax)
  {
    double ptMin = leading ? min(kDisplayJetPtMin, leading->pt()) * (1 - 1e-9) : kDisplayJetPtMin;
    vector<PseudoJet> jets = fastjet::sorted_by_pt((!fastjet::SelectorIsPureGhost())(clusterSeq.inclusive_jets(ptMin)));
    for (auto &j : jets) {
      int index = (int)jetPt.size(), nConst = 0;
      if (leading && j.cluster_hist_index() == leading->cluster_hist_index()) leadJet = index;
      for (auto &c : j.constituents()) {
        if (c.is_pure_ghost()) { ghostEta.push_back(c.eta()); ghostPhi.push_back(c.phi_std()); ghostJet.push_back(index); }
        else { trkJet[c.user_index()] = index; nConst++; }
      }
      jetPt.push_back(j.pt()); jetEta.push_back(j.eta()); jetPhi.push_back(j.phi_std());
      jetArea.push_back(j.has_area() ? j.area() : 0.); jetNConst.push_back(nConst);
      jetFiducial.push_back(fabs(j.eta()) < jetEtaMax);
    }
  }
};
// display state of one variant in one event
struct VariantDisplay {
  vector<int> origin;                                      // per particle of the variant's list
  DisplaySlotRecord slot[kNumAccR];
  bool wants(int k) const { for (int s : kDisplaySlots) if (s == k) return true; return false; }
};

// ===========================================================================
//  Jets and background
// ===========================================================================
// Leading anti-kT jet (inclusive jets with pT > 1 GeV/c) with |eta_jet| < etaMax - R;
// false if none. withArea: active area with explicit ghosts up to |eta| = etaMax
// (ghost area 0.05; pure-ghost jets removed); otherwise no ghosts and the area is
// set to pi R^2. The ghosts draw from FastJet's global random generator.
static bool findLeadingJet(const vector<PseudoJet> &tracks, double R, double etaMax, bool withArea,
                           PseudoJet &jet, double &area, unique_ptr<fastjet::ClusterSequence> &clusterSeq)
{
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);
  if (withArea) {
    fastjet::GhostedAreaSpec ghosts(etaMax, 1, kJetGhostArea);
    fastjet::AreaDefinition areaDef(fastjet::active_area_explicit_ghosts, ghosts);
    clusterSeq.reset(new fastjet::ClusterSequenceArea(tracks, jetDef, areaDef));
  } else {
    clusterSeq.reset(new fastjet::ClusterSequence(tracks, jetDef));
  }
  fastjet::Selector select = fastjet::SelectorAbsEtaMax(etaMax - R);   // jet fully inside the acceptance
  if (withArea) select = select * !fastjet::SelectorIsPureGhost();
  auto jets = fastjet::sorted_by_pt(select(clusterSeq->inclusive_jets(1.0)));
  if (jets.empty()) return false;
  jet = jets[0]; area = withArea ? jet.area() : M_PI * R * R;
  return true;
}

// Perpendicular-cone background: two cones of radius R at (eta_jet, phi_jet +- pi/2).
// Tracks inside the jet cone are never counted. Marks the cone tracks in
// inPerpCone, adds 0.5 per track to perpCounts[probe bin] (average of the two
// cones) and returns rho = sum pT / (2 pi R^2).
static double fillPerpendicularCones(const vector<PseudoJet> &tracks, double eta, double phi, double R,
                                     float *perpCounts, vector<char> &inPerpCone)
{
  double sumPt = 0;
  for (int side = -1; side <= 1; side += 2) {
    double conePhi = phi + side * M_PI / 2;
    for (size_t i = 0; i < tracks.size(); i++) {
      if (insideCone(tracks[i], eta, phi, R) || !insideCone(tracks[i], eta, conePhi, R)) continue;
      inPerpCone[i] = 1;
      sumPt += tracks[i].pt();
      if (perpCounts) { int b = probePtBin(tracks[i].pt()); if (b >= 0) perpCounts[b] += 0.5f; }
    }
  }
  return sumPt / (2 * M_PI * R * R);
}

// Reference-region sums for one (acceptance, R) slot k. Soft tracks
// 0.2 <= pT < 3 GeV/c that are neither jet constituents nor within R of the jet
// axis enter:
//   - if outside both perpendicular cones: ref (full, a: eta<0, b: eta>=0),
//     refAw if |dphi| <= pi - 1 (away side |dphi - pi| < 1 excluded),
//     and exactly one of near (|dphi| < pi/3), trans (< 2pi/3), away (>= 2pi/3);
//   - independently of the perpendicular cones: gap if |eta - eta_jet| > R + 0.5.
// displayFlags (event display only): per track, bits kDisplayFlagRef/RefAway/Gap.
static void fillReferenceSums(const vector<PseudoJet> &tracks, const vector<char> &inJet,
                              const vector<char> &inPerpCone, double jetEta, double jetPhi, double R,
                              int k, VariantTree &o, vector<int> *displayFlags = nullptr)
{
  for (size_t i = 0; i < tracks.size(); i++) {
    double pt = tracks[i].pt();
    if (pt < kRefPtMin || pt >= kRefPtMax) continue;
    if (inJet[i] || insideCone(tracks[i], jetEta, jetPhi, R)) continue;
    if (!inPerpCone[i]) {
      bool subeventA = tracks[i].eta() < 0;
      addToMoments(o.refN, o.refS1, o.refS2, k, pt);
      if (displayFlags) (*displayFlags)[i] |= kDisplayFlagRef;
      if (subeventA) addToMoments(o.refNa, o.refS1a, o.refS2a, k, pt);
      else           addToMoments(o.refNb, o.refS1b, o.refS2b, k, pt);
      double dphiToJet = deltaPhi(tracks[i].phi(), jetPhi);
      if (dphiToJet <= M_PI - kAwayHalfWidth) {
        addToMoments(o.refAwN, o.refAwS1, o.refAwS2, k, pt);
        if (displayFlags) (*displayFlags)[i] |= kDisplayFlagRefAway;
        if (subeventA) addToMoments(o.refAwNa, o.refAwS1a, o.refAwS2a, k, pt);
        else           addToMoments(o.refAwNb, o.refAwS1b, o.refAwS2b, k, pt);
      }
      if (dphiToJet < M_PI / 3)          addToMoments(o.nearN, o.nearS1, o.nearS2, k, pt);
      else if (dphiToJet < 2 * M_PI / 3) addToMoments(o.transN, o.transS1, o.transS2, k, pt);
      else                               addToMoments(o.awayN, o.awayS1, o.awayS2, k, pt);
    }
    if (fabs(tracks[i].eta() - jetEta) > R + kEtaGap) {
      addToMoments(o.gapN, o.gapS1, o.gapS2, k, pt);
      if (displayFlags) (*displayFlags)[i] |= kDisplayFlagGap;
    }
  }
}

// ===========================================================================
//  Correlators (-corr), see header. R = 0.4 / 0.8 (radius index 2 / 3), jet
//  classes 20-40 / 40-60 GeV/c in jetPtCorr.
// ===========================================================================
static constexpr int kNumCorrRadii = 2, kNumCorrJetClasses = 2;
static const int kCorrRadiusIndex[kNumCorrRadii] = {2, 3};
static const double kCorrJetPtLo[kNumCorrJetClasses] = {20., 40.}, kCorrJetPtHi[kNumCorrJetClasses] = {40., 60.};

struct CorrelatorHistograms {
  using Grid = TH1D *[kNumAcceptances][kNumCorrRadii][kNumCorrJetClasses];
  bool on = false;
  Grid eec, flowXS, flowXN, flowS, flowN;
  Grid eq4nS, eq4nNr, eq4n, eq4fr, eq4f, scal;
  void book(TDirectory *d, const string &variant)
  {
    on = true; d->cd();
    double logEdges[21]; for (int i = 0; i <= 20; i++) logEdges[i] = pow(10., -2 + 2. * i / 20);   // dR 0.01-1
    for (int ia = 0; ia < kNumAcceptances; ia++) for (int r = 0; r < kNumCorrRadii; r++) for (int c = 0; c < kNumCorrJetClasses; c++) {
      string p = Form("%s_a%d_R%d_j%.0f_%.0f_", variant.c_str(), ia, kCorrRadiusIndex[r], kCorrJetPtLo[c], kCorrJetPtHi[c]);
      eec[ia][r][c]    = new TH1D((p + "eec").c_str(), ";#DeltaR_{ij};#Sigma p_{T,i}p_{T,j}/p_{T,jet}^{2}", 20, logEdges);
      flowXS[ia][r][c] = new TH1D((p + "flowXS").c_str(), ";#Delta#theta to jet axis;", 25, 0, 2.5);
      flowXN[ia][r][c] = new TH1D((p + "flowXN").c_str(), ";#Delta#theta to jet axis;", 25, 0, 2.5);
      flowS[ia][r][c]  = new TH1D((p + "flowS").c_str(),  ";#Delta#theta to jet axis;", 25, 0, 2.5);
      flowN[ia][r][c]  = new TH1D((p + "flowN").c_str(),  ";#Delta#theta to jet axis;", 25, 0, 2.5);
      eq4nS[ia][r][c]  = new TH1D((p + "eq4_nS").c_str(),  ";probe bin;", kNumProbeBins, kProbePtEdges);
      eq4nNr[ia][r][c] = new TH1D((p + "eq4_nNr").c_str(), ";probe bin;", kNumProbeBins, kProbePtEdges);
      eq4n[ia][r][c]   = new TH1D((p + "eq4_n").c_str(),   ";probe bin;", kNumProbeBins, kProbePtEdges);
      eq4fr[ia][r][c]  = new TH1D((p + "eq4_fr").c_str(),  ";probe bin;", kNumProbeBins, kProbePtEdges);
      eq4f[ia][r][c]   = new TH1D((p + "eq4_f").c_str(),   ";probe bin;", kNumProbeBins, kProbePtEdges);
      scal[ia][r][c]   = new TH1D((p + "scal").c_str(), "w,NS,NNr,N,S,Nr,r,X,njet", 9, 0, 9);
    }
  }
  void write(TDirectory *d)
  {
    if (!on) return;
    d->cd();
    for (int ia = 0; ia < kNumAcceptances; ia++) for (int r = 0; r < kNumCorrRadii; r++) for (int c = 0; c < kNumCorrJetClasses; c++)
      for (TH1D *h : {eec[ia][r][c], flowXS[ia][r][c], flowXN[ia][r][c], flowS[ia][r][c], flowN[ia][r][c], eq4nS[ia][r][c],
                      eq4nNr[ia][r][c], eq4n[ia][r][c], eq4fr[ia][r][c], eq4f[ia][r][c], scal[ia][r][c]}) h->Write();
  }
};

// Fill the correlators for one leading jet (only R = 0.4, 0.8 and jetPtCorr in a
// jet class). Reference tracks: 0.2 <= pT < 3, not constituents, not in the
// perpendicular cones, not within R of the jet axis (same as refN).
static void fillCorrelators(CorrelatorHistograms &co, int ia, int iR, double w, const PseudoJet &jet, double jetPtCorr,
                            const vector<PseudoJet> &tracks, const vector<char> &inJet, const vector<char> &inPerpCone,
                            const float *inJetCounts)
{
  int r = (iR == kCorrRadiusIndex[0]) ? 0 : (iR == kCorrRadiusIndex[1] ? 1 : -1);
  if (r < 0) return;
  int c = -1; for (int i = 0; i < kNumCorrJetClasses; i++) if (jetPtCorr >= kCorrJetPtLo[i] && jetPtCorr < kCorrJetPtHi[i]) c = i;
  if (c < 0) return;
  double R = kJetRadii[iR], eta = jet.eta(), phi = jet.phi_std(), pj = jet.pt();
  vector<const PseudoJet *> constituents;
  for (size_t i = 0; i < tracks.size(); i++) if (inJet[i]) constituents.push_back(&tracks[i]);
  // X = fraction of the jet pT carried by charged-track constituents
  double X = 0; for (auto p : constituents) X += p->pt(); X /= pj;
  // energy-energy correlator inside the jet (dR floored at 0.0101: first log bin)
  for (size_t i = 0; i < constituents.size(); i++) for (size_t j = i + 1; j < constituents.size(); j++) {
    double dR = hypot(constituents[i]->eta() - constituents[j]->eta(), deltaPhi(constituents[i]->phi(), constituents[j]->phi()));
    co.eec[ia][r][c]->Fill(max(dR, 0.0101), w * constituents[i]->pt() * constituents[j]->pt() / (pj * pj));
  }
  // jet-UE energy flow vs angle of the reference tracks to the jet axis
  double S = 0, Nr = 0;
  for (size_t i = 0; i < tracks.size(); i++) {
    double pt = tracks[i].pt();
    if (pt < kRefPtMin || pt >= kRefPtMax || inJet[i] || inPerpCone[i] || insideCone(tracks[i], eta, phi, R)) continue;
    double theta = hypot(tracks[i].eta() - eta, deltaPhi(tracks[i].phi(), phi));
    co.flowXS[ia][r][c]->Fill(theta, w * X * pt); co.flowXN[ia][r][c]->Fill(theta, w * X);
    co.flowS[ia][r][c]->Fill(theta, w * pt);      co.flowN[ia][r][c]->Fill(theta, w);
    S += pt; Nr += 1;
  }
  // Eq. (4) count moments per probe bin (raw in-cone counts n_b, N = sum_b n_b)
  double N = 0; for (int b = 0; b < kNumProbeBins; b++) N += inJetCounts[b];
  double rr = Nr > 0 ? S / Nr : 0;
  for (int b = 0; b < kNumProbeBins; b++) {
    co.eq4nS[ia][r][c]->AddBinContent(b + 1, w * inJetCounts[b] * S);
    co.eq4nNr[ia][r][c]->AddBinContent(b + 1, w * inJetCounts[b] * Nr);
    co.eq4n[ia][r][c]->AddBinContent(b + 1, w * inJetCounts[b]);
    if (N > 0 && Nr > 0) { co.eq4fr[ia][r][c]->AddBinContent(b + 1, w * inJetCounts[b] / N * rr); co.eq4f[ia][r][c]->AddBinContent(b + 1, w * inJetCounts[b] / N); }
  }
  double scalars[9] = {w, w * N * S, w * N * Nr, w * N, w * S, w * Nr, Nr > 0 ? w * rr : 0, w * X, 1};
  for (int i = 0; i < 9; i++) co.scal[ia][r][c]->AddBinContent(i + 1, scalars[i]);
}

// ===========================================================================
//  Per-variant analysis
// ===========================================================================
// Event-level charged multiplicities (status 1, charge != 0).
static void countEventMultiplicities(const vector<Particle> &event, VariantTree &o)
{
  o.nchMid = o.nchFwdA = o.nchFwdC = 0;
  for (auto &p : event) {
    if (p.status != 1 || threeTimesCharge(p.pid) == 0) continue;
    double pt = transverseMomentum(p), eta = pseudorapidity(p.px, p.py, p.pz);
    if (fabs(eta) < 0.8 && pt >= 0.2 && pt < 5.0) o.nchMid++;   // mid-rapidity, 0.2-5 GeV/c
    // forward centrality estimators: charged particles with 0.2 < pT < 10 GeV/c
    // (author 2026-10-29; the 2026-10-28 production trees have no pT cut here)
    const bool fwdPt = (pt > kFwdPtMin && pt < kFwdPtMax);
    if (fwdPt && eta > 3.5 && eta < 4.9) o.nchFwdA++;           // FT0A-like
    if (fwdPt && eta > -3.3 && eta < -2.1) o.nchFwdC++;         // FT0C-like
  }
}

// Raw jets for every R in one acceptance: leading jet with active area,
// perpendicular-cone background (rho, nPerp), jetPtCorr = pT - rho A, in-cone
// counts nIn, reference sums and correlators. Returns the event rho for the
// constituent subtraction (rho of the R = 0.4 perpendicular cones; 0 if there
// is no R = 0.4 jet). display: event-display record (nullptr in production).
static double analyseRawJets(const vector<PseudoJet> &tracks, int ia, VariantTree &o, CorrelatorHistograms *co,
                             VariantDisplay *display = nullptr)
{
  double etaMax = kTrackEtaMax[ia];
  double rhoForCS = 0;
  for (int iR = 0; iR < kNumJetRadii; iR++) {
    double R = kJetRadii[iR];
    int k = ia * kNumJetRadii + iR;
    PseudoJet jet; double area; unique_ptr<fastjet::ClusterSequence> clusterSeq;
    const bool haveJet = findLeadingJet(tracks, R, etaMax, true, jet, area, clusterSeq);
    DisplaySlotRecord *dsp = (display && display->wants(k)) ? &display->slot[k] : nullptr;
    if (dsp) dsp->recordJets(*clusterSeq, haveJet ? &jet : nullptr, etaMax - R);
    if (!haveJet) continue;
    double eta = jet.eta(), phi = jet.phi_std();

    vector<char> inPerpCone(tracks.size(), 0);
    double rho = fillPerpendicularCones(tracks, eta, phi, R, &o.nPerp[k * kNumProbeBins], inPerpCone);
    if (fabs(R - kRhoConeR) < 1e-6) rhoForCS = rho;

    o.jetPt[k] = jet.pt(); o.jetEta[k] = eta; o.jetPhi[k] = phi; o.jetArea[k] = area;
    o.rho[k] = rho; o.jetPtCorr[k] = jet.pt() - rho * area;          // area-based UE subtraction
    vector<char> inJet(tracks.size(), 0);
    for (auto &c : jet.constituents()) {
      if (c.is_pure_ghost()) continue;
      inJet[c.user_index()] = 1;
      o.leadTrk[k] = max<float>(o.leadTrk[k], c.pt());
      int b = probePtBin(c.pt()); if (b >= 0) o.nIn[k * kNumProbeBins + b] += 1;
    }
    if (dsp)
      for (size_t i = 0; i < tracks.size(); i++) {
        if (inJet[i])      dsp->trkFlags[i] |= kDisplayFlagLeadJet;
        if (inPerpCone[i]) dsp->trkFlags[i] |= kDisplayFlagPerpCone;
        if (insideCone(tracks[i], eta, phi, R)) dsp->trkFlags[i] |= kDisplayFlagJetCone;
      }

    fillReferenceSums(tracks, inJet, inPerpCone, eta, phi, R, k, o, dsp ? &dsp->trkFlags : nullptr);
    if (co && co->on) fillCorrelators(*co, ia, iR, o.weight, jet, o.jetPtCorr[k], tracks, inJet, inPerpCone, &o.nIn[k * kNumProbeBins]);
  }
  return rhoForCS;
}

// Event-wise constituent subtraction (fjcontrib): distance deltaR, dR_max 0.25,
// alpha 0, ghost area 0.01 up to |eta| = acceptance, scalar background rhoForCS
// (rho_m = 0). Then the leading jet (no area) per R and its constituents nCS.
// The subtracted particles keep the user_index of their input track.
static void analyseSubtractedJets(const vector<PseudoJet> &tracks, int ia, double rhoForCS, VariantTree &o,
                                  VariantDisplay *display = nullptr)
{
  double etaMax = kTrackEtaMax[ia];
  fastjet::contrib::ConstituentSubtractor subtractor;
  subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
  subtractor.set_max_distance(kCSMaxDeltaR);
  subtractor.set_alpha(kCSAlpha);
  subtractor.set_ghost_area(kCSGhostArea);
  subtractor.set_max_eta(etaMax);
  subtractor.set_scalar_background_density(rhoForCS, 0.);
  subtractor.initialize();
  vector<PseudoJet> subtractedTracks = subtractor.subtract_event(tracks, etaMax);
  if (display)
    for (int iR = 0; iR < kNumJetRadii; iR++) {
      int k = ia * kNumJetRadii + iR;
      if (!display->wants(k)) continue;
      DisplaySlotRecord &dsp = display->slot[k];
      for (auto &p : subtractedTracks) {
        int i = p.user_index();
        if (i < 0 || i >= (int)tracks.size() || deltaPhi(p.phi(), tracks[i].phi()) > 1e-6)
          { cerr << "event display: CS output particle without its input track\n"; continue; }
        dsp.trkPtCS[i] = p.pt();
      }
      for (auto &g : subtractor.get_ghosts()) { dsp.csGhostEta.push_back(g.eta()); dsp.csGhostPhi.push_back(g.phi_std()); }
    }
  for (int iR = 0; iR < kNumJetRadii && !subtractedTracks.empty(); iR++) {
    double R = kJetRadii[iR];
    int k = ia * kNumJetRadii + iR;
    PseudoJet jet; double area; unique_ptr<fastjet::ClusterSequence> clusterSeq;
    if (findLeadingJet(subtractedTracks, R, etaMax, false, jet, area, clusterSeq)) {
      o.csJetPt[k] = jet.pt(); o.csJetEta[k] = jet.eta(); o.csJetPhi[k] = jet.phi_std();
      for (auto &c : jet.constituents()) {
        o.csLeadTrk[k] = max<float>(o.csLeadTrk[k], c.pt());
        int b = probePtBin(c.pt()); if (b >= 0) o.nCS[k * kNumProbeBins + b] += 1;
        if (display && display->wants(k)) display->slot[k].trkFlags[c.user_index()] |= kDisplayFlagCSJet;
      }
    }
  }
}

// All jet quantities of one variant (weight/centr/ueMult are set by the caller).
// display: event-display record of this variant (nullptr in production).
static void analyseVariant(const vector<Particle> &event, VariantTree &o, CorrelatorHistograms *co = nullptr,
                           VariantDisplay *display = nullptr)
{
  o.reset();
  countEventMultiplicities(event, o);
  for (int ia = 0; ia < kNumAcceptances; ia++) {
    vector<int> sourceIndex;
    vector<PseudoJet> tracks = selectChargedTracks(event, kTrackEtaMax[ia], display ? &sourceIndex : nullptr);
    if (display)
      for (int iR = 0; iR < kNumJetRadii; iR++) {
        int k = ia * kNumJetRadii + iR;
        if (!display->wants(k)) continue;
        display->slot[k].clear();
        display->slot[k].recordTracks(tracks, sourceIndex, display->origin);
      }
    if (tracks.empty()) continue;
    double rhoForCS = analyseRawJets(tracks, ia, o, co, display);
    o.rhoEv[ia] = rhoForCS;
    analyseSubtractedJets(tracks, ia, rhoForCS, o, display);
  }
}

// Truth matching: copy the leading-jet pT of the no-UE variant (nom, or cs10 for
// the recoil sample; for that variant itself truth == o) and the (eta, phi)
// distance between the two leading jets, then fill the tree.
static void fillWithTruthMatch(VariantTree &o, const VariantTree &truth)
{
  for (int k = 0; k < kNumAccR; k++) {
    o.jetPtTruth[k] = truth.jetPt[k];
    if (o.jetPt[k] > 0 && truth.jetPt[k] > 0)
      o.dRTruth[k] = hypot(o.jetEta[k] - truth.jetEta[k], deltaPhi(o.jetPhi[k], truth.jetPhi[k]));
  }
  o.tree->Fill();
}

// ===========================================================================
//  Truth-level spectra (-spectra), see header. Classes in centrality percent.
// ===========================================================================
static constexpr int kNumSpecClasses = 4, kNumSpecRadii = 2;
static const double kSpecCentLo[kNumSpecClasses] = {0, 0, 40, 0}, kSpecCentHi[kNumSpecClasses] = {10, 40, 80, 100};
static const double kSpecJetRadii[kNumSpecRadii] = {0.2, 0.4};

struct TruthSpectra {
  bool on = false;
  TH1D *track[kNumSpecClasses], *jet[kNumSpecRadii][kNumSpecClasses], *norm[kNumSpecClasses];
  void book(TDirectory *d)
  {
    on = true; d->cd();
    double logEdges[61]; for (int i = 0; i <= 60; i++) logEdges[i] = 0.2 * pow(500., i / 60.);   // 0.2-100 GeV/c
    for (int c = 0; c < kNumSpecClasses; c++) {
      string s = Form("c%.0f_%.0f", kSpecCentLo[c], kSpecCentHi[c]);
      track[c] = new TH1D(("trk_" + s).c_str(), ";p_{T} (GeV/c);#Sigma w dN/dp_{T} (|#eta|<0.9)", 60, logEdges);
      for (int r = 0; r < kNumSpecRadii; r++)
        jet[r][c] = new TH1D(Form("jet_R%02.0f_%s", 10 * kSpecJetRadii[r], s.c_str()),
                             ";p_{T,ch jet} (GeV/c);#Sigma w dN/dp_{T}", 50, 0, 250);
      norm[c] = new TH1D(("norm_" + s).c_str(), "sum w, n events, sum w^2", 3, 0, 3);
      track[c]->Sumw2(); for (int r = 0; r < kNumSpecRadii; r++) jet[r][c]->Sumw2();
    }
  }
  // event: truth particle list (no added UE); centrality < 0: only class 0-100
  void fill(const vector<Particle> &event, double w, double centrality)
  {
    vector<int> classes;
    for (int c = 0; c < kNumSpecClasses; c++)
      if (c == kNumSpecClasses - 1 || (centrality >= 0 && centrality >= kSpecCentLo[c] && centrality < kSpecCentHi[c])) classes.push_back(c);
    for (int c : classes) { norm[c]->AddBinContent(1, w); norm[c]->AddBinContent(2, 1); norm[c]->AddBinContent(3, w * w); }
    vector<PseudoJet> tracks = selectChargedTracks(event, 0.9);
    for (auto &t : tracks) if (t.pt() >= 0.2) for (int c : classes) track[c]->Fill(t.pt(), w);
    for (int r = 0; r < kNumSpecRadii; r++) {
      fastjet::ClusterSequence cs(tracks, fastjet::JetDefinition(fastjet::antikt_algorithm, kSpecJetRadii[r]));
      for (auto &j : (fastjet::SelectorAbsEtaMax(0.9 - kSpecJetRadii[r]) * fastjet::SelectorPtMin(5.))(cs.inclusive_jets()))
        for (int c : classes) jet[r][c]->Fill(j.pt(), w);
    }
  }
  void write(TDirectory *d)
  {
    if (!on) return;
    d->cd();
    for (int c = 0; c < kNumSpecClasses; c++) { track[c]->Write(); norm[c]->Write(); for (int r = 0; r < kNumSpecRadii; r++) jet[r][c]->Write(); }
  }
};

// ===========================================================================
//  Event-display file (-display), see header
// ===========================================================================
class EventDisplayWriter {
public:
  bool open(const string &fileName, bool recoil)
  {
    TDirectory *previous = gDirectory;
    file.reset(new TFile(fileName.c_str(), "RECREATE"));
    if (file->IsZombie()) return false;
    tree = new TTree("display", "jetv0ana event display: one entry per (event, variant, acceptance, R)");
    tree->Branch("event", &event, "event/I");
    tree->Branch("variant", &variant);
    tree->Branch("iacc", &iacc, "iacc/I"); tree->Branch("iR", &iR, "iR/I");
    for (auto &b : vector<pair<const char *, float *>>{
           {"R", &R}, {"trackEtaMax", &trackEtaMax}, {"jetEtaMax", &jetEtaMax}, {"weight", &weight}, {"centr", &centr},
           {"ueMult", &ueMult}, {"leadPt", &leadPt}, {"leadPtCorr", &leadPtCorr}, {"leadEta", &leadEta},
           {"leadPhi", &leadPhi}, {"leadArea", &leadArea}, {"leadTrk", &leadTrk}, {"rho", &rho}, {"rhoCS", &rhoCS},
           {"csJetPt", &csJetPt}, {"csJetEta", &csJetEta}, {"csJetPhi", &csJetPhi}, {"jetPtTruth", &jetPtTruth},
           {"dRTruth", &dRTruth}, {"refN", &refN}, {"refS1", &refS1}, {"jetGhostArea", &jetGhostArea},
           {"csGhostArea", &csGhostArea}, {"csMaxDR", &csMaxDR}})
      tree->Branch(b.first, b.second, (string(b.first) + "/F").c_str());
    tree->Branch("leadJet", &rec.leadJet, "leadJet/I");
    tree->Branch("perpPhi", perpPhi, "perpPhi[2]/F");
    tree->Branch("probeEdges", probeEdges, Form("probeEdges[%d]/F", kNumProbeBins + 1));
    for (auto &b : vector<pair<const char *, float *>>{{"nIn", nIn}, {"nPerp", nPerp}, {"nSub", nSub}, {"nCS", nCS},
                                                       {"fRaw", fRaw}, {"fSub", fSub}, {"fCS", fCS}})
      tree->Branch(b.first, b.second, Form("%s[%d]/F", b.first, kNumProbeBins));
    for (auto &b : vector<pair<const char *, vector<float> *>>{
           {"trkPt", &rec.trkPt}, {"trkEta", &rec.trkEta}, {"trkPhi", &rec.trkPhi}, {"trkPtCS", &rec.trkPtCS},
           {"jetPt", &rec.jetPt}, {"jetEta", &rec.jetEta}, {"jetPhi", &rec.jetPhi}, {"jetArea", &rec.jetArea},
           {"ghostEta", &rec.ghostEta}, {"ghostPhi", &rec.ghostPhi}, {"csGhostEta", &rec.csGhostEta}, {"csGhostPhi", &rec.csGhostPhi}})
      tree->Branch(b.first, b.second);
    for (auto &b : vector<pair<const char *, vector<int> *>>{
           {"trkOrigin", &rec.trkOrigin}, {"trkJet", &rec.trkJet}, {"trkFlags", &rec.trkFlags},
           {"jetNConst", &rec.jetNConst}, {"jetFiducial", &rec.jetFiducial}, {"ghostJet", &rec.ghostJet}})
      tree->Branch(b.first, b.second);
    if (recoil) {
      recoilTree = new TTree("recoil", "jetv0ana event display: recoil thermal-momentum subtraction");
      recoilTree->Branch("event", &event, "event/I");
      recoilTree->Branch("dRmax", &dRmax, "dRmax/F");
      recoilTree->Branch("nPairs", &nPairs, "nPairs/L");
      for (auto &b : vector<pair<const char *, vector<float> *>>{
             {"finPt", &recoilRec.finalPt}, {"finY", &recoilRec.finalY}, {"finPhi", &recoilRec.finalPhi},
             {"finPtAfter", &recoilRec.finalPtAfter}, {"thPt", &recoilRec.thermalPt}, {"thY", &recoilRec.thermalY},
             {"thPhi", &recoilRec.thermalPhi}, {"thPtAfter", &recoilRec.thermalPtAfter}})
        recoilTree->Branch(b.first, b.second);
      recoilTree->Branch("finCharged", &recoilRec.finalCharged);
    }
    for (int b = 0; b <= kNumProbeBins; b++) probeEdges[b] = kProbePtEdges[b];
    if (previous) previous->cd();
    return true;
  }
  // one slot k of one variant; o is the variant's tree record after the truth match
  void fill(long eventNumber, const VariantTree &o, int k, const DisplaySlotRecord &slot)
  {
    event = (int)eventNumber; variant = o.name; iacc = k / kNumJetRadii; iR = k % kNumJetRadii;
    R = kJetRadii[iR]; trackEtaMax = kTrackEtaMax[iacc]; jetEtaMax = trackEtaMax - R;
    weight = o.weight; centr = o.centr; ueMult = o.ueMult;
    leadPt = o.jetPt[k]; leadPtCorr = o.jetPtCorr[k]; leadEta = o.jetEta[k]; leadPhi = o.jetPhi[k];
    leadArea = o.jetArea[k]; leadTrk = o.leadTrk[k]; rho = o.rho[k]; rhoCS = o.rhoEv[iacc];
    csJetPt = o.csJetPt[k]; csJetEta = o.csJetEta[k]; csJetPhi = o.csJetPhi[k];
    jetPtTruth = o.jetPtTruth[k]; dRTruth = o.dRTruth[k]; refN = o.refN[k]; refS1 = o.refS1[k];
    jetGhostArea = kJetGhostArea; csGhostArea = kCSGhostArea; csMaxDR = kCSMaxDeltaR;
    perpPhi[0] = TVector2::Phi_mpi_pi(leadPhi - M_PI / 2); perpPhi[1] = TVector2::Phi_mpi_pi(leadPhi + M_PI / 2);
    // probe spectra and fractions as in the offline analysis (Phase1_mean.C)
    double sumRaw = 0, sumSub = 0, sumCS = 0;
    for (int b = 0; b < kNumProbeBins; b++) {
      nIn[b] = o.nIn[k * kNumProbeBins + b]; nPerp[b] = o.nPerp[k * kNumProbeBins + b]; nCS[b] = o.nCS[k * kNumProbeBins + b];
      nSub[b] = nIn[b] - nPerp[b]; sumRaw += nIn[b]; sumSub += nSub[b]; sumCS += nCS[b];
    }
    for (int b = 0; b < kNumProbeBins; b++) {
      fRaw[b] = sumRaw > 0 ? nIn[b] / sumRaw : 0; fSub[b] = sumSub > 0 ? nSub[b] / sumSub : 0; fCS[b] = sumCS > 0 ? nCS[b] / sumCS : 0;
    }
    rec = slot;
    tree->Fill();
  }
  void fillRecoil(long eventNumber, const RecoilSubtractionRecord &r)
  {
    if (!recoilTree) return;
    event = (int)eventNumber; dRmax = r.dRmax; nPairs = r.nPairs; recoilRec = r;
    recoilTree->Fill();
  }
  void close()
  {
    if (!file) return;
    file->cd(); tree->Write(); if (recoilTree) recoilTree->Write();
    printf("event display: %lld entries (display) and %lld (recoil) -> %s\n", tree->GetEntries(),
           recoilTree ? recoilTree->GetEntries() : 0LL, file->GetName());
    file->Close(); file.reset();
  }
private:
  unique_ptr<TFile> file;
  TTree *tree = nullptr, *recoilTree = nullptr;
  int event = 0, iacc = 0, iR = 0; string variant;
  float R, trackEtaMax, jetEtaMax, weight, centr, ueMult, leadPt, leadPtCorr, leadEta, leadPhi, leadArea, leadTrk;
  float rho, rhoCS, csJetPt, csJetEta, csJetPhi, jetPtTruth, dRTruth, refN, refS1, jetGhostArea, csGhostArea, csMaxDR;
  float perpPhi[2], probeEdges[kNumProbeBins + 1];
  float nIn[kNumProbeBins], nPerp[kNumProbeBins], nSub[kNumProbeBins], nCS[kNumProbeBins];
  float fRaw[kNumProbeBins], fSub[kNumProbeBins], fCS[kNumProbeBins];
  DisplaySlotRecord rec;
  float dRmax = 0; Long64_t nPairs = 0;
  RecoilSubtractionRecord recoilRec;
};

// ===========================================================================
//  Command line
// ===========================================================================
struct UEEmbedSpec { string tag, listFile; double lo, hi; };
struct Options {
  string hepmcFile, outFile = "jetv0.root";
  long nPythia = 0, maxEvents = -1;
  bool recoil = false, noMPI = false, thermCent = false, doCorr = false;
  bool thermFlow = false, doSpectra = false;
  vector<double> thermPtScan;                               // window <pT> values
  double thermPtMatch = -1;                                 // window <pT>, < 0: off
  vector<UEEmbedSpec> ueEmbeds;
  unsigned seed = 1;
  double ecm = 5360., ptHatMin = 5.;
  string displayFile;                                       // event display, "": off
  long displayEvents = 5;
  vector<string> displayVariants;                           // empty: no-UE variant + embMB
};

static const char *kUsage =
  "usage: jetv0ana (-hepmc file [-recoil] | -pythia N [-ecm 5360] [-pthatmin 5] [-nompi]) -out f.root -seed s [-nev N]\n"
  "       [-thermcent [-thermflow]] [-thermpt X] [-thermptscan X]... [-ueembed TAG list lo hi]... [-corr] [-spectra]\n"
  "       [-display file.root [-displayEvents 5] [-displayVariants A,B,...]]\n";

// returns false (after printing the reason) on invalid input
static bool parseOptions(int argc, char **argv, Options &opt)
{
  for (int i = 1; i < argc; i++) {
    string a = argv[i];
    auto nextArg = [&]() { return string(i + 1 < argc ? argv[++i] : ""); };
    if (a == "-hepmc") opt.hepmcFile = nextArg();
    else if (a == "-pythia") opt.nPythia = atol(nextArg().c_str());
    else if (a == "-out") opt.outFile = nextArg();
    else if (a == "-recoil") opt.recoil = true;
    else if (a == "-nompi") opt.noMPI = true;
    else if (a == "-seed") opt.seed = atoi(nextArg().c_str());
    else if (a == "-nev") opt.maxEvents = atol(nextArg().c_str());
    else if (a == "-ecm") opt.ecm = atof(nextArg().c_str());
    else if (a == "-pthatmin") opt.ptHatMin = atof(nextArg().c_str());
    else if (a == "-thermcent") opt.thermCent = true;
    else if (a == "-thermpt") opt.thermPtMatch = atof(nextArg().c_str());
    else if (a == "-corr") opt.doCorr = true;
    else if (a == "-thermflow") opt.thermFlow = true;
    else if (a == "-thermptscan") opt.thermPtScan.push_back(atof(nextArg().c_str()));
    else if (a == "-spectra") opt.doSpectra = true;
    else if (a == "-display") opt.displayFile = nextArg();
    else if (a == "-displayEvents") opt.displayEvents = atol(nextArg().c_str());
    else if (a == "-displayVariants") {
      stringstream list(nextArg()); string v;
      while (getline(list, v, ',')) if (!v.empty()) opt.displayVariants.push_back(v);
    }
    else if (a == "-ueembed") {
      UEEmbedSpec u; u.tag = nextArg(); u.listFile = nextArg(); u.lo = atof(nextArg().c_str()); u.hi = atof(nextArg().c_str());
      opt.ueEmbeds.push_back(u);
    }
    else { cerr << "unknown option " << a << endl; return false; }
  }
  if (opt.hepmcFile.empty() == (opt.nPythia == 0)) { cerr << kUsage; return false; }
  if (opt.thermFlow && !opt.thermCent) { cerr << "-thermflow needs -thermcent\n"; return false; }
  return true;
}

// ===========================================================================
//  main
// ===========================================================================
int main(int argc, char **argv)
{
  Options opt;
  if (!parseOptions(argc, argv, opt)) return 1;
  const bool recoil = opt.recoil;
  fastjet::ClusterSequence::set_fastjet_banner_stream(nullptr);

  // Two generators: rng for the original (v1/v2) thermal variants and the drawn
  // centrality, rngV3 for embFlow and the pT scan. Keeping them separate leaves
  // the random streams of the older variants unchanged when v3 options are added.
  TRandom3 rng(opt.seed);
  TRandom3 rngV3(opt.seed + 7777777u);

  // ---- variant list (tree order = order of the particle lists built per event)
  const string embPrefix = recoil ? "cs10_" : "";            // added UE is built on cs10 for recoils
  vector<VariantTree> variants;
  auto addVariant = [&](const string &name) { VariantTree v; v.name = name; variants.push_back(v); };
  for (auto &b : recoil ? vector<string>{"nosub", "cs10", "cs05"} : vector<string>{"nom"}) addVariant(b);
  for (auto &th : kThermalConfigs) addVariant(embPrefix + th.tag);

  OpticalGlauberOxygen glauber;
  if (opt.thermCent) {
    glauber.build();
    printf("Glauber 16O: dN/deta = %.4f * Npart^%.4f ; <Npart> 0-5%% %.2f, 0-40%% %.2f, 40-80%% %.2f, 0-100%% %.2f ;"
           " <dN/deta> 0-5%% %.1f, 0-40%% %.1f, 40-80%% %.1f, 0-100%% %.1f\n", glauber.A, glauber.beta,
           glauber.classMean(0, 5, false), glauber.classMean(0, 40, false), glauber.classMean(40, 80, false), glauber.classMean(0, 100, false),
           glauber.classMean(0, 5, true), glauber.classMean(0, 40, true), glauber.classMean(40, 80, true), glauber.classMean(0, 100, true));
    addVariant(embPrefix + "embNp");
  }
  double matchedSlope = 0;
  if (opt.thermPtMatch > 0) {
    matchedSlope = slopeForWindowMeanPt(opt.thermPtMatch, kRefPtMin, kRefPtMax);
    printf("thermal slope T = %.4f GeV: window mean %.4f (target %.4f); default embMB/embC window mean %.4f\n",
           matchedSlope, gammaWindowMeanPt(matchedSlope, kRefPtMin, kRefPtMax), opt.thermPtMatch,
           gammaWindowMeanPt(kThermalDefaultSlope, kRefPtMin, kRefPtMax));
    for (auto &th : kThermalConfigs) addVariant(embPrefix + th.tag + "match");
  }
  vector<unique_ptr<UEEventSource>> ueSources;
  for (auto &u : opt.ueEmbeds) {
    unique_ptr<UEEventSource> s(new UEEventSource); s->tag = u.tag; s->lo = u.lo; s->hi = u.hi;
    if (!s->open(u.listFile, opt.seed)) { cerr << "cannot open UE list " << u.listFile << endl; return 1; }
    ueSources.push_back(move(s));
    addVariant(embPrefix + "emb" + u.tag);
  }
  // v3 variants (embFlow, pT scan) come last; see the ghost-state note in the event loop
  const size_t firstV3Variant = variants.size();
  FlowingUEModel flowModel;
  if (opt.thermFlow) {
    flowModel.build(glauber);
    printf("flow UE: window <pT> = %.3f + %.4f*(ln dN/deta - %.4f); class means 0-5%% %.3f, 0-10%% %.3f, 0-40%% %.3f, "
           "40-80%% %.3f, 70-80%% %.3f, 0-100%% %.3f GeV/c; T fluct %.0f%%, mult fluct 1/sqrt(Npart)\n",
           flowModel.M_MB, flowModel.K, flowModel.L0, flowModel.classWindowMeanPt(glauber, 0, 5),
           flowModel.classWindowMeanPt(glauber, 0, 10), flowModel.classWindowMeanPt(glauber, 0, 40),
           flowModel.classWindowMeanPt(glauber, 40, 80), flowModel.classWindowMeanPt(glauber, 70, 80),
           flowModel.classWindowMeanPt(glauber, 0, 100), 100 * flowModel.TFLUC);
    addVariant(embPrefix + "embFlow");
  }
  vector<double> scanSlopes;
  for (double x : opt.thermPtScan) {
    scanSlopes.push_back(slopeForWindowMeanPt(x, kRefPtMin, kRefPtMax));
    printf("pT scan: window mean %.3f -> T = %.4f GeV\n", x, scanSlopes.back());
    addVariant(embPrefix + Form("embMBpt%03ld", lround(100 * x)));
  }

  // ---- output file (trees booked only after the variant vector is final:
  //      the branch addresses point into its elements)
  TFile fout(opt.outFile.c_str(), "RECREATE");
  fout.SetCompressionSettings(505);                          // ZSTD level 5
  for (auto &v : variants) v.book(&fout);
  vector<CorrelatorHistograms> correlators(variants.size());
  TDirectory *corrDir = nullptr;
  if (opt.doCorr) {
    corrDir = fout.mkdir("corr");
    for (size_t i = 0; i < variants.size(); i++) correlators[i].book(corrDir, variants[i].name);
    fout.cd();
  }
  TruthSpectra spectra;
  TDirectory *specDir = nullptr;
  if (opt.doSpectra) { specDir = fout.mkdir("spectra"); spectra.book(specDir); fout.cd(); }

  // ---- event display (separate file; reads the objects of the analysis only)
  unique_ptr<EventDisplayWriter> display;
  vector<char> displayVariant(variants.size(), 0);
  if (!opt.displayFile.empty()) {
    vector<string> names = opt.displayVariants;
    if (names.empty()) names = {recoil ? "cs10" : "nom", embPrefix + "embMB"};
    for (auto &n : names) {
      size_t v = 0; while (v < variants.size() && variants[v].name != n) v++;
      if (v == variants.size()) { cerr << "-displayVariants: no variant " << n << endl; return 1; }
      displayVariant[v] = 1;
    }
    display.reset(new EventDisplayWriter);
    if (!display->open(opt.displayFile, recoil)) { cerr << "cannot open " << opt.displayFile << endl; return 1; }
    fout.cd();
  }
  vector<VariantDisplay> variantDisplays(variants.size());

  // ---- input
  unique_ptr<HepMC2Reader> reader;
#ifdef WITH_PYTHIA8
  unique_ptr<Pythia8::Pythia> pythia;
#endif
  if (!opt.hepmcFile.empty()) {
    reader.reset(new HepMC2Reader(opt.hepmcFile));
    if (!reader->ok()) { cerr << "cannot open " << opt.hepmcFile << endl; return 1; }
  } else {
#ifdef WITH_PYTHIA8
    pythia.reset(new Pythia8::Pythia("", false));
    pythia->readString("Beams:eCM = " + to_string(opt.ecm));
    pythia->readString("HardQCD:all = on");
    pythia->readString("PhaseSpace:pTHatMin = " + to_string(opt.ptHatMin));
    pythia->readString("PhaseSpace:bias2Selection = on");       // pT-hat^4.5 weighting, as JEWEL WEXPO 4.5
    pythia->readString("PhaseSpace:bias2SelectionPow = 4.5");
    pythia->readString("PhaseSpace:bias2SelectionRef = " + to_string(opt.ptHatMin));
    pythia->readString("Random:setSeed = on");
    pythia->readString("Random:seed = " + to_string(opt.seed % 900000000));
    pythia->readString("Next:numberCount = 0");
    if (opt.noMPI) pythia->readString("PartonLevel:MPI = off");
    if (!pythia->init()) { cerr << "pythia init failed\n"; return 1; }
#else
    cerr << "built without PYTHIA 8\n"; return 1;
#endif
  }

  // ---- event loop
  vector<Particle> event;
  long nEvents = 0;
  while (opt.maxEvents < 0 || nEvents < opt.maxEvents) {
    double weight = 1, centrality = -1;
    if (reader) {
      if (!reader->next(event, weight, centrality)) break;
    } else {
#ifdef WITH_PYTHIA8
      if (nEvents >= opt.nPythia) break;
      if (!pythia->next()) continue;
      event.clear();
      for (int i = 0; i < pythia->event.size(); i++) {
        auto &p = pythia->event[i];
        if (!p.isFinal()) continue;
        event.push_back({p.px(), p.py(), p.pz(), p.e(), p.id(), 1});
      }
      weight = pythia->info.weight();
#endif
    }
    nEvents++;

    // Particle lists, one per variant and in the same order. ueDndeta and
    // listCentrality are written to ueMult and centr of the variant.
    vector<vector<Particle>> lists;
    vector<Particle> noUE;                                   // base of every added-UE list (nom or cs10)
    const bool displayThisEvent = display && nEvents <= opt.displayEvents;
    RecoilSubtractionRecord recoilLarge, recoilSmall;        // event display only
    if (recoil) {
      lists.push_back(finalStateParticles(event));                               // nosub
      lists.push_back(subtractRecoilThermalMomenta(event, kRecoilSubDRmaxLarge, displayThisEvent ? &recoilLarge : nullptr)); // cs10
      lists.push_back(subtractRecoilThermalMomenta(event, kRecoilSubDRmaxSmall, displayThisEvent ? &recoilSmall : nullptr)); // cs05
      noUE = lists[1];
    } else {
      lists.push_back(finalStateParticles(event));                               // nom
      noUE = lists[0];
    }
    // embMB, embC (rng)
    for (auto &th : kThermalConfigs) lists.push_back(addThermalPions(noUE, th.dndeta, th.meanPt / 2., rng));
    vector<float> ueDndeta(lists.size(), 0.f), listCentrality(lists.size(), (float)centrality);
    for (size_t v = 0; v < kThermalConfigs.size(); v++) ueDndeta[lists.size() - kThermalConfigs.size() + v] = kThermalConfigs[v].dndeta;

    // embNp (rng): event centrality, or a uniform draw for events without one
    double ueCentrality = -1;
    if (opt.thermCent) {
      double c = (centrality > 0 && centrality <= 100) ? centrality : 100. * rng.Rndm();
      ueCentrality = c;
      double dndeta = glauber.dndetaAt(c);
      lists.push_back(addThermalPions(noUE, dndeta, kThermalDefaultSlope, rng)); ueDndeta.push_back(dndeta); listCentrality.push_back(c);
    }
    // embMBmatch, embCmatch (rng)
    if (opt.thermPtMatch > 0)
      for (auto &th : kThermalConfigs) {
        lists.push_back(addThermalPions(noUE, th.dndeta, matchedSlope, rng)); ueDndeta.push_back(th.dndeta); listCentrality.push_back(centrality);
      }
    // emb<TAG>: generator UE events (no random numbers)
    for (auto &s : ueSources) {
      vector<Particle> ue, combined = noUE;
      if (!s->next(ue)) { cerr << "UE source " << s->tag << " has no event in window\n"; return 1; }
      int nChargedMid = 0;
      for (auto &p : ue) if (threeTimesCharge(p.pid) != 0 && fabs(pseudorapidity(p.px, p.py, p.pz)) < 0.5) nChargedMid++;
      combined.insert(combined.end(), ue.begin(), ue.end());
      lists.push_back(combined); ueDndeta.push_back(nChargedMid); listCentrality.push_back(centrality);
    }
    // embFlow (rngV3): one Gaussian z shared by multiplicity and slope fluctuation
    if (opt.thermFlow) {
      double z = rngV3.Gaus();
      double dndeta = max(0., glauber.dndetaAt(ueCentrality) * (1 + z / sqrt(max(glauber.npartAt(ueCentrality), 2.))));
      double T = flowModel.inverse.slopeFor(flowModel.windowMeanPt(glauber, ueCentrality)) * (1 + flowModel.TFLUC * z);
      lists.push_back(addThermalPions(noUE, dndeta, max(T, 0.05), rngV3)); ueDndeta.push_back(dndeta); listCentrality.push_back(ueCentrality);
    }
    // embMBptNNN (rngV3): MB density, scanned slope
    for (double T : scanSlopes) {
      lists.push_back(addThermalPions(noUE, kThermalConfigs[0].dndeta, T, rngV3));
      ueDndeta.push_back(kThermalConfigs[0].dndeta); listCentrality.push_back(centrality);
    }

    // Analyse every variant. The active-area ghosts use FastJet's global random
    // generator: its state is saved before the first v3 variant and restored
    // after the last, so the v1/v2 variants see the same ghost sequence as in
    // a run without v3 options (their jet areas and branches stay unchanged).
    // Event display: particle origin per list (0 generator, 1 added UE, 2 reduced
    // by the recoil thermal subtraction); lists after the base lists are noUE + UE.
    if (displayThisEvent) {
      const size_t nBase = recoil ? 3 : 1;
      vector<int> noUEOrigin = recoil ? recoilLarge.outputOrigin() : vector<int>(noUE.size(), 0);
      for (size_t v = 0; v < variants.size(); v++) {
        vector<int> &origin = variantDisplays[v].origin;
        if (v >= nBase)       { origin.assign(lists[v].size(), 1); copy(noUEOrigin.begin(), noUEOrigin.end(), origin.begin()); }
        else if (!recoil || v == 0) origin.assign(lists[v].size(), 0);
        else                  origin = (v == 1 ? recoilLarge : recoilSmall).outputOrigin();
      }
    }
    vector<int> ghostRandomState;
    for (size_t v = 0; v < variants.size(); v++) {
      if (v == firstV3Variant) fastjet::GhostedAreaSpec().get_random_status(ghostRandomState);
      variants[v].weight = weight; variants[v].centr = listCentrality[v]; variants[v].ueMult = ueDndeta[v];
      analyseVariant(lists[v], variants[v], opt.doCorr ? &correlators[v] : nullptr,
                     displayThisEvent && displayVariant[v] ? &variantDisplays[v] : nullptr);
    }
    if (!ghostRandomState.empty()) fastjet::GhostedAreaSpec().set_random_status(ghostRandomState);
    const VariantTree &truth = variants[recoil ? 1 : 0];     // no-UE variant: nom, or cs10 for recoils
    for (auto &v : variants) fillWithTruthMatch(v, truth);
    if (displayThisEvent) {
      for (size_t v = 0; v < variants.size(); v++)
        if (displayVariant[v]) for (int k : kDisplaySlots) display->fill(nEvents, variants[v], k, variantDisplays[v].slot[k]);
      if (recoil) { display->fillRecoil(nEvents, recoilLarge); display->fillRecoil(nEvents, recoilSmall); }
    }
    if (opt.doSpectra) spectra.fill(noUE, weight, centrality);
    if (nEvents % 1000 == 0) cout << "event " << nEvents << endl;
  }

  // ---- write
  if (display) display->close();
  fout.cd();
  for (auto &v : variants) v.tree->Write();
  if (corrDir) for (auto &c : correlators) c.write(corrDir);
  spectra.write(specDir);
  for (auto &s : ueSources) printf("UE source %s: scanned %ld, used %ld (window %.0f-%.0f)\n", s->tag.c_str(), s->scanned, s->used, s->lo, s->hi);
  fout.cd();
  fout.Close();
  cout << "jetv0ana: processed " << nEvents << " events -> " << opt.outFile << endl;
  return 0;
}
