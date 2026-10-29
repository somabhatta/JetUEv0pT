// ============================================================================
//  Phase2.C -- combine the fine-Nch Phase1 moments of three samples into
//  v0(pT), with Poisson-bootstrap errors, and form sample differences.
//
//  Run   : root -l -b -q 'Phase2.C+O("<dir>","<tag0>,<tag1>,<tag2>","<out>.root")'
//          (or gSystem->Load("Phase2_C.so"); Phase2(...)). Compile with the same
//          -DV0_NBOOT as the Phase1 build that wrote the inputs, or larger:
//          replicas absent from the input are skipped.
//  Input : <dir>/<tag>_mean.root and <dir>/<tag>_corr.root (merged Phase1) for
//          each of the three tags; optional INPUT/cutNch_<tag>.txt (wide
//          centrality class in Nch, relative to the working directory).
//          The first tag is the reference of the stored differences.
//  Env   : V0_PTGRID=coarse -> output pT grid COARSE[] below; otherwise the
//          probe grid (Rule.h TREE_PTEDGE).
//
//  Method, per sample / Nch class / R / jet class / probe variant:
//   * moments per fine Nch bin j (event count N_j = entries of the reference
//     profile, mu_j = <[pT]_ref>, C2_j = <c2>, sigma_j = sqrt(max(C2_j,0))):
//       cov_j(Q)   = sum_{b in Q} ( <probe_b ref>_j - <probe_b>_j mu_j )
//       <probe>_j(Q) = sum_{b in Q} <probe_b>_j
//     for an output pT bin Q (fine probe bins whose centre lies in Q). Summing
//     cov and <probe> over the fine bins makes a coarse bin exact.
//   * average of ratios over the Nch range of the class (AVG_RATIO = true):
//       v0(Q)        = sum_j N_j cov_j / (<probe>_j sigma_j)  /  sum_j N_j
//       v0refnorm(Q) = sum_j N_j v0_j mu_j / sigma_j          /  sum_j N_j
//     (bins with <probe>_j <= 0 or sigma_j <= 0 do not enter the numerator).
//     This removes the bias from multiplicity fluctuations within the class.
//     AVG_RATIO = false would average the moments first, then take the ratio.
//   * matched reference per probe (see Rule.h):
//       f     : mean pRef_Nch,     C2 pC2_Nch
//       f_sub : mean pRefFsub_Nch, C2 pC2Sub_Nch
//       f_cs  : mean pRefFcs_Nch,  C2 pC2Sub_Nch
//   * bootstrap: the same combination is redone for each replica k = 1..NBOOT
//     with ALL four inputs taken from replica k (mean and corr). A replica with
//     any input missing is dropped (counted in nBoot_*). The error is CalErr()
//     of the replica values; the central value is the nominal (k = 0).
//   * differences d = v0(tag1 or tag2) - v0(tag0): error = spread of the
//     replica-wise differences (ReplicaDiffError). Replica k carries the same
//     per-event weights in every sample (seeded by the global event index), so
//     for samples built from the same events (tree variants of one sample) the
//     error includes the correlation; independent samples give the usual
//     quadrature-like error.
//
//  Output <dir>/<out>.root; <tag> = <sample>_c<ic>_R<i>_<SELN>, ic = 1 is the
//  inclusive class (ic = 0 only with a cutNch file), <var> = "", "_fsub", "_fcs":
//    v0<var>_<tag>             TH1D v0(pT), bootstrap error (analytic fallback)
//    v0refnorm<var>_<tag>      TH1D v0(pT)/v0_ref, v0_ref = sqrt(C2)/<[pT]_ref>
//    cov<var>_<tag>            TH1D class-averaged Cov(probe,[pT]_ref)
//    meanN<var>_<tag>          TH1D class-averaged <probe>
//    v0errAnalytic<var>_<tag>  TH1D analytic error of v0 (comparison only)
//    v0boot<var>_<tag>         TH2D v0 per (output pT bin, replica): y bin k+1 = replica k
//    v0refnormboot<var>_<tag>  TH2D same for v0refnorm
//    C2_<tag>, refbar_<tag>    TH1D (1 bin) class-averaged C2 and <[pT]_ref> (f reference)
//    nBoot_<tag>               TH1D (1 bin) number of replicas used
//    meanPtVsNch_<sample>_R<i>_<SELN>      copy of pRef_Nch (TProfile)
//    meanPtVsNchBoot_<sample>_R<i>_<SELN>  same as TH1D with bootstrap errors
//    dv0<var>_<tag1>m<tag0>_c<ic>_R<i>_<SELN>, dv0refnorm<var>_... (and tag2)
//    c_jet, c_jet_refnorm      quick-look canvases (also <out>_jet.pdf, <out>_refnorm_jet.pdf)
// ============================================================================
#include "Rule.h"
#include "Tool.C"
#include "TProfile2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TObjString.h"

static const bool AVG_RATIO = true;   // true: average of ratios over Nch (nominal); false: ratio of averages

// coarse output pT grid (V0_PTGRID=coarse)
static const int    NCOARSE = 8;
static const double COARSE[NCOARSE+1] = {0.2, 0.5, 1., 2., 3., 6., 12., 30., 100.};

// The four fine-Nch moment profiles of one (sample, R, jet class, probe, replica).
struct MomentBlock {
  TProfile2D *probe;       // <probe_b>(Nch,pT)          mean file
  TProfile2D *probeRef;    // <probe_b [pT]_ref>(Nch,pT) corr file
  TProfile   *refMean;     // <[pT]_ref>(Nch)            mean file (matched to the probe)
  TProfile   *c2;          // <c2>(Nch)                  corr file
};

// Fetch the block for replica k (k = 0: nominal names) of probe variant ivar.
// All four come from the same replica; a missing one makes the replica unusable.
static MomentBlock GetBlock(TFile* fMean, TFile* fCorr, int iR, int c, int k = 0, int ivar = IVAR_F){
  const char* vs      = VARSUF[ivar];
  const char* refStem = (ivar == IVAR_FSUB) ? "pRefFsub_Nch" : (ivar == IVAR_FCS) ? "pRefFcs_Nch" : "pRef_Nch";
  const char* c2Stem  = (ivar == IVAR_F) ? "pC2_Nch" : "pC2Sub_Nch";
  TString rep = (k == 0) ? TString("") : TString(Form("_boot%d", k));
  MomentBlock b;
  b.probe    = (TProfile2D*)fMean->Get(Form("pN_NchPt%s%s_R%d_%s",  vs, rep.Data(), iR, SELN[c]));
  b.refMean  = (TProfile*)  fMean->Get(Form("%s%s_R%d_%s", refStem, rep.Data(), iR, SELN[c]));
  b.probeRef = (TProfile2D*)fCorr->Get(Form("pNR_NchPt%s%s_R%d_%s", vs, rep.Data(), iR, SELN[c]));
  b.c2       = (TProfile*)  fCorr->Get(Form("%s%s_R%d_%s",  c2Stem, rep.Data(), iR, SELN[c]));
  return b;
}
static bool Complete(const MomentBlock& b){ return b.probe && b.probeRef && b.refMean && b.c2; }

// a - b with errors added in quadrature (replaced by ReplicaDiffError when replicas exist)
static TH1D* DiffHist(TH1D* a, TH1D* b, const char* name, const char* title){
  if (!a || !b) return 0;
  TH1D* d = (TH1D*)a->Clone(name); d->SetTitle(title);
  for (int i = 1; i <= a->GetNbinsX(); i++){
    d->SetBinContent(i, a->GetBinContent(i) - b->GetBinContent(i));
    d->SetBinError(i, TMath::Hypot(a->GetBinError(i), b->GetBinError(i)));
  }
  return d;
}

// Class averages over Nch bins [jlo,jhi] weighted by the event count N_j:
// C2bar = <C2> (clamped at 0), sigbar = sqrt(C2bar), refbar = <[pT]_ref>.
static void ClassAverageRef(const MomentBlock& b, int jlo, int jhi, double& C2bar, double& sigbar, double& refbar){
  double sumC2 = 0, sumMu = 0, sumN = 0;
  for (int j = jlo; j <= jhi; j++){
    double Nj = b.refMean->GetBinEntries(j); if (Nj <= 0) continue;
    sumC2 += Nj*b.c2->GetBinContent(j); sumMu += Nj*b.refMean->GetBinContent(j); sumN += Nj;
  }
  C2bar = (sumN > 0) ? sumC2/sumN : 0; if (C2bar < 0) C2bar = 0; sigbar = sqrt(C2bar);
  refbar = (sumN > 0) ? sumMu/sumN : 0;
}

struct BinResult { double v0, v0refnorm, cov, meanProbe, sumErr2, sumN; };

// v0 etc. for output pT bin [ptLo,ptHi) over Nch bins [jlo,jhi] (see header).
static void CombineNchBins(const MomentBlock& b, int jlo, int jhi, double ptLo, double ptHi, int nFinePt,
                           double C2bar, double sigbar, double refbar, BinResult& r){
  double sumCov = 0, sumProbe = 0, sumN = 0, sumErr2 = 0, sumV0 = 0, sumV0n = 0;
  for (int j = jlo; j <= jhi; j++){
    double Nj = b.refMean->GetBinEntries(j); if (Nj <= 0) continue;
    double muj = b.refMean->GetBinContent(j);
    double c2j = b.c2->GetBinContent(j); if (c2j < 0) c2j = 0; double sigj = sqrt(c2j);
    double covj = 0, probej = 0, err2j = 0;
    for (int ip = 1; ip <= nFinePt; ip++){
      double pt = b.probe->GetYaxis()->GetBinCenter(ip);
      if (pt < ptLo || pt >= ptHi) continue;
      double probe = b.probe->GetBinContent(j, ip), probeRef = b.probeRef->GetBinContent(j, ip);
      covj += probeRef - probe*muj; probej += probe;
      double e = b.probeRef->GetBinError(j, ip); err2j += e*e;
    }
    sumCov += Nj*covj; sumProbe += Nj*probej; sumN += Nj; sumErr2 += Nj*Nj*err2j;
    if (AVG_RATIO && probej > 0 && sigj > 0){
      double v0j = covj/(probej*sigj);
      sumV0 += Nj*v0j;
      if (muj > 0) sumV0n += Nj*(v0j*muj/sigj);   // v0_j / v0ref_j, v0ref_j = sigma_j/mu_j
    }
  }
  double meanProbe = (sumN > 0) ? sumProbe/sumN : 0, cov = (sumN > 0) ? sumCov/sumN : 0;
  double v0, v0n;
  if (AVG_RATIO){ v0 = (sumN > 0) ? sumV0/sumN : 0; v0n = (sumN > 0) ? sumV0n/sumN : 0; }
  else { v0  = (meanProbe > 0 && sigbar > 0) ? cov/(meanProbe*sigbar) : 0;
         v0n = (meanProbe > 0 && C2bar > 0 && refbar > 0) ? cov*refbar/(meanProbe*C2bar) : 0; }
  r.v0 = v0; r.v0refnorm = v0n; r.cov = cov; r.meanProbe = meanProbe; r.sumErr2 = sumErr2; r.sumN = sumN;
}

// Bootstrap spread of replica values (port of Phase3::cal_err):
//  - ignore non-finite and exactly-zero values (0 = quantity empty in that replica)
//  - population SD s about the mean of the rest
//  - drop values more than 3 s from that mean and return the population SD of
//    the remaining ones (s itself if fewer than 2 remain)
// Returns -1 with fewer than 2 usable values (caller falls back).
static double CalErr(const vector<double>& v){
  int n = 0; double m = 0;
  for (size_t i = 0; i < v.size(); i++) if (TMath::Finite(v[i]) && v[i] != 0.0){ m += v[i]; n++; }
  if (n < 2) return -1; m /= n;
  double s = 0;
  for (size_t i = 0; i < v.size(); i++) if (TMath::Finite(v[i]) && v[i] != 0.0){ double d = v[i] - m; s += d*d; }
  s = sqrt(s/n); if (s <= 0) return 0.;
  int n2 = 0; double m2 = 0;
  for (size_t i = 0; i < v.size(); i++) if (TMath::Finite(v[i]) && v[i] != 0.0 && fabs(v[i] - m) < 3*s){ m2 += v[i]; n2++; }
  if (n2 < 2) return s; m2 /= n2;
  double s2 = 0;
  for (size_t i = 0; i < v.size(); i++) if (TMath::Finite(v[i]) && v[i] != 0.0 && fabs(v[i] - m) < 3*s){ double d = v[i] - m2; s2 += d*d; }
  return sqrt(s2/n2);
}

// Error of d = a - b from replica-wise differences a_k - b_k (k >= 1, both
// non-zero and finite): 3-sigma-trimmed population SD as in CalErr, but exact
// zeros are kept (a difference of present replicas can be exactly 0).
static void ReplicaDiffError(TH1D* d, TH2D* aBoot, TH2D* bBoot){
  if (!d || !aBoot || !bBoot) return;
  for (int q = 1; q <= d->GetNbinsX(); q++){
    vector<double> v;
    for (int k = 2; k <= aBoot->GetNbinsY(); k++){   // y bin k = replica k-1; bin 1 = nominal
      double x = aBoot->GetBinContent(q, k), y = bBoot->GetBinContent(q, k);
      if (x != 0.0 && y != 0.0 && TMath::Finite(x) && TMath::Finite(y)) v.push_back(x - y);
    }
    int n = v.size(); if (n < 2) continue;
    double m = 0; for (double x : v) m += x; m /= n;
    double sd = 0; for (double x : v) sd += (x - m)*(x - m); sd = sqrt(sd/n);
    double m2 = 0, s2 = 0; int n2 = 0;
    for (double x : v) if (sd == 0 || fabs(x - m) < 3*sd){ m2 += x; n2++; }
    if (n2 >= 2){ m2 /= n2; for (double x : v) if (sd == 0 || fabs(x - m) < 3*sd) s2 += (x - m2)*(x - m2); sd = sqrt(s2/n2); }
    d->SetBinError(q, sd);
  }
}

void Phase2(const char* dir = "OUTPUT", const char* samples = "OO_vac_nom,OO_NR_nom,OO_RC_cs10", const char* outname = "v0pt_final.root"){
  gStyle->SetOptStat(0);
  TString sampleStr[3];
  { TString all(samples); TObjArray* tok = all.Tokenize(",");
    if (tok->GetEntries() != 3){ printf("Phase2: need exactly 3 comma-separated samples, got '%s'\n", samples); return; }
    for (int i = 0; i < 3; i++) sampleStr[i] = ((TObjString*)tok->At(i))->GetString(); }
  const char* S[3] = {sampleStr[0].Data(), sampleStr[1].Data(), sampleStr[2].Data()};
  const int color[3] = {kBlack, kRed+1, kAzure+2};

  // output pT grid
  const bool coarse = getenv("V0_PTGRID") && TString(getenv("V0_PTGRID")) == "coarse";
  const int nQ = coarse ? NCOARSE : NCEDGE-1;
  double qEdges[NCEDGE]; for (int i = 0; i <= nQ; i++) qEdges[i] = coarse ? COARSE[i] : CEDGES[i];

  TFile *fMean[3], *fCorr[3];
  for (int s = 0; s < 3; s++){
    fMean[s] = TFile::Open(Form("%s/%s_mean.root", dir, S[s]));
    fCorr[s] = TFile::Open(Form("%s/%s_corr.root", dir, S[s]));
    if (!fMean[s] || !fCorr[s] || fMean[s]->IsZombie() || fCorr[s]->IsZombie()){ printf("missing merged mean/corr for %s\n", S[s]); return; }
  }
  // axes from the first sample's nominal f profile
  TProfile2D* axisRef = (TProfile2D*)fMean[0]->Get(Form("pN_NchPt_R0_%s", SELN[ISEL_JET]));
  if (!axisRef){ printf("Phase2: pN_NchPt_R0_%s missing in %s mean\n", SELN[ISEL_JET], S[0]); return; }
  int nFinePt = axisRef->GetNbinsY(), nNchBins = axisRef->GetNbinsX();
  int nR = 0; while (fMean[0]->Get(Form("pN_NchPt_R%d_%s", nR, SELN[ISEL_JET]))) nR++;
  printf("Phase2: nR=%d fine-pT=%d Nch-bins=%d coarse-pT=%d NCENT=%d  AVG_RATIO=%d  NBOOT=%d\n",
         nR, nFinePt, nNchBins, nQ, NCENT, (int)AVG_RATIO, NBOOT);

  TFile* fout = new TFile(Form("%s/%s", dir, outname), "RECREATE");
  // results kept for the differences: [variant][sample][Nch class][R][jet class]
  TH1D* hV0       [NVAR][3][NCENTX][RMAX][NSEL];
  TH1D* hV0n      [NVAR][3][NCENTX][RMAX][NSEL];
  TH2D* hV0Boot   [NVAR][3][NCENTX][RMAX][NSEL];
  TH2D* hV0nBoot  [NVAR][3][NCENTX][RMAX][NSEL];
  for (int iv = 0; iv < NVAR; iv++) for (int s = 0; s < 3; s++) for (int ic = 0; ic < NCENTX; ic++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++){
    hV0[iv][s][ic][iR][c] = 0; hV0n[iv][s][ic][iR][c] = 0; hV0Boot[iv][s][ic][iR][c] = 0; hV0nBoot[iv][s][ic][iR][c] = 0;
  }

  for (int s = 0; s < 3; s++){
    // optional wide class: INPUT/cutNch_<sample>.txt, Nch edges in DESCENDING
    // order, class ic = [cut[ic+1], cut[ic]) (ic = 0 most central)
    double cut[NCENT+1];
    bool haveCut = Tool::readCuts(Form("INPUT/cutNch_%s.txt", S[s]), cut, NCENT+1);
    if (!haveCut) printf("  WARN: no Nch cuts for %s -> inclusive only\n", S[s]);
    for (int ic = 0; ic < NCENTX; ic++){
      if (ic != ICENT_INCL && !haveCut) continue;
      int jlo, jhi;   // Nch bin range of the class
      if (ic == ICENT_INCL){ jlo = 1; jhi = nNchBins; }
      else {
        jlo = axisRef->GetXaxis()->FindBin(cut[ic+1]);
        jhi = axisRef->GetXaxis()->FindBin(cut[ic]) - 1;
        if (jlo < 1) jlo = 1; if (jhi > nNchBins) jhi = nNchBins;
      }

      for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++){
        MomentBlock fRef = GetBlock(fMean[s], fCorr[s], iR, c, 0, IVAR_F);
        if (!fRef.refMean || !fRef.c2) continue;
        TString tag = Form("%s_c%d_R%d_%s", S[s], ic, iR, SELN[c]);

        for (int ivar = 0; ivar < NVAR; ivar++){
          MomentBlock nominal = GetBlock(fMean[s], fCorr[s], iR, c, 0, ivar);
          if (!Complete(nominal)) continue;
          double C2bar, sigbar, refbar; ClassAverageRef(nominal, jlo, jhi, C2bar, sigbar, refbar);

          const char* vs = VARSUF[ivar]; const char* vl = VARLAB[ivar];
          TH1D* hv0   = new TH1D(Form("v0%s_%s", vs, tag.Data()), Form("%s c%d R%d %s (%s);p_{T}^{probe} (GeV/c);v_{0}", S[s], ic, iR, SELN[c], vl), nQ, qEdges);
          TH1D* hv0n  = new TH1D(Form("v0refnorm%s_%s", vs, tag.Data()), ";p_{T};v_{0}/v_{0}^{ref}", nQ, qEdges);
          TH1D* hcov  = new TH1D(Form("cov%s_%s", vs, tag.Data()), Form(";p_{T};cov(%s,ref)", vl), nQ, qEdges);
          TH1D* hmean = new TH1D(Form("meanN%s_%s", vs, tag.Data()), Form(";p_{T};<%s(p_{T})>", vl), nQ, qEdges);
          TH1D* hv0eA = new TH1D(Form("v0errAnalytic%s_%s", vs, tag.Data()), ";p_{T};analytic err(v0)", nQ, qEdges);
          TH2D* hv0b  = new TH2D(Form("v0boot%s_%s", vs, tag.Data()), ";p_{T} bin;replica", nQ, 0, nQ, NBOOT+1, -0.5, NBOOT+0.5);
          TH2D* hv0nb = new TH2D(Form("v0refnormboot%s_%s", vs, tag.Data()), ";p_{T} bin;replica", nQ, 0, nQ, NBOOT+1, -0.5, NBOOT+0.5);

          // bootstrap replicas (C2 and refbar replicas: f reference only)
          vector<vector<double> > repV0(nQ), repV0n(nQ), repCov(nQ), repMean(nQ);
          vector<double> repC2, repRef;
          int nBootUsed = 0;
          for (int k = 1; k <= NBOOT; k++){
            MomentBlock rep = GetBlock(fMean[s], fCorr[s], iR, c, k, ivar);
            if (!Complete(rep)) continue;
            double C2k, sigk, refk; ClassAverageRef(rep, jlo, jhi, C2k, sigk, refk);
            nBootUsed++;
            if (ivar == IVAR_F){ repC2.push_back(C2k); repRef.push_back(refk); }
            for (int q = 0; q < nQ; q++){
              BinResult r; CombineNchBins(rep, jlo, jhi, qEdges[q], qEdges[q+1], nFinePt, C2k, sigk, refk, r);
              repV0[q].push_back(r.v0); repV0n[q].push_back(r.v0refnorm);
              hv0b->SetBinContent(q+1, k+1, r.v0); hv0nb->SetBinContent(q+1, k+1, r.v0refnorm);
              repCov[q].push_back(r.cov); repMean[q].push_back(r.meanProbe);
            }
          }

          // nominal values with bootstrap errors
          for (int q = 0; q < nQ; q++){
            BinResult r; CombineNchBins(nominal, jlo, jhi, qEdges[q], qEdges[q+1], nFinePt, C2bar, sigbar, refbar, r);
            double v0errAnalytic = (r.meanProbe > 0 && sigbar > 0 && r.sumN > 0) ? sqrt(r.sumErr2)/r.sumN/(r.meanProbe*sigbar) : 0;
            double v0err   = CalErr(repV0[q]);
            double v0nerr  = CalErr(repV0n[q]);
            double coverr  = CalErr(repCov[q]);
            double meanerr = CalErr(repMean[q]);
            hv0->SetBinContent(q+1, r.v0); hv0b->SetBinContent(q+1, 1, r.v0); hv0nb->SetBinContent(q+1, 1, r.v0refnorm);
            hv0->SetBinError(q+1, (v0err >= 0) ? v0err : v0errAnalytic);
            hv0n->SetBinContent(q+1, r.v0refnorm);
            hv0n->SetBinError(q+1, (v0nerr >= 0) ? v0nerr : 0.);
            hcov->SetBinContent(q+1, r.cov);
            hcov->SetBinError(q+1, (coverr >= 0) ? coverr : 0.);
            hmean->SetBinContent(q+1, r.meanProbe);
            hmean->SetBinError(q+1, (meanerr >= 0) ? meanerr : 0.);
            hv0eA->SetBinContent(q+1, v0errAnalytic);
          }
          fout->cd(); hv0->Write(); hv0n->Write(); hcov->Write(); hmean->Write(); hv0eA->Write();
          hv0b->Write(); hv0nb->Write();

          if (ivar == IVAR_F){   // reference summaries, written once per tag
            double c2err = CalErr(repC2), referr = CalErr(repRef);
            TH1D* hC2 = new TH1D("C2_"+tag, ";;C2=Var([p_{T}]_{ref})", 1, 0, 1);
            hC2->SetBinContent(1, C2bar); hC2->SetBinError(1, (c2err >= 0) ? c2err : 0.); hC2->Write();
            TH1D* hRef = new TH1D("refbar_"+tag, ";;<[p_{T}]_{ref}>", 1, 0, 1);
            hRef->SetBinContent(1, refbar); hRef->SetBinError(1, (referr >= 0) ? referr : 0.); hRef->Write();
            TH1D* hNb = new TH1D("nBoot_"+tag, ";;#replicas used", 1, 0, 1); hNb->SetBinContent(1, nBootUsed); hNb->Write();
          }
          hV0[ivar][s][ic][iR][c] = hv0; hV0n[ivar][s][ic][iR][c] = hv0n;
          hV0Boot[ivar][s][ic][iR][c] = hv0b; hV0nBoot[ivar][s][ic][iR][c] = hv0nb;
        }
      }
    }

    // <[pT]_ref> vs fine Nch: passthrough (TProfile) and a TH1D copy whose bin
    // errors are the bootstrap spread over the pRef_Nch_boot<k> replicas
    // (profile error where fewer than 2 replicas have entries)
    for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++){
      TProfile* refVsNch = (TProfile*)fMean[s]->Get(Form("pRef_Nch_R%d_%s", iR, SELN[c]));
      if (!refVsNch) continue;
      fout->cd(); refVsNch->Clone(Form("meanPtVsNch_%s_R%d_%s", S[s], iR, SELN[c]))->Write();
      int nb = refVsNch->GetNbinsX();
      TH1D* hBoot = refVsNch->ProjectionX(Form("meanPtVsNchBoot_%s_R%d_%s", S[s], iR, SELN[c]));
      hBoot->SetTitle(";N_{ch};#LT[p_{T}]_{ref}#GT");
      vector<vector<double> > repBin(nb+1);
      for (int k = 1; k <= NBOOT; k++){
        TProfile* rep = (TProfile*)fMean[s]->Get(Form("pRef_Nch_boot%d_R%d_%s", k, iR, SELN[c]));
        if (!rep) continue;
        for (int j = 1; j <= nb; j++) if (rep->GetBinEntries(j) > 0) repBin[j].push_back(rep->GetBinContent(j));
      }
      for (int j = 1; j <= nb; j++){
        hBoot->SetBinContent(j, refVsNch->GetBinContent(j));
        double sd = CalErr(repBin[j]); hBoot->SetBinError(j, (sd >= 0) ? sd : refVsNch->GetBinError(j));
      }
      hBoot->Write();
    }
  }

  // differences to the first sample, replica-wise errors
  for (int iv = 0; iv < NVAR; iv++) for (int ic = 0; ic < NCENTX; ic++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++){
    const char* vs = VARSUF[iv];
    for (int s = 1; s <= 2; s++){
      const char* lab = (s == 1) ? "NR-vac" : "RC-vac";
      TH1D* d = DiffHist(hV0[iv][s][ic][iR][c], hV0[iv][0][ic][iR][c],
                         Form("dv0%s_%sm%s_c%d_R%d_%s", vs, S[s], S[0], ic, iR, SELN[c]),
                         Form("%s c%d R%d %s;p_{T};#Deltav_{0}", lab, ic, iR, SELN[c]));
      ReplicaDiffError(d, hV0Boot[iv][s][ic][iR][c], hV0Boot[iv][0][ic][iR][c]);
      if (d) d->Write();
    }
    for (int s = 1; s <= 2; s++){
      const char* lab = (s == 1) ? "NR-vac" : "RC-vac";
      TH1D* d = DiffHist(hV0n[iv][s][ic][iR][c], hV0n[iv][0][ic][iR][c],
                         Form("dv0refnorm%s_%sm%s_c%d_R%d_%s", vs, S[s], S[0], ic, iR, SELN[c]),
                         Form("%s c%d R%d %s;p_{T};#Delta(v_{0}/v_{0})", lab, ic, iR, SELN[c]));
      ReplicaDiffError(d, hV0nBoot[iv][s][ic][iR][c], hV0nBoot[iv][0][ic][iR][c]);
      if (d) d->Write();
    }
  }

  // quick-look canvases: f probe, inclusive class, jet class ISEL_JET
  int ic = ICENT_INCL, c = ISEL_JET;
  TString outStem = TString(outname).ReplaceAll(".root", "");
  TCanvas* cv = new TCanvas("c_jet", "", 1500, 900); cv->Divide(3, 2);
  for (int iR = 0; iR < nR; iR++){
    cv->cd(iR+1); gPad->SetLogx();
    TLegend* lg = new TLegend(0.55, 0.72, 0.88, 0.88);
    for (int s = 0; s < 3; s++){
      TH1D* h = hV0[IVAR_F][s][ic][iR][c]; if (!h) continue;
      h->SetLineColor(color[s]); h->SetMarkerColor(color[s]); h->SetMarkerStyle(20+s);
      h->GetYaxis()->SetRangeUser(-0.12, 0.15);
      h->Draw(s == 0 ? "E1" : "E1 SAME"); lg->AddEntry(h, S[s], "lp");
    }
    lg->Draw();
  }
  cv->SaveAs(Form("%s/%s_jet.pdf", dir, outStem.Data())); cv->Write();

  TCanvas* cvn = new TCanvas("c_jet_refnorm", "", 1500, 900); cvn->Divide(3, 2);
  for (int iR = 0; iR < nR; iR++){
    cvn->cd(iR+1); gPad->SetLogx();
    double ylo = 1e30, yhi = -1e30;
    for (int s = 0; s < 3; s++){
      TH1D* h = hV0n[IVAR_F][s][ic][iR][c]; if (!h) continue;
      for (int b = 1; b <= h->GetNbinsX(); b++){
        double v = h->GetBinContent(b), e = h->GetBinError(b);
        if (v == 0 && e == 0) continue; if (v < ylo) ylo = v; if (v > yhi) yhi = v;
      }
    }
    if (ylo > yhi){ ylo = -1; yhi = 1; }
    double pad = 0.25*(yhi - ylo + 1e-9); ylo -= pad; yhi += pad;
    TLegend* lg = new TLegend(0.55, 0.72, 0.88, 0.88);
    for (int s = 0; s < 3; s++){
      TH1D* h = hV0n[IVAR_F][s][ic][iR][c]; if (!h) continue;
      h->SetLineColor(color[s]); h->SetMarkerColor(color[s]); h->SetMarkerStyle(20+s);
      h->GetYaxis()->SetRangeUser(ylo, yhi);
      h->Draw(s == 0 ? "E1" : "E1 SAME"); lg->AddEntry(h, S[s], "lp");
    }
    lg->Draw();
  }
  cvn->SaveAs(Form("%s/%s_refnorm_jet.pdf", dir, outStem.Data())); cvn->Write();

  fout->Write(); fout->Close();
  printf("DONE -> %s/%s  (v0, v0refnorm, cov, meanN, v0boot, C2, refbar, nBoot, meanPtVsNch[Boot], dv0, dv0refnorm)\n", dir, outname);
}
