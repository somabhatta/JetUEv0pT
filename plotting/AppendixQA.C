// ============================================================================
//  AppendixQA.C -- cross-check figures for the appendix of the jet--UE v0(pT)
//  paper.
//
//  Two kinds of input:
//   * tree-level checks (qa_jets, qa_spectra) read the first <nfiles> files of
//     the jetv0ana trees, <treeDir>/<sample>/*.root, tree t_<variant>. The same
//     generated events exist with and without embedded UE, so jet-level
//     subtraction can be closed event by event. Both are skipped when
//     <treeDir> does not exist.
//   * observable-level checks read the Phase2 outputs of the first (coarse)
//     production, <resultsDir>/<ref>/acc<i>/v0pt_<file>.root with
//     ref = refnom (reference outside jet + perpendicular cones) or refgap
//     (reference at |Delta eta| > R + 0.5). Objects used (TH1D, x = probe pT
//     unless stated):
//       v0<probe>_<tag>_c1_R<iR>_<jetClass>      v0(pT), bootstrap errors
//       meanN<probe>_<tag>_c1_R<iR>_<jetClass>   <f(pT)>
//       refbar_<tag>_c1_R<iR>_<jetClass>         <[pT]_ref> (bin 1)
//       C2_<tag>_c1_R<iR>_<jetClass>             C2 of the reference (bin 1)
//     <probe> = "" (raw f), "_fsub" (perp. cone subtracted), "_fcs" (constituent
//     subtracted); iR = 0..3 for R = 0.1, 0.2, 0.4, 0.8.
//  Acceptance: R <= 0.4 with |eta_trk| < 0.9 (acc0), R = 0.8 with |eta_trk| < 2.5 (acc1).
//
//  Figures (PDF and PNG) in <outDir>:
//   qa_jets       anti-kT jets: UE-corrected jet pT spectra, area/(pi R^2), and
//                 the event-matched residual delta pT (perp-cone rho*A and
//                 constituent subtraction) for PYTHIA + UE
//   qa_spectra    in-cone spectra: truth n vs embedded raw n, n - b, n_cs (closure)
//   qa_reference  <[pT]_ref> and C2,ref vs R for several samples and UE levels
//   qa_fmean      <f(pT)> for raw / perp / CS probes, pp and O+O
//   qa_closure    v0(pT) PYTHIA: no UE vs +UE (MB, central), f_sub and f_cs
//   qa_recoilsub  Delta v0(pT): thermal-subtraction variants (O+O recoils)
//   qa_refgap     signal vs R: nominal vs eta-gap reference
//   qa_annulus    R-differential signal S(R2) - S(R1)
//
//  Run: root -l -b -q 'AppendixQA.C+("/dcache/alice/sbhatta/JewelV0pT/trees","/data/alice/sbhatta/JewelV0pT/results_coarse","qa",20)'
// ============================================================================
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include <algorithm>
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TMath.h"
#include <vector>

static TString gTreeDir, gResultsDir, gOutDir;
static int gNumFiles = 20;     // tree files read per sample
static int gUniqueId = 0;      // suffix for unique object names
static const int kNumRadii = 4, kNumPtBins = 25, kNumJetDefs = 8;   // jet definitions stored in the trees
static const double kRadius[kNumRadii] = {0.1, 0.2, 0.4, 0.8};
static const int kAcceptanceOfRadius[kNumRadii] = {0, 0, 0, 1};
static const double kPtEdges[kNumPtBins + 1] = {0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10., 12., 15., 20., 25., 30., 40., 50., 60., 80., 100., 150., 300.};
static const int kCurveColor[6] = {kBlack, kRed + 1, kAzure + 2, kGreen + 2, kOrange + 7, kViolet + 1};
static const int kCurveMarker[6] = {20, 21, 22, 24, 25, 26};
static TLatex gText;

// probe-pT ranges of the range averages: soft and hard probes
static const double kProbeRange[2][2] = {{1, 6}, {12, 100}};
static const char* kProbeRangeLabel[2] = {"soft 1-6 GeV/#it{c}", "hard 12-100 GeV/#it{c}"};

// ---------------------------------------------------------------------------
// trees

// chain of tree t_<variant> over the first gNumFiles files (sorted by name) in <treeDir>/<sample>
static TChain* OpenTreeChain(const char* sample, const char* variant)
{
  TChain* c = new TChain(Form("t_%s", variant));
  TString dirName = gTreeDir + "/" + sample;
  void* dir = gSystem->OpenDirectory(dirName); const char* entry; int n = 0;
  std::vector<TString> files;
  while (dir && (entry = gSystem->GetDirEntry(dir))) { TString s(entry); if (s.EndsWith(".root")) files.push_back(dirName + "/" + s); }
  if (dir) gSystem->FreeDirectory(dir);
  std::sort(files.begin(), files.end());
  for (auto& f : files) { if (n++ >= gNumFiles) break; c->Add(f); }
  return c;
}

// per-event branches; arrays are indexed by jet definition k (and probe-pT bin b: k * kNumPtBins + b)
struct EventBranches {
  float weight, jetPt[kNumJetDefs], jetPtCorr[kNumJetDefs], jetArea[kNumJetDefs], rho[kNumJetDefs], csJetPt[kNumJetDefs];
  float nIn[kNumJetDefs * kNumPtBins];     // particles in the cone
  float nPerp[kNumJetDefs * kNumPtBins];   // particles in the perpendicular cones (background b)
  float nCS[kNumJetDefs * kNumPtBins];     // particles after constituent subtraction
};
static void BindBranches(TChain* c, EventBranches& e)
{
  c->SetBranchStatus("*", 0);
  for (auto b : {"weight", "jetPt", "jetPtCorr", "jetArea", "rho", "csJetPt", "nIn", "nPerp", "nCS"}) c->SetBranchStatus(b, 1);
  c->SetBranchAddress("weight", &e.weight); c->SetBranchAddress("jetPt", e.jetPt); c->SetBranchAddress("jetPtCorr", e.jetPtCorr);
  c->SetBranchAddress("jetArea", e.jetArea); c->SetBranchAddress("rho", e.rho); c->SetBranchAddress("csJetPt", e.csJetPt);
  c->SetBranchAddress("nIn", e.nIn); c->SetBranchAddress("nPerp", e.nPerp); c->SetBranchAddress("nCS", e.nCS);
}

// ---------------------------------------------------------------------------
// drawing helpers

static void SetPadMargins(double left = 0.15)
{ gPad->SetTicks(1, 1); gPad->SetLeftMargin(left); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); }

// colour/marker number i of the palette; lineOnly hides the markers
static void StyleHist(TH1* h, int i, bool lineOnly = false)
{
  h->SetLineColor(kCurveColor[i]); h->SetMarkerColor(kCurveColor[i]); h->SetMarkerStyle(kCurveMarker[i]); h->SetMarkerSize(0.9); h->SetLineWidth(2);
  h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05); h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
  if (lineOnly) h->SetMarkerSize(0);
}
static void StyleGraph(TGraphErrors* g, int color, int marker)
{ g->SetLineColor(color); g->SetMarkerColor(color); g->SetMarkerStyle(marker); g->SetLineWidth(2); }

static TLegend* MakeLegend(double x0, double y0, double x1, double y1)
{ TLegend* l = new TLegend(x0, y0, x1, y1); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.042); return l; }

// empty axis frame; a dashed line at y = 0 unless the y axis is logarithmic
static TH1D* DrawFrame(double xMin, double xMax, double yMin, double yMax, const char* xTitle, const char* yTitle, bool logX = false, bool logY = false)
{
  TH1D* f = new TH1D(Form("fr%d", gUniqueId++), "", 10, xMin, xMax); f->SetDirectory(0); f->GetYaxis()->SetRangeUser(yMin, yMax);
  f->GetXaxis()->SetTitle(xTitle); f->GetYaxis()->SetTitle(yTitle); StyleHist(f, 0); f->Draw("AXIS");
  if (logX) { gPad->SetLogx(); f->GetXaxis()->SetMoreLogLabels(); f->GetXaxis()->SetNoExponent(); }
  if (logY) gPad->SetLogy();
  if (!logY) { TLine* zero = new TLine(xMin, 0, xMax, 0); zero->SetLineStyle(2); zero->SetLineColor(kGray + 1); zero->Draw(); }
  return f;
}

// y range that follows the drawn points: the data fill the pad between a bottom
// and a top fraction reserved for labels and legends; zero is always included.
// growExtent widens [lo, hi] to the drawn points (bins without an estimate are ignored).
static void growExtent(TH1* h, double& lo, double& hi)
{
  if (!h) return;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double v = h->GetBinContent(i), e = h->GetBinError(i);
    if (e <= 0 || v < -900) continue;
    lo = std::min(lo, v - e); hi = std::max(hi, v + e);
  }
}
static void growExtent(TGraphErrors* g, double& lo, double& hi)
{
  if (!g) return;
  for (int i = 0; i < g->GetN(); i++) {
    double v = g->GetY()[i], e = g->GetEY()[i];
    lo = std::min(lo, v - e); hi = std::max(hi, v + e);
  }
}
static void fitYRange(TH1* frame, double lo, double hi, double topFrac, double bottomFrac)
{
  if (lo > hi) return;   // nothing drawn: keep the default range
  lo = std::min(lo, 0.0); hi = std::max(hi, 0.0);
  double span = hi - lo; if (span <= 0) span = 0.1;
  double total = span / (1.0 - topFrac - bottomFrac);
  frame->GetYaxis()->SetRangeUser(lo - bottomFrac * total, hi + topFrac * total); gPad->Modified();
}

// ---------------------------------------------------------------------------
// Phase2 results

// one histogram from <resultsDir>/<ref>/acc<acc>/v0pt_<file>.root (a detached copy), or 0
static TH1D* LoadResult(int acc, const char* file, const char* name, const char* ref = "refnom")
{
  TFile* f = TFile::Open(Form("%s/%s/acc%d/v0pt_%s.root", gResultsDir.Data(), ref, acc, file));
  if (!f || f->IsZombie()) return 0;
  TH1D* h = (TH1D*)f->Get(name);
  TH1D* copy = h ? (TH1D*)h->Clone(Form("r%d", gUniqueId++)) : 0;
  if (copy) copy->SetDirectory(0);
  f->Close();
  return copy;
}

// a - b per bin; errors added in quadrature (0 if either error is missing)
static TH1D* Difference(TH1D* a, TH1D* b)
{
  if (!a || !b) return 0;
  TH1D* d = (TH1D*)a->Clone(Form("d%d", gUniqueId++)); d->SetDirectory(0);
  for (int i = 1; i <= a->GetNbinsX(); i++) {
    d->SetBinContent(i, a->GetBinContent(i) - b->GetBinContent(i));
    d->SetBinError(i, (a->GetBinError(i) > 0 && b->GetBinError(i) > 0) ? hypot(a->GetBinError(i), b->GetBinError(i)) : 0);
  }
  return d;
}

// hide bins without an estimate (no error: no events); content -999, below every frame
static TH1D* HideUnreliableBins(TH1D* h)
{
  if (!h) return 0;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinError(i) <= 0) { h->SetBinContent(i, -999); h->SetBinError(i, 0); }
  return h;
}

// Range average over probe pT: inverse-variance weighted mean of the bins whose
// lower edge lies in [ptLow, ptHigh); bins without an estimate (no error) are
// skipped. Returns false if no bin contributes.
static bool RangeAverage(TH1D* h, double ptLow, double ptHigh, double& mean, double& error)
{
  if (!h) return false;
  double sumW = 0, sumWX = 0;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double lowEdge = h->GetXaxis()->GetBinLowEdge(i), err = h->GetBinError(i);
    if (lowEdge < ptLow - 1e-6 || lowEdge >= ptHigh - 1e-6 || err <= 0) continue;
    sumW += 1 / (err * err); sumWX += h->GetBinContent(i) / (err * err);
  }
  if (sumW <= 0) return false;
  mean = sumWX / sumW; error = 1 / sqrt(sumW); return true;
}

// signal S(R) = recoils - vacuum, both with O+O MB UE embedded (f_sub, jets 40-60),
// for radius index iR and reference ref
static TH1D* RecoilSignal(int iR, const char* ref = "refnom")
{
  return Difference(LoadResult(kAcceptanceOfRadius[iR], "OO_embMB", Form("v0_fsub_OO_RC_cs10_embMB_c1_R%d_jet40_60", iR), ref),
                    LoadResult(kAcceptanceOfRadius[iR], "OO_embMB", Form("v0_fsub_OO_vac_embMB_c1_R%d_jet40_60", iR), ref));
}

// ---------------------------------------------------------------------------
// qa_jets: jet pT spectra, jet area, and event-matched residual delta pT.
// Top row R = 0.4 (jet definition k = 2), bottom row R = 0.8 (k = 7).
static void PlotJetQA()
{
  TCanvas c("cj", "", 1500, 1000); c.Divide(3, 2, 0.001, 0.001);
  int jetDefOfRow[2] = {2, 7}; const char* rowLabel[2] = {"R = 0.4, |#eta_{trk}|<0.9", "R = 0.8, |#eta_{trk}|<2.5"};
  struct Sample { const char* sample; const char* variant; const char* label; };
  Sample samples[4] = {{"pp_vac", "nom", "JEWEL pp"}, {"py_pp", "nom", "PYTHIA pp"}, {"OO_RC", "cs10", "JEWEL O+O recoils"}, {"OO_RC", "cs10_embMB", "O+O recoils + UE(MB)"}};
  for (int row = 0; row < 2; row++) {
    int k = jetDefOfRow[row], iR = (k == 2) ? 2 : 3; double R = kRadius[iR];
    TH1D* hJetPt[4]; TH1D* hArea[4];
    for (int s = 0; s < 4; s++) {
      hJetPt[s] = new TH1D(Form("hpt%d%d", row, s), ";#it{p}_{T}^{jet,corr} (GeV/#it{c});normalized", 60, 0, 150);
      hArea[s] = new TH1D(Form("har%d%d", row, s), ";A_{jet}/(#pi R^{2});normalized", 60, 0, 1.5);
      TChain* ch = OpenTreeChain(samples[s].sample, samples[s].variant); EventBranches e; BindBranches(ch, e);
      for (Long64_t i = 0; i < ch->GetEntries(); i++) {
        ch->GetEntry(i); if (e.jetPt[k] <= 0) continue;   // no jet in this event
        hJetPt[s]->Fill(e.jetPtCorr[k], e.weight);
        if (e.jetPtCorr[k] > 20) hArea[s]->Fill(e.jetArea[k] / (TMath::Pi() * R * R), e.weight);
      }
      if (hJetPt[s]->Integral() > 0) hJetPt[s]->Scale(1. / hJetPt[s]->Integral());
      if (hArea[s]->Integral() > 0) hArea[s]->Scale(1. / hArea[s]->Integral());
      StyleHist(hJetPt[s], s, true); StyleHist(hArea[s], s, true); delete ch;
    }
    // (a, d) jet pT spectra
    c.cd(row * 3 + 1); SetPadMargins(); gPad->SetLogy(); hJetPt[0]->GetYaxis()->SetRangeUser(1e-5, 1); hJetPt[0]->Draw("HIST");
    TLegend* legendPt = MakeLegend(0.45, 0.66, 0.97, 0.93);
    for (int s = 0; s < 4; s++) { hJetPt[s]->Draw("HIST SAME"); legendPt->AddEntry(hJetPt[s], samples[s].label, "l"); }
    legendPt->Draw();
    gText.SetTextSize(0.045); gText.DrawLatex(0.2, 0.2, rowLabel[row]);
    // (b, e) jet area for jets above 20 GeV/c
    c.cd(row * 3 + 2); SetPadMargins(); hArea[0]->GetYaxis()->SetRangeUser(0, 1.3 * hArea[0]->GetMaximum()); hArea[0]->Draw("HIST");
    for (int s = 1; s < 4; s++) hArea[s]->Draw("HIST SAME");
    gText.DrawLatex(0.2, 0.86, Form("%s, #it{p}_{T}^{jet,corr}>20 GeV/#it{c}", rowLabel[row]));
    // (c, f) residual delta pT: the same PYTHIA events without UE (truth) and with
    // UE embedded (MB, 0-5%); subtracted jet pT minus truth jet pT, truth jets > 20 GeV/c
    c.cd(row * 3 + 3); SetPadMargins();
    TChain* tTruth = OpenTreeChain("py_pp", "nom"); TChain* tEmbMB = OpenTreeChain("py_pp", "embMB"); TChain* tEmbC = OpenTreeChain("py_pp", "embC");
    EventBranches eTruth, eEmbMB, eEmbC; BindBranches(tTruth, eTruth); BindBranches(tEmbMB, eEmbMB); BindBranches(tEmbC, eEmbC);
    TH1D* hDeltaPt[4]; const char* deltaLabel[4] = {"#rho A, UE MB", "CS, UE MB", "#rho A, UE 0-5%", "CS, UE 0-5%"};
    for (int j = 0; j < 4; j++) {
      hDeltaPt[j] = new TH1D(Form("hd%d%d", row, j), ";#delta#it{p}_{T} = #it{p}_{T}^{sub}(+UE) - #it{p}_{T}(no UE) (GeV/#it{c});normalized", 80, -20, 20);
      StyleHist(hDeltaPt[j], j, true);
    }
    Long64_t n = std::min(tTruth->GetEntries(), std::min(tEmbMB->GetEntries(), tEmbC->GetEntries()));
    for (Long64_t i = 0; i < n; i++) {
      tTruth->GetEntry(i); tEmbMB->GetEntry(i); tEmbC->GetEntry(i);
      if (eTruth.jetPt[k] < 20) continue;
      if (eEmbMB.jetPt[k] > 0) { hDeltaPt[0]->Fill(eEmbMB.jetPtCorr[k] - eTruth.jetPt[k], eTruth.weight); hDeltaPt[1]->Fill(eEmbMB.csJetPt[k] - eTruth.jetPt[k], eTruth.weight); }
      if (eEmbC.jetPt[k] > 0) { hDeltaPt[2]->Fill(eEmbC.jetPtCorr[k] - eTruth.jetPt[k], eTruth.weight); hDeltaPt[3]->Fill(eEmbC.csJetPt[k] - eTruth.jetPt[k], eTruth.weight); }
    }
    TLegend* legendDelta = MakeLegend(0.17, 0.62, 0.62, 0.93);
    double yMax = 0;
    for (int j = 0; j < 4; j++) { if (hDeltaPt[j]->Integral() > 0) hDeltaPt[j]->Scale(1. / hDeltaPt[j]->Integral()); yMax = std::max(yMax, hDeltaPt[j]->GetMaximum()); }
    for (int j = 0; j < 4; j++) {
      hDeltaPt[j]->GetYaxis()->SetRangeUser(0, 1.35 * yMax); hDeltaPt[j]->Draw(j ? "HIST SAME" : "HIST");
      legendDelta->AddEntry(hDeltaPt[j], Form("%s: #mu=%.2f, #sigma=%.2f", deltaLabel[j], hDeltaPt[j]->GetMean(), hDeltaPt[j]->GetRMS()), "l");
    }
    legendDelta->Draw(); gText.DrawLatex(0.6, 0.2, Form("PYTHIA, %s", row ? "R = 0.8" : "R = 0.4"));
    delete tTruth; delete tEmbMB; delete tEmbC;
  }
  c.Print(gOutDir + "/qa_jets.pdf"); c.Print(gOutDir + "/qa_jets.png");
}

// ---------------------------------------------------------------------------
// qa_spectra: in-cone spectra closure (PYTHIA, truth jets 20-40 GeV/c), per selected jet:
// truth n (no UE) vs embedded raw n, n - b (perpendicular cones) and n_cs
static void PlotInConeSpectra()
{
  TCanvas c("cs", "", 1500, 520); c.Divide(3, 1, 0.001, 0.001);
  struct Panel { int k; const char* embedding; const char* label; };
  Panel panels[3] = {{2, "embMB", "R = 0.4, UE O+O MB"}, {2, "embC", "R = 0.4, UE O+O 0-5%"}, {7, "embMB", "R = 0.8, UE O+O MB"}};
  for (int p = 0; p < 3; p++) {
    int k = panels[p].k;
    TChain* tTruth = OpenTreeChain("py_pp", "nom"); TChain* tEmb = OpenTreeChain("py_pp", panels[p].embedding);
    EventBranches eTruth, eEmb; BindBranches(tTruth, eTruth); BindBranches(tEmb, eEmb);
    TH1D* h[4]; const char* labels[4] = {"truth n (no UE)", "+UE: raw n", "+UE: n - b (perp)", "+UE: n_{cs}"};
    for (int j = 0; j < 4; j++) { h[j] = new TH1D(Form("hs%d%d", p, j), ";#it{p}_{T} (GeV/#it{c});#LTdN/d#it{p}_{T}#GT per jet", kNumPtBins, kPtEdges); StyleHist(h[j], j); }
    double sumWeights = 0; Long64_t n = std::min(tTruth->GetEntries(), tEmb->GetEntries());
    for (Long64_t i = 0; i < n; i++) {
      tTruth->GetEntry(i); tEmb->GetEntry(i);
      if (eTruth.jetPt[k] < 20 || eTruth.jetPt[k] >= 40) continue;   // truth-jet selection
      sumWeights += eTruth.weight;
      for (int b = 0; b < kNumPtBins; b++) {
        int ix = k * kNumPtBins + b;
        h[0]->Fill(h[0]->GetBinCenter(b + 1), eTruth.weight * eTruth.nIn[ix]); h[1]->Fill(h[1]->GetBinCenter(b + 1), eTruth.weight * eEmb.nIn[ix]);
        h[2]->Fill(h[2]->GetBinCenter(b + 1), eTruth.weight * (eEmb.nIn[ix] - eEmb.nPerp[ix])); h[3]->Fill(h[3]->GetBinCenter(b + 1), eTruth.weight * eEmb.nCS[ix]);
      }
    }
    c.cd(p + 1); SetPadMargins(); gPad->SetLogx(); gPad->SetLogy();
    TLegend* legend = MakeLegend(0.45, 0.66, 0.97, 0.93);
    for (int j = 0; j < 4; j++) {
      if (sumWeights > 0) h[j]->Scale(1. / sumWeights, "width");   // per jet, per GeV/c
      for (int b = 1; b <= kNumPtBins; b++) h[j]->SetBinError(b, 0);
      h[j]->GetYaxis()->SetRangeUser(1e-4, 50); h[j]->GetXaxis()->SetRangeUser(0.15, 100); h[j]->Draw(j ? "P SAME" : "P"); legend->AddEntry(h[j], labels[j], "p");
    }
    legend->Draw(); gText.SetTextSize(0.045); gText.DrawLatex(0.2, 0.2, Form("PYTHIA pp, %s, truth jet 20-40 GeV/#it{c}", panels[p].label));
    delete tTruth; delete tEmb;
  }
  c.Print(gOutDir + "/qa_spectra.pdf"); c.Print(gOutDir + "/qa_spectra.png");
}

// ---------------------------------------------------------------------------
// qa_reference: <[pT]_ref> (refbar_) and C2,ref (C2_) vs R, jets 40-60 GeV/c
static void PlotReferenceQA()
{
  TCanvas c("cr", "", 1200, 520); c.Divide(2, 1, 0.001, 0.001);
  struct Sample { const char* file; const char* tag; const char* label; };
  Sample samples[6] = {{"pp_py_closure", "py_pp_nom", "PYTHIA pp"}, {"pp_py_closure", "py_pp_embMB", "PYTHIA + UE(MB)"}, {"pp_py_closure", "py_pp_embC", "PYTHIA + UE(0-5%)"},
                       {"OO_main", "OO_vac_nom", "JEWEL O+O vac"}, {"OO_main", "OO_RC_cs10", "JEWEL O+O recoils"}, {"OO_embMB", "OO_RC_cs10_embMB", "O+O recoils + UE(MB)"}};
  const char* quantity[2] = {"refbar", "C2"};
  const char* yTitle[2] = {"#LT[#it{p}_{T}]_{ref}#GT (GeV/#it{c})", "#it{C}_{2,ref} (GeV^{2}/#it{c}^{2})"};
  for (int q = 0; q < 2; q++) {
    c.cd(q + 1); SetPadMargins(); if (q == 1) gPad->SetLogy();
    DrawFrame(0, 0.9, q ? 1e-4 : 0.6, q ? 0.1 : 1.0, "#it{R}", yTitle[q], false, q == 1);
    TLegend* legend = MakeLegend(0.45, 0.62, 0.97, 0.93);
    for (int s = 0; s < 6; s++) {
      TGraphErrors* g = new TGraphErrors();
      for (int iR = 0; iR < kNumRadii; iR++) {
        TH1D* h = LoadResult(kAcceptanceOfRadius[iR], samples[s].file, Form("%s_%s_c1_R%d_jet40_60", quantity[q], samples[s].tag, iR));
        if (!h) continue;
        int n = g->GetN(); g->SetPoint(n, kRadius[iR] + 0.008 * (s - 2.5), h->GetBinContent(1)); g->SetPointError(n, 0, h->GetBinError(1));
      }
      StyleGraph(g, kCurveColor[s], kCurveMarker[s]); g->Draw("PL SAME"); legend->AddEntry(g, samples[s].label, "lp");
    }
    legend->Draw(); gText.SetTextSize(0.045); gText.DrawLatex(0.18, 0.2, "jet 40-60 GeV/#it{c}, reference 0.2-3 GeV/#it{c}");
  }
  c.Print(gOutDir + "/qa_reference.pdf"); c.Print(gOutDir + "/qa_reference.png");
}

// ---------------------------------------------------------------------------
// qa_fmean: <f(pT)> (meanN<probe>_) at R = 0.4 and 0.8, jets 20-40 GeV/c
static void PlotMeanFraction()
{
  TCanvas c("cf", "", 1200, 520); c.Divide(2, 1, 0.001, 0.001);
  int radiusIndexOfPanel[2] = {2, 3};
  for (int p = 0; p < 2; p++) {
    c.cd(p + 1); SetPadMargins(); int iR = radiusIndexOfPanel[p];
    DrawFrame(0.2, 100, 1e-4, 1, "#it{p}_{T}^{probe} (GeV/#it{c})", "#LT#it{f}(#it{p}_{T})#GT", true, true);
    struct Curve { const char* file; const char* probe; const char* tag; const char* label; };
    Curve curves[5] = {{"pp_main", "", "py_pp_nom", "PYTHIA pp, raw f"}, {"pp_main", "_fsub", "py_pp_nom", "PYTHIA pp, f_{sub}"},
                       {"OO_embMB", "", "OO_RC_cs10_embMB", "O+O rec.+UE, raw f"}, {"OO_embMB", "_fsub", "OO_RC_cs10_embMB", "O+O rec.+UE, f_{sub}"},
                       {"OO_embMB", "_fcs", "OO_RC_cs10_embMB", "O+O rec.+UE, f_{cs}"}};
    TLegend* legend = MakeLegend(0.18, 0.18, 0.6, 0.45);
    for (int s = 0; s < 5; s++) {
      TH1D* h = LoadResult(kAcceptanceOfRadius[iR], curves[s].file, Form("meanN%s_%s_c1_R%d_jet20_40", curves[s].probe, curves[s].tag, iR));
      if (!h) continue;
      StyleHist(h, s); h->Draw("P SAME"); legend->AddEntry(h, curves[s].label, "p");
    }
    legend->Draw(); gText.SetTextSize(0.045); gText.DrawLatex(0.5, 0.86, Form("R = %s, jet 20-40 GeV/#it{c}", p ? "0.8" : "0.4"));
  }
  c.Print(gOutDir + "/qa_fmean.pdf"); c.Print(gOutDir + "/qa_fmean.png");
}

// ---------------------------------------------------------------------------
// qa_<name>: v0(pT) curves at R = 0.4 and 0.8 from one Phase2 file. Each curve is
// {probe, tag, tagMinus, label}. differenceMode: draw v0[tag] - v0[tagMinus] (jets
// 40-60); otherwise v0[tag] (jets 20-40, tagMinus unused).
static void PlotV0vsPt(const char* name, const char* file, std::vector<std::vector<const char*>> curves, const char* title,
                       bool differenceMode, double yMin, double yMax)
{
  TCanvas c(Form("c%s", name), "", 1200, 520); c.Divide(2, 1, 0.001, 0.001);
  int radiusIndexOfPanel[2] = {2, 3};
  for (int p = 0; p < 2; p++) {
    c.cd(p + 1); SetPadMargins(); int iR = radiusIndexOfPanel[p], acc = kAcceptanceOfRadius[iR];
    // x axis starts at 0.2 GeV/c so that every probe bin (first bin 0.2-0.5) lies inside the frame
    TH1D* frame = DrawFrame(0.2, 100, yMin, yMax, "#it{p}_{T}^{probe} (GeV/#it{c})", differenceMode ? "#Delta#it{v}_{0}(#it{p}_{T})" : "#it{v}_{0}(#it{p}_{T})^{jet-UE}", true);
    double yLow = 1e9, yHigh = -1e9;
    TLegend* legend = MakeLegend(0.18, 0.72, 0.75, 0.93);
    for (size_t s = 0; s < curves.size(); s++) {
      const char* probe = curves[s][0]; const char* tag = curves[s][1]; const char* tagMinus = curves[s][2]; const char* label = curves[s][3];
      TH1D* h = differenceMode ? Difference(LoadResult(acc, file, Form("v0_%s_%s_c1_R%d_jet40_60", probe, tag, iR)), LoadResult(acc, file, Form("v0_%s_%s_c1_R%d_jet40_60", probe, tagMinus, iR)))
                               : LoadResult(acc, file, Form("v0_%s_%s_c1_R%d_jet20_40", probe, tag, iR));
      if (!HideUnreliableBins(h)) continue;
      StyleHist(h, s); h->Draw("E1 X0 SAME"); legend->AddEntry(h, label, "lp"); growExtent(h, yLow, yHigh);
    }
    fitYRange(frame, yLow, yHigh, 0.32, 0.14);   // legend on top, label at the bottom
    legend->Draw(); gText.SetTextSize(0.045); gText.DrawLatex(0.2, 0.2, Form("%s, R = %s", title, p ? "0.8" : "0.4"));
  }
  c.Print(Form("%s/qa_%s.pdf", gOutDir.Data(), name)); c.Print(Form("%s/qa_%s.png", gOutDir.Data(), name));
}

// ---------------------------------------------------------------------------
// qa_refgap: range-averaged signal S(R) for the nominal and the eta-gap reference;
// qa_annulus: S(R_i) - S(R_{i-1}) for successive radii (nominal reference).
static void PlotSignalVsR()
{
  TCanvas c("cg", "", 1200, 520); c.Divide(2, 1, 0.001, 0.001);
  for (int p = 0; p < 2; p++) {   // panel 0: soft probes, 1: hard probes
    c.cd(p + 1); SetPadMargins(); DrawFrame(0, 0.9, -0.25, 0.2, "#it{R}", "#LT#Delta#it{v}_{0}#GT recoils - vacuum");
    TLegend* legend = MakeLegend(0.18, 0.75, 0.8, 0.93);
    const char* refs[2] = {"refnom", "refgap"}; const char* labels[2] = {"reference outside jet + perp. cones", "reference |#Delta#eta| > R + 0.5"};
    for (int r = 0; r < 2; r++) {
      TGraphErrors* g = new TGraphErrors();
      for (int iR = 0; iR < kNumRadii; iR++) {
        double mean, error;
        if (!RangeAverage(RecoilSignal(iR, refs[r]), kProbeRange[p][0], kProbeRange[p][1], mean, error)) continue;
        int n = g->GetN(); g->SetPoint(n, kRadius[iR] + 0.01 * (r - 0.5), mean); g->SetPointError(n, 0, error);
      }
      StyleGraph(g, kCurveColor[r], kCurveMarker[r]); g->Draw("PL SAME"); legend->AddEntry(g, labels[r], "lp");
    }
    legend->Draw(); gText.SetTextSize(0.045); gText.DrawLatex(0.2, 0.2, Form("O+O, UE O+O MB, jet 40-60 GeV/#it{c}, %s", kProbeRangeLabel[p]));
  }
  c.Print(gOutDir + "/qa_refgap.pdf"); c.Print(gOutDir + "/qa_refgap.png");

  // annulus: S(R2) - S(R1) for successive radii; errors in quadrature (same events: an overestimate)
  TCanvas a("ca", "", 1200, 520); a.Divide(2, 1, 0.001, 0.001);
  for (int p = 0; p < 2; p++) {
    a.cd(p + 1); SetPadMargins(); DrawFrame(0, 0.9, -0.2, 0.2, "outer radius #it{R}_{2}", "#LT#Delta#it{v}_{0}#GT(#it{R}_{2}) - #LT#Delta#it{v}_{0}#GT(#it{R}_{1})");
    TLegend* legend = MakeLegend(0.18, 0.78, 0.8, 0.93);
    double signal[kNumRadii], error[kNumRadii]; bool ok[kNumRadii];
    for (int iR = 0; iR < kNumRadii; iR++) ok[iR] = RangeAverage(RecoilSignal(iR), kProbeRange[p][0], kProbeRange[p][1], signal[iR], error[iR]);
    TGraphErrors* g = new TGraphErrors();
    for (int iR = 1; iR < kNumRadii; iR++)
      if (ok[iR] && ok[iR - 1]) { int n = g->GetN(); g->SetPoint(n, kRadius[iR], signal[iR] - signal[iR - 1]); g->SetPointError(n, 0, hypot(error[iR], error[iR - 1])); }
    StyleGraph(g, kRed + 1, 20); g->Draw("PL SAME");
    legend->AddEntry(g, "recoils - vacuum, successive radii (R_{1} = previous)", "lp"); legend->Draw();
    gText.SetTextSize(0.045); gText.DrawLatex(0.2, 0.2, Form("O+O, UE O+O MB, jet 40-60 GeV/#it{c}, %s", kProbeRangeLabel[p]));
  }
  a.Print(gOutDir + "/qa_annulus.pdf"); a.Print(gOutDir + "/qa_annulus.png");
}

// ---------------------------------------------------------------------------
void AppendixQA(const char* treeDir = "/dcache/alice/sbhatta/JewelV0pT/trees", const char* resultsDir = "/data/alice/sbhatta/JewelV0pT/results_coarse",
                const char* outDir = "qa", int nfiles = 20)
{
  gTreeDir = treeDir; gResultsDir = resultsDir; gOutDir = outDir; gNumFiles = nfiles; gSystem->mkdir(outDir, true);
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0); gStyle->SetTextFont(42); gStyle->SetLegendFont(42);
  gText.SetNDC(); gText.SetTextFont(42);
  // the tree-based checks need the per-event trees; skip them if the tree directory is gone
  if (gSystem->AccessPathName(treeDir)) printf("trees %s not found: qa_jets and qa_spectra not remade\n", treeDir);
  else { PlotJetQA(); PlotInConeSpectra(); }
  PlotReferenceQA(); PlotMeanFraction();
  // curves: {probe, tag, tagMinus (difference mode only), label}
  PlotV0vsPt("closure", "pp_py_closure", {{"fsub", "py_pp_nom", "", "PYTHIA, no UE (truth), f_{sub}"}, {"fsub", "py_pp_embMB", "", "+UE O+O MB, f_{sub}"},
                                         {"fsub", "py_pp_embC", "", "+UE O+O 0-5%, f_{sub}"}, {"fcs", "py_pp_embMB", "", "+UE O+O MB, f_{cs}"}},
             "PYTHIA pp, jet 20-40 GeV/#it{c}", false, -0.4, 1.2);
  PlotV0vsPt("recoilsub", "OO_recoilsub", {{"fsub", "OO_RC_nosub", "OO_RC_cs10", "no thermal subtraction - CS(#DeltaR<1)"}, {"fsub", "OO_RC_cs05", "OO_RC_cs10", "CS(#DeltaR<0.5) - CS(#DeltaR<1)"}},
             "O+O recoils, jet 40-60 GeV/#it{c}", true, -0.3, 0.3);
  PlotSignalVsR();
}
