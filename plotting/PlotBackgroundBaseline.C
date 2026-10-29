// ============================================================================
//  PlotBackgroundBaseline.C -- medium modification S and background-induced
//  change F of v0^{jet-UE}, versus jet radius, as soft- and hard-probe range
//  summaries (1-6 and 12-100 GeV/c, inverse-variance weighted over bins).
//
//  Definitions (paper notation):
//    S(U)   = v0[medium + U] - v0[vacuum + U]      (same embedded UE U)
//    F_g(U) = v0[g + U]      - v0[g]               (vacuum generator g)
//  with g = J (JEWEL vacuum, O+O) or P (PYTHIA 8 pp with MPI). Both terms of F
//  come from the same generated events, so its replica-wise error includes
//  their correlation.
//
//  Pages (minimum-bias configuration refnom_call_jt0):
//    1. thermal UE at MB density, and the same UE with its reference-window
//       <pT> matched to the jet-correlated reference particles (embMBmatch):
//       the composition term of the UE-only change should vanish at matching;
//    2. flow-correlated thermal UE (embFlow): multiplicity and <pT> follow the
//       centrality JEWEL used for the medium of the same event;
//    3. UE from Angantyr and HIJING events (embAng, embAngC, embHij).
//  A text table with all numbers is written next to the plots.
//
//  Run: root -l -b -q 'PlotBackgroundBaseline.C+("<resultsDir>","<outDir>")'
// ============================================================================
#include "JetUEResults.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include <fstream>
#include <string>
#include <vector>

using namespace jetue;

namespace {

struct Curve {
  std::string label;
  std::string plusTag, minusTag;   // difference plusTag - minusTag
  const char* color;
  int marker;
};

struct ProbeRange { const char* name; double low, high; };
const ProbeRange kSoft = {"soft probes, 1 < p_{T} < 6 GeV/#it{c}", 1.0, 6.0};
const ProbeRange kHard = {"hard probes, 12 < p_{T} < 100 GeV/#it{c}", 12.0, 100.0};

// S(U) for the recoil calculation; F for JEWEL vacuum and PYTHIA
Curve SignalRecoil(const std::string& emb, const std::string& label, const char* color, int marker) {
  return {label, "OO_RC_cs10_" + emb, "OO_vac_" + emb, color, marker};
}
Curve SignalNoRecoil(const std::string& emb, const std::string& label, const char* color, int marker) {
  return {label, "OO_NR_" + emb, "OO_vac_" + emb, color, marker};
}
Curve ChangeJewelVacuum(const std::string& emb, const std::string& label, const char* color, int marker) {
  return {label, "OO_vac_" + emb, "OO_vac_nom", color, marker};
}
Curve ChangePythia(const std::string& emb, const std::string& label, const char* color, int marker) {
  return {label, "py_pp_" + emb, "py_pp_nom", color, marker};
}

void DrawRange(const std::vector<ResultSet*>& resultsPerAcc, const std::vector<Curve>& curves,
               const ProbeRange& range, const std::string& jetClass, const std::string& probe,
               TLegend* legend, std::ofstream& table, const std::string& page) {
  TH1F* frame = new TH1F(Form("frame_%s_%s_%g", page.c_str(), jetClass.c_str(), range.low),
                         ";R (cone radius);#Deltav_{0}^{jet-UE}", 100, 0.0, 0.92);
  frame->SetDirectory(nullptr);
  // y range from the drawn points
  double yMin = 0, yMax = 0;
  std::vector<TGraphErrors*> graphs;
  for (size_t ic = 0; ic < curves.size(); ic++) {
    const Curve& cv = curves[ic];
    TGraphErrors* g = new TGraphErrors();
    for (int iR = 0; iR < kNumRadii; iR++) {
      const ResultSet* rs = resultsPerAcc[kNominalAcceptance[iR]];
      V0Curve d = Difference(rs->Get(cv.plusTag, iR, jetClass, probe), rs->Get(cv.minusTag, iR, jetClass, probe),
                             Form("d_%s_%zu_%d", page.c_str(), ic, iR));
      Summary s = RangeSummary(d, range.low, range.high);
      if (!s.ok()) continue;
      double x = kRadius[iR] + (ic - 0.5 * (curves.size() - 1)) * 0.012;
      int n = g->GetN();
      g->SetPoint(n, x, s.value); g->SetPointError(n, 0, s.error);
      yMin = std::min(yMin, s.value - s.error); yMax = std::max(yMax, s.value + s.error);
      table << Form("%-10s %-9s %-26s %-44s R=%s  %+.4f +- %.4f\n", page.c_str(), jetClass.c_str(),
                    range.low < 5 ? "soft(1-6)" : "hard(12-100)", cv.label.c_str(), kRadiusLabel[iR], s.value, s.error);
    }
    int color = TColor::GetColor(cv.color);
    g->SetMarkerStyle(cv.marker); g->SetMarkerColor(color); g->SetLineColor(color); g->SetMarkerSize(1.1);
    graphs.push_back(g);
    if (legend) legend->AddEntry(g, cv.label.c_str(), "p");
  }
  double pad = 0.15 * (yMax - yMin + 0.02);
  frame->SetMinimum(yMin - pad); frame->SetMaximum(yMax + 2.5 * pad);
  frame->Draw("AXIS");
  TLine* zero = new TLine(0.0, 0.0, 0.92, 0.0); zero->SetLineStyle(2); zero->Draw();
  for (TGraphErrors* g : graphs) if (g->GetN() > 0) g->Draw("P same");
  TLatex tl; tl.SetNDC(); tl.SetTextSize(0.045);
  tl.DrawLatex(0.15, 0.88, range.name);
}

}  // namespace

void PlotBackgroundBaseline(const char* resultsDir = "/data/alice/sbhatta/JewelV0pT/results_v3",
                            const char* outDir = "plots_v3", const char* probe = "_fsub") {
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, true);
  const std::string cfg = "refnom_call_jt0";
  std::vector<ResultSet*> resultsPerAcc = {new ResultSet(resultsDir, cfg, 0), new ResultSet(resultsDir, cfg, 1)};

  // Okabe-Ito colours; filled markers = JEWEL, open = PYTHIA, stars = signal
  const std::vector<std::pair<std::string, std::vector<Curve>>> pages = {
    {"thermal", {SignalRecoil("embMB", "S: recoils - vacuum, thermal MB", "#D55E00", 29),
                 ChangeJewelVacuum("embMB", "F_{J}: JEWEL vac., thermal MB", "#0072B2", 20),
                 ChangePythia("embMB", "F_{P}: PYTHIA, thermal MB", "#0072B2", 24),
                 ChangeJewelVacuum("embMBmatch", "F_{J}: JEWEL vac., matched #LTp_{T}#GT", "#009E73", 21),
                 ChangePythia("embMBmatch", "F_{P}: PYTHIA, matched #LTp_{T}#GT", "#009E73", 25)}},
    {"flow",    {SignalRecoil("embFlow", "S: recoils - vacuum, flow-correlated UE", "#D55E00", 29),
                 SignalNoRecoil("embFlow", "no recoils - vacuum, flow-correlated UE", "#E69F00", 33),
                 ChangeJewelVacuum("embFlow", "F_{J}: JEWEL vac., flow-correlated UE", "#0072B2", 20),
                 ChangePythia("embFlow", "F_{P}: PYTHIA, flow-correlated UE", "#0072B2", 24)}},
    {"generator", {ChangeJewelVacuum("embAng", "F_{J}: Angantyr UE (dN/d#eta 30-60)", "#CC79A7", 20),
                   ChangePythia("embAng", "F_{P}: Angantyr UE (dN/d#eta 30-60)", "#CC79A7", 24),
                   ChangeJewelVacuum("embHij", "F_{J}: HIJING UE (dN/d#eta 30-60)", "#56B4E9", 21),
                   ChangePythia("embHij", "F_{P}: HIJING UE (dN/d#eta 30-60)", "#56B4E9", 25),
                   ChangeJewelVacuum("embMB", "F_{J}: thermal MB (for comparison)", "#000000", 22),
                   ChangePythia("embMB", "F_{P}: thermal MB (for comparison)", "#000000", 26)}},
  };
  const std::vector<std::string> jetClasses = {"jet20_40", "jet40_60"};

  std::ofstream table(Form("%s/background_baseline_%s.txt", outDir, probe[0] ? probe + 1 : "f"));
  table << "# page  jetClass  range  difference  R  value +- replica-wise bootstrap error\n";
  for (const auto& page : pages) {
    for (const std::string& jc : jetClasses) {
      TCanvas* c = new TCanvas(Form("c_%s_%s", page.first.c_str(), jc.c_str()), "", 1100, 470);
      c->Divide(3, 1, 0.002, 0.002);
      TLegend* legend = new TLegend(0.02, 0.15, 0.98, 0.85);
      legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.055);
      c->cd(1); gPad->SetLeftMargin(0.15);
      DrawRange(resultsPerAcc, page.second, kSoft, jc, probe, legend, table, page.first);
      c->cd(2); gPad->SetLeftMargin(0.15);
      DrawRange(resultsPerAcc, page.second, kHard, jc, probe, nullptr, table, page.first);
      c->cd(3);
      legend->SetHeader(Form("O+O MB, %s, probe %s", jc == "jet20_40" ? "20-40 GeV/#it{c} jets" : "40-60 GeV/#it{c} jets",
                             probe[0] ? probe + 1 : "f"));
      legend->Draw();
      TString base = Form("%s/background_baseline_%s_%s", outDir, page.first.c_str(), jc.c_str());
      c->Print(base + ".pdf"); c->Print(base + ".png");
    }
  }
  printf("table -> %s/background_baseline_%s.txt\n", outDir, probe[0] ? probe + 1 : "f");
}
