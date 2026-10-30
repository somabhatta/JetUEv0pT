// ============================================================================
//  PlotDeltaV0vsR.C -- medium modification of v0^{jet-UE}(pT) versus jet radius,
//  with the probe-pT bins drawn side by side at each radius.
//
//  For one analysis configuration and one jet class, each page shows
//    columns: no recoils - vacuum,  recoils - vacuum
//    rows   : the embedded-UE scenarios of that configuration
//  where every difference is taken at the SAME embedded UE (S(b;U) of the
//  paper): e.g. "thermal MB" row = OO_NR_embMB - OO_vac_embMB.
//  R <= 0.4 uses |eta_trk| < 0.9, R = 0.8 uses |eta_trk| < 2.5.
//  Errors: replica-wise bootstrap (see JetUEResults.h). Every probe bin with an
//  estimate is drawn; the vertical range of a page covers all its points.
//
//  Run (from this directory):
//    root -l -b -q 'PlotDeltaV0vsR.C+("<resultsDir>","<outDir>")'
//  with <resultsDir> = .../results_v3 (contains <cfg>/acc<i>/v0pt_*.root).
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
#include <string>
#include <vector>

using namespace jetue;

namespace {

struct UEScenario {
  std::string label;       // row title
  std::string embedding;   // "" (no added UE) or the variant suffix, e.g. "embMB"
};

struct Configuration {
  std::string cfg;         // results subdirectory
  std::string title;       // event-class title on the page
  std::vector<UEScenario> scenarios;
};

// sample tags at a given embedding: vacuum, medium without recoils, with recoils
std::string VacuumTag(const std::string& emb)    { return emb.empty() ? "OO_vac_nom" : "OO_vac_" + emb; }
std::string NoRecoilTag(const std::string& emb)  { return emb.empty() ? "OO_NR_nom"  : "OO_NR_" + emb; }
std::string RecoilTag(const std::string& emb)    { return emb.empty() ? "OO_RC_cs10" : "OO_RC_cs10_" + emb; }

// colour-blind-safe palette (Okabe-Ito) and markers, one per probe-pT bin
const char* kBinColor[8] = {"#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#0072B2", "#999999"};
const int   kBinMarker[8] = {20, 21, 22, 23, 33, 29, 34, 47};

void DrawPanel(const std::vector<ResultSet*>& resultsPerAcc, const std::string& minusTag,
               const std::string& plusTag, const std::string& jetClass, const std::string& probe,
               const std::string& title, double yMin, double yMax, TLegend* legend) {
  TH1F* frame = new TH1F(Form("frame_%s_%s_%s", plusTag.c_str(), minusTag.c_str(), jetClass.c_str()),
                         Form(";R (cone radius);#Deltav_{0}^{jet-UE}"), 100, 0.0, 0.92);
  frame->SetMinimum(yMin); frame->SetMaximum(yMax);
  frame->SetDirectory(nullptr);
  // text sizes chosen so that the page stays legible when printed at page width
  frame->GetXaxis()->SetLabelSize(0.065); frame->GetYaxis()->SetLabelSize(0.065);
  frame->GetXaxis()->SetTitleSize(0.07); frame->GetYaxis()->SetTitleSize(0.07);
  frame->GetXaxis()->SetTitleOffset(0.9); frame->GetYaxis()->SetTitleOffset(0.95);
  frame->GetYaxis()->SetNdivisions(505);
  frame->Draw("AXIS");
  TLine* zero = new TLine(0.0, 0.0, 0.92, 0.0); zero->SetLineStyle(2); zero->Draw();

  const double spacing = 0.012;   // horizontal offset between probe-pT bins
  for (int iR = 0; iR < kNumRadii; iR++) {
    const ResultSet* rs = resultsPerAcc[kNominalAcceptance[iR]];
    V0Curve plus = rs->Get(plusTag, iR, jetClass, probe);
    V0Curve minus = rs->Get(minusTag, iR, jetClass, probe);
    V0Curve diff = Difference(plus, minus, Form("d_%s_%d", plusTag.c_str(), iR));
    if (!diff.ok()) continue;
    const int nBins = diff.central->GetNbinsX();
    for (int q = 1; q <= nBins; q++) {
      double err = diff.central->GetBinError(q);
      if (err <= 0) continue;   // bin without an estimate (no events)
      double x = kRadius[iR] + (q - 0.5 * (nBins + 1)) * spacing;
      TGraphErrors* g = new TGraphErrors(1);
      g->SetPoint(0, x, diff.central->GetBinContent(q));
      g->SetPointError(0, 0, err);
      int color = TColor::GetColor(kBinColor[(q - 1) % 8]);
      g->SetMarkerStyle(kBinMarker[(q - 1) % 8]); g->SetMarkerColor(color); g->SetLineColor(color);
      g->SetMarkerSize(1.2); g->SetLineWidth(2);
      g->Draw("P same");
      if (legend && iR == kNumRadii - 1)
        legend->AddEntry(g, Form("%g-%g GeV/#it{c}", diff.central->GetBinLowEdge(q),
                                  diff.central->GetBinLowEdge(q) + diff.central->GetBinWidth(q)), "p");
    }
  }
  TLatex tl; tl.SetNDC(); tl.SetTextSize(0.068);
  tl.DrawLatex(0.16, 0.86, title.c_str());
}

// smallest and largest value -/+ error of a difference over all radii and bins
void GrowExtent(const std::vector<ResultSet*>& resultsPerAcc, const std::string& minusTag, const std::string& plusTag,
                const std::string& jetClass, const std::string& probe, double& lo, double& hi) {
  for (int iR = 0; iR < kNumRadii; iR++) {
    const ResultSet* rs = resultsPerAcc[kNominalAcceptance[iR]];
    V0Curve d = Difference(rs->Get(plusTag, iR, jetClass, probe), rs->Get(minusTag, iR, jetClass, probe), Form("ext_%s_%d", plusTag.c_str(), iR));
    if (!d.ok()) continue;
    for (int q = 1; q <= d.central->GetNbinsX(); q++) {
      double v = d.central->GetBinContent(q), e = d.central->GetBinError(q);
      if (e <= 0) continue;
      lo = std::min(lo, v - e); hi = std::max(hi, v + e);
    }
  }
}

}  // namespace

void PlotDeltaV0vsR(const char* resultsDir = "/data/alice/sbhatta/JewelV0pT/results_v3",
                    const char* outDir = "plots_v3", const char* probe = "_fsub",
                    double unusedMaxError = 0) {   // kept for call compatibility; no points are hidden
  (void)unusedMaxError;
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, true);

  const std::vector<Configuration> configurations = {
    {"refnom_call_jt0", "O+O 0-100%",
     {{"no added UE", ""}, {"thermal UE, MB density", "embMB"}, {"flow-correlated UE", "embFlow"}}},
    {"refnom_c0-40_jt0", "O+O 0-40%",
     {{"no added UE", ""}, {"thermal UE, N_{part}-scaled", "embNp"}, {"flow-correlated UE", "embFlow"}}},
    {"refnom_c40-80_jt0", "O+O 40-80%",
     {{"no added UE", ""}, {"thermal UE, N_{part}-scaled", "embNp"}, {"flow-correlated UE", "embFlow"}}},
  };
  const std::vector<std::pair<std::string, std::string>> jetClasses = {
    {"jet20_40", "20 < p_{T}^{jet,corr} < 40 GeV/#it{c}"}, {"jet40_60", "40 < p_{T}^{jet,corr} < 60 GeV/#it{c}"}};

  for (const Configuration& conf : configurations) {
    std::vector<ResultSet*> resultsPerAcc = {new ResultSet(resultsDir, conf.cfg, 0), new ResultSet(resultsDir, conf.cfg, 1)};
    if (resultsPerAcc[0]->empty() && resultsPerAcc[1]->empty()) { printf("no results for %s\n", conf.cfg.c_str()); continue; }
    for (const auto& jc : jetClasses) {
      const int nRows = conf.scenarios.size();
      TCanvas* c = new TCanvas(Form("c_%s_%s", conf.cfg.c_str(), jc.first.c_str()), "", 1000, 330 * nRows + 60);
      // header strip on top, panel grid below
      // header: title line and the probe-pT key
      const double headerFrac = 110.0 / (330 * nRows + 60);
      TPad* header = new TPad("header", "", 0.0, 1.0 - headerFrac, 1.0, 1.0);
      TPad* body = new TPad("body", "", 0.0, 0.0, 1.0, 1.0 - headerFrac);
      header->Draw(); body->Draw();
      body->Divide(2, nRows, 0.002, 0.002);
      TLegend* legend = new TLegend(0.01, 0.02, 0.99, 0.58);
      legend->SetNColumns(4); legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.20);
      // common vertical range for the page, covering every point
      double yMin = 0, yMax = 0;
      for (const UEScenario& sc : conf.scenarios) {
        GrowExtent(resultsPerAcc, VacuumTag(sc.embedding), NoRecoilTag(sc.embedding), jc.first, probe, yMin, yMax);
        GrowExtent(resultsPerAcc, VacuumTag(sc.embedding), RecoilTag(sc.embedding), jc.first, probe, yMin, yMax);
      }
      { double pad = 0.05 * (yMax - yMin); yMin -= pad; yMax += 0.25 * (yMax - yMin); }   // room for the panel title
      for (int row = 0; row < nRows; row++) {
        const UEScenario& sc = conf.scenarios[row];
        body->cd(2 * row + 1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.02);
        DrawPanel(resultsPerAcc, VacuumTag(sc.embedding), NoRecoilTag(sc.embedding), jc.first, probe,
                  "no recoils - vacuum, " + sc.label, yMin, yMax, nullptr);
        body->cd(2 * row + 2); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.02);
        DrawPanel(resultsPerAcc, VacuumTag(sc.embedding), RecoilTag(sc.embedding), jc.first, probe,
                  "recoils - vacuum, " + sc.label, yMin, yMax, row == 0 ? legend : nullptr);
      }
      header->cd();
      TLatex head; head.SetNDC(); head.SetTextSize(0.24);
      head.DrawLatex(0.01, 0.68, Form("JEWEL %s, #sqrt{s_{NN}} = 5.36 TeV, %s, probe %s",
                                       conf.title.c_str(), jc.second.c_str(), probe[0] ? probe + 1 : "f"));
      legend->Draw();   // probe-pT key, filled while drawing the first recoil panel
      TString base = Form("%s/deltaV0_vsR_%s_%s", outDir, conf.cfg.c_str(), jc.first.c_str());
      c->cd(); c->Print(base + ".pdf"); c->Print(base + ".png");
      // probe-pT legend on its own canvas
      TCanvas* cl = new TCanvas(Form("cl_%s_%s", conf.cfg.c_str(), jc.first.c_str()), "", 900, 90);
      legend->Draw();
      cl->Print(base + "_legend.pdf"); cl->Print(base + "_legend.png");
    }
  }
}
