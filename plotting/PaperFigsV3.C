// ============================================================================
//  PaperFigsV3.C -- figures and the summary table of the paper from the
//  2026-10-28 production (Phase2 results_v3/, Phase1 phase1_v3/).
//
//  All curves use the perpendicular-cone-subtracted probe f_sub, the nominal
//  reference (0.2 < pT < 3 GeV/c outside the jet and perpendicular cones), and
//  R <= 0.4 with |eta_trk| < 0.9, R = 0.8 with |eta_trk| < 2.5. Range summaries:
//  soft 1-6 and hard 12-100 GeV/c, inverse-variance weighted over probe bins,
//  errors from the bootstrap replicas (JetUEResults.h).
//
//  Notation (paper Eqs. for S and F): at embedded background U,
//    S(U)   = v0[medium + U] - v0[vacuum + U]   (medium = recoils or no recoils)
//    F_g(U) = v0[g + U] - v0[g],  g = J (JEWEL O+O vacuum) or P (PYTHIA 8 pp)
//  Background scenarios U: none; independent thermal pions at the O+O MB density
//  (embMB); thermal pions with multiplicity and <pT> that follow the event
//  centrality (embFlow); N_part-scaled thermal pions (embNp); thermal pions with
//  matched reference <pT> (embMBmatch); Angantyr / HIJING events (embAng, embHij).
//
//  Outputs (in <outDir>):
//    figV3_signal_vsR.pdf     S vs R, soft and hard, minimum bias: recoils for
//                             no / thermal / flow-correlated background, and the
//                             no-recoil control with the flow-correlated background
//    figV3_centrality.pdf     S vs R for O+O 0-40% and 40-80%, recoils, for no /
//                             N_part-scaled / flow-correlated background
//    figV3_baseline.pdf       F_J and F_P vs R for the background scenarios
//    figV3_c2_vs_nch.pdf      per-multiplicity-bin sqrt(C2,ref) and v0 (R = 0.4,
//                             hard probes) for thermal vs flow-correlated background
//    tableV3_signals.tex      S, its error and Z per scenario, R = 0.4 and 0.8
//
//  Run (stbc, LCG_107 view):
//    root -l -b -q 'PaperFigsV3.C+("/data/alice/sbhatta/JewelV0pT/results_v3",
//                   "/dcache/alice/sbhatta/JewelV0pT/phase1_v3","figs_v3")'
// ============================================================================
#include "JetUEResults.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include <fstream>
#include <string>
#include <vector>

using namespace jetue;

namespace {

const char* kProbe = "_fsub";
const std::string kJetClass = "jet40_60";
struct ProbeRange { const char* label; double low, high; };
const ProbeRange kSoft = {"1 < #it{p}_{T} < 6 GeV/#it{c}", 1.0, 6.0};
const ProbeRange kHard = {"12 < #it{p}_{T} < 100 GeV/#it{c}", 12.0, 100.0};
const double kMaxError = 0.15;   // range summaries with a larger error are not drawn

// one difference "plusTag - minusTag" drawn as a curve vs R
struct Curve { std::string label, plusTag, minusTag, color; int marker; };

std::string Tag(const std::string& sample, const std::string& embedding) {
  if (sample == "OO_RC") return embedding.empty() ? "OO_RC_cs10" : "OO_RC_cs10_" + embedding;
  return embedding.empty() ? sample + "_nom" : sample + "_" + embedding;
}
Curve Signal(const std::string& medium, const std::string& emb, const std::string& label, const char* color, int marker) {
  return {label, Tag(medium, emb), Tag("OO_vac", emb), color, marker};
}
Curve Baseline(const std::string& vacuum, const std::string& emb, const std::string& label, const char* color, int marker) {
  return {label, Tag(vacuum, emb), Tag(vacuum, ""), color, marker};
}

void SetPadMargins() {
  gPad->SetTicks(1, 1); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.04);
}

TH1D* DrawFrame(double x0, double x1, double y0, double y1, const char* xTitle, const char* yTitle) {
  static int id = 0;
  TH1D* frame = new TH1D(Form("frameV3_%d", id++), "", 10, x0, x1);
  frame->SetDirectory(nullptr);
  frame->GetYaxis()->SetRangeUser(y0, y1);
  frame->GetXaxis()->SetTitle(xTitle); frame->GetYaxis()->SetTitle(yTitle);
  frame->GetXaxis()->SetTitleSize(0.055); frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetXaxis()->SetLabelSize(0.05); frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame->Draw("AXIS");
  TLine* zero = new TLine(x0, 0, x1, 0); zero->SetLineStyle(2); zero->SetLineColor(kGray + 1); zero->Draw();
  return frame;
}

// summary of one curve at every radius; returns the graph (offset in x by slot)
TGraphErrors* SummaryVsR(const std::vector<ResultSet*>& resultsPerAcc, const Curve& c, const ProbeRange& range,
                         int slot, int nSlots, double& yLow, double& yHigh) {
  TGraphErrors* g = new TGraphErrors();
  for (int iR = 0; iR < kNumRadii; iR++) {
    const ResultSet* rs = resultsPerAcc[kNominalAcceptance[iR]];
    V0Curve d = Difference(rs->Get(c.plusTag, iR, kJetClass, kProbe), rs->Get(c.minusTag, iR, kJetClass, kProbe),
                           Form("dV3_%s_%s_%d", c.plusTag.c_str(), c.minusTag.c_str(), iR));
    Summary s = RangeSummary(d, range.low, range.high);
    if (!s.ok() || s.error > kMaxError) continue;   // hide points with large statistical error
    double x = kRadius[iR] + (slot - 0.5 * (nSlots - 1)) * 0.018;
    int n = g->GetN(); g->SetPoint(n, x, s.value); g->SetPointError(n, 0, s.error);
    yLow = std::min(yLow, s.value - s.error); yHigh = std::max(yHigh, s.value + s.error);
  }
  int color = TColor::GetColor(c.color.c_str());
  g->SetMarkerStyle(c.marker); g->SetMarkerColor(color); g->SetLineColor(color); g->SetLineWidth(2); g->SetMarkerSize(1.2);
  return g;
}

// two panels (soft, hard) of curves vs R; legend in the left panel
void TwoRangeFigure(const std::vector<ResultSet*>& resultsPerAcc, const std::vector<Curve>& curves,
                    const std::string& header, const char* yTitle, const std::string& fileBase) {
  static int canvasId = 0;
  TCanvas c(Form("cTwoRange%d", canvasId++), "", 1200, 560); c.Divide(2, 1, 0.001, 0.001);
  const ProbeRange* ranges[2] = {&kSoft, &kHard};
  for (int p = 0; p < 2; p++) {
    c.cd(p + 1); SetPadMargins();
    double yLow = 0, yHigh = 0;
    std::vector<TGraphErrors*> graphs;
    for (size_t i = 0; i < curves.size(); i++)
      graphs.push_back(SummaryVsR(resultsPerAcc, curves[i], *ranges[p], i, curves.size(), yLow, yHigh));
    // data fill the lower part; top 40% for text and legend
    double span = std::max(yHigh - yLow, 0.05), total = span / 0.52;
    TH1D* frame = DrawFrame(0, 0.9, yLow - 0.08 * total, yHigh + 0.40 * total, "#it{R}", yTitle);
    (void)frame;
    TLegend* legend = new TLegend(0.19, 0.89 - 0.055 * curves.size() / (curves.size() > 4 ? 2 : 1), 0.97, 0.89);
    if (curves.size() > 4) legend->SetNColumns(2);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.036);
    for (size_t i = 0; i < graphs.size(); i++) { graphs[i]->Draw("P SAME"); legend->AddEntry(graphs[i], curves[i].label.c_str(), "p"); }
    if (p == 0) legend->Draw();
    TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(0.042);
    text.DrawLatex(0.19, 0.91, (header + ", " + ranges[p]->label).c_str());
  }
  c.Print((fileBase + ".pdf").c_str()); c.Print((fileBase + ".png").c_str());
}

// sqrt(C2,ref) and v0 per fine multiplicity bin from the Phase1 moments
// (f_sub probe, R = 0.4, jets 40-60, probe pT range): v0_j = cov_j / (<f>_j sqrt(C2_j))
void PerNchCurves(const std::string& phase1Dir, const std::string& tag, double ptLow, double ptHigh,
                  TGraph* sigma, TGraph* v0) {
  TFile fm((phase1Dir + "/refnom_call_jt0/acc0/" + tag + "_mean.root").c_str());
  TFile fc((phase1Dir + "/refnom_call_jt0/acc0/" + tag + "_corr.root").c_str());
  auto fraction = (TProfile2D*)fm.Get("pN_NchPt_fsub_R2_jet40_60");
  auto refMean = (TProfile*)fm.Get("pRefFsub_Nch_R2_jet40_60");
  auto fractionTimesRef = (TProfile2D*)fc.Get("pNR_NchPt_fsub_R2_jet40_60");
  auto c2 = (TProfile*)fc.Get("pC2Sub_Nch_R2_jet40_60");
  if (!fraction || !refMean || !fractionTimesRef || !c2) return;
  // merge 5 fine bins (width 25 in N_ch) for display, event-weighted
  const int group = 5;
  for (int j0 = 1; j0 + group - 1 <= refMean->GetNbinsX(); j0 += group) {
    double w = 0, wSigma = 0, wV0 = 0;
    for (int j = j0; j < j0 + group; j++) {
      double n = refMean->GetBinEntries(j); if (n <= 0) continue;
      double mu = refMean->GetBinContent(j), s = std::sqrt(std::max(c2->GetBinContent(j), 0.0)), cov = 0, f = 0;
      for (int b = 1; b <= fraction->GetNbinsY(); b++) {
        double lo = fraction->GetYaxis()->GetBinLowEdge(b);
        if (lo < ptLow - 1e-6 || lo >= ptHigh - 1e-6) continue;
        cov += fractionTimesRef->GetBinContent(j, b) - fraction->GetBinContent(j, b) * mu;
        f += fraction->GetBinContent(j, b);
      }
      if (s <= 0 || f == 0) continue;
      w += n; wSigma += n * s; wV0 += n * cov / (f * s);
    }
    if (w <= 0) continue;
    double x = 0.5 * (refMean->GetBinLowEdge(j0) + refMean->GetBinLowEdge(j0 + group));
    sigma->SetPoint(sigma->GetN(), x, wSigma / w);
    v0->SetPoint(v0->GetN(), x, wV0 / w);
  }
}

}  // namespace

void PaperFigsV3(const char* resultsDir = "/data/alice/sbhatta/JewelV0pT/results_v3",
                 const char* phase1Dir = "/dcache/alice/sbhatta/JewelV0pT/phase1_v3",
                 const char* outDir = "figs_v3") {
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0);
  gStyle->SetTextFont(42); gStyle->SetLegendFont(42);
  gSystem->mkdir(outDir, true);
  auto results = [&](const char* cfg) {
    return std::vector<ResultSet*>{new ResultSet(resultsDir, cfg, 0), new ResultSet(resultsDir, cfg, 1)};
  };
  std::vector<ResultSet*> mb = results("refnom_call_jt0");

  // ---- S vs R, minimum bias ----
  TwoRangeFigure(mb,
    {Signal("OO_RC", "", "recoils, no added background", "#000000", 20),
     Signal("OO_RC", "embMB", "recoils, thermal background", "#D55E00", 21),
     Signal("OO_RC", "embFlow", "recoils, flow-correlated background", "#0072B2", 33),
     Signal("OO_NR", "embFlow", "no recoils, flow-correlated background", "#56B4E9", 24)},
    "O+O 0-100%", "#it{S} = #LT#Delta#it{v}_{0}#GT", std::string(outDir) + "/figV3_signal_vsR");

  // ---- S vs R in centrality classes ----
  for (const char* cls : {"0-40", "40-80"}) {
    std::vector<ResultSet*> rs = results(Form("refnom_c%s_jt0", cls));
    TwoRangeFigure(rs,
      {Signal("OO_RC", "", "recoils, no added background", "#000000", 20),
       Signal("OO_RC", "embNp", "recoils, N_{part}-scaled background", "#009E73", 22),
       Signal("OO_RC", "embFlow", "recoils, flow-correlated background", "#0072B2", 33),
       Signal("OO_NR", "embFlow", "no recoils, flow-correlated background", "#56B4E9", 24)},
      std::string("O+O ") + cls + "%", "#it{S} = #LT#Delta#it{v}_{0}#GT", std::string(outDir) + "/figV3_centrality_" + cls);
  }

  // ---- background-induced change F ----
  TwoRangeFigure(mb,
    {Baseline("OO_vac", "embMB", "#it{F}_{J}, thermal", "#D55E00", 20),
     Baseline("py_pp", "embMB", "#it{F}_{P}, thermal", "#D55E00", 24),
     Baseline("OO_vac", "embMBmatch", "#it{F}_{J}, matched #LT#it{p}_{T}#GT", "#009E73", 21),
     Baseline("py_pp", "embMBmatch", "#it{F}_{P}, matched #LT#it{p}_{T}#GT", "#009E73", 25),
     Baseline("OO_vac", "embFlow", "#it{F}_{J}, flow-correlated", "#0072B2", 33),
     Baseline("py_pp", "embFlow", "#it{F}_{P}, flow-correlated", "#0072B2", 27),
     Baseline("OO_vac", "embAng", "#it{F}_{J}, Angantyr", "#CC79A7", 22),
     Baseline("py_pp", "embAng", "#it{F}_{P}, Angantyr", "#CC79A7", 26)},
    "background-induced change", "#it{F} = #LT#Delta#it{v}_{0}#GT", std::string(outDir) + "/figV3_baseline");

  // ---- per-multiplicity sqrt(C2) and v0, hard probes, R = 0.4 ----
  {
    TCanvas c("cNch", "", 1200, 560); c.Divide(2, 1, 0.001, 0.001);
    struct S { const char* tag; const char* label; const char* color; int marker; };
    std::vector<S> samples = {{"OO_vac_embMB", "vacuum, thermal", "#D55E00", 24}, {"OO_RC_cs10_embMB", "recoils, thermal", "#D55E00", 20},
                              {"OO_vac_embFlow", "vacuum, flow-correlated", "#0072B2", 25}, {"OO_RC_cs10_embFlow", "recoils, flow-correlated", "#0072B2", 21}};
    std::vector<TGraph*> sig, v0;
    for (auto& s : samples) {
      TGraph* gs = new TGraph(); TGraph* gv = new TGraph();
      PerNchCurves(phase1Dir, s.tag, kHard.low, kHard.high, gs, gv);
      for (TGraph* g : {gs, gv}) { int col = TColor::GetColor(s.color); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(s.marker); g->SetMarkerSize(1.1); }
      sig.push_back(gs); v0.push_back(gv);
    }
    c.cd(1); SetPadMargins(); gPad->SetLogy();
    TH1D* f1 = DrawFrame(0, 400, 0.008, 0.3, "#it{N}_{ch} (|#it{#eta}| < 0.8, 0.2 < #it{p}_{T} < 5 GeV/#it{c})", "#sqrt{#it{C}_{2,ref}} (GeV/#it{c})");
    (void)f1;
    TLegend* lg = new TLegend(0.45, 0.68, 0.97, 0.94); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04);
    for (size_t i = 0; i < samples.size(); i++) { sig[i]->Draw("P SAME"); lg->AddEntry(sig[i], samples[i].label, "p"); }
    lg->Draw();
    c.cd(2); SetPadMargins();
    DrawFrame(0, 400, -0.6, 0.8, "#it{N}_{ch} (|#it{#eta}| < 0.8, 0.2 < #it{p}_{T} < 5 GeV/#it{c})", "#it{v}_{0} per #it{N}_{ch} bin, 12 < #it{p}_{T} < 100 GeV/#it{c}");
    for (TGraph* g : v0) g->Draw("P SAME");
    TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(0.042);
    text.DrawLatex(0.19, 0.91, "O+O 0-100%, #it{R} = 0.4, 40 < #it{p}_{T}^{jet,corr} < 60 GeV/#it{c}");
    c.Print(Form("%s/figV3_c2_vs_nch.pdf", outDir)); c.Print(Form("%s/figV3_c2_vs_nch.png", outDir));
  }

  // ---- table: S, sigma_S and Z per scenario (recoils), R = 0.4 and 0.8 ----
  // Z = S / sqrt(sigma_S^2 + sigma_Fbar^2 + delta_F^2), Fbar = (F_J + F_P)/2,
  // sigma_Fbar^2 = (sigma_FJ^2 + sigma_FP^2)/4, delta_F = |F_J - F_P|/2 (same background U)
  {
    std::ofstream tex(Form("%s/tableV3_signals.tex", outDir));
    struct Row { const char* cfg; const char* label; const char* emb; };
    std::vector<Row> rows = {{"refnom_call_jt0", "0--100\\%, thermal", "embMB"}, {"refnom_call_jt0", "0--100\\%, flow-correlated", "embFlow"},
                             {"refnom_c0-40_jt0", "0--40\\%, $N_{\\mathrm{part}}$-scaled", "embNp"}, {"refnom_c0-40_jt0", "0--40\\%, flow-correlated", "embFlow"},
                             {"refnom_c40-80_jt0", "40--80\\%, $N_{\\mathrm{part}}$-scaled", "embNp"}, {"refnom_c40-80_jt0", "40--80\\%, flow-correlated", "embFlow"}};
    tex << "% generated by PaperFigsV3.C: S (recoils - vacuum, same background), statistical error, Z\n";
    for (const Row& r : rows) {
      std::vector<ResultSet*> rs = results(r.cfg);
      // F is taken from the minimum-bias configuration (vacuum samples have no centrality)
      tex << r.label;
      for (int iR : {2, 3}) {
        const ResultSet* set = rs[kNominalAcceptance[iR]];
        const ResultSet* setMB = mb[kNominalAcceptance[iR]];
        for (const ProbeRange* range : {&kSoft, &kHard}) {
          Summary s = RangeSummary(Difference(set->Get(Tag("OO_RC", r.emb), iR, kJetClass, kProbe), set->Get(Tag("OO_vac", r.emb), iR, kJetClass, kProbe), "tS"), range->low, range->high);
          Summary fj = RangeSummary(Difference(setMB->Get(Tag("OO_vac", r.emb), iR, kJetClass, kProbe), setMB->Get(Tag("OO_vac", ""), iR, kJetClass, kProbe), "tFJ"), range->low, range->high);
          Summary fp = RangeSummary(Difference(setMB->Get(Tag("py_pp", r.emb), iR, kJetClass, kProbe), setMB->Get(Tag("py_pp", ""), iR, kJetClass, kProbe), "tFP"), range->low, range->high);
          if (!s.ok() || !fj.ok() || !fp.ok()) { tex << " & -- & --"; continue; }
          double sigmaF2 = 0.25 * (fj.error * fj.error + fp.error * fp.error), deltaF = 0.5 * std::fabs(fj.value - fp.value);
          double z = s.value / std::sqrt(s.error * s.error + sigmaF2 + deltaF * deltaF);
          tex << Form(" & $%+.3f\\pm%.3f$ & %.1f", s.value, s.error, z);
        }
      }
      tex << " \\\\\n";
    }
  }
  printf("figures and table -> %s\n", outDir);
}
