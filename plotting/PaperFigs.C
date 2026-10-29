// ============================================================================
//  PaperFigs.C -- figures of the jet--UE v0(pT) model paper (draft v1), made
//  from the Phase2 outputs of the first (coarse) production.
//
//  Nominal choices in every figure:
//    reference  : particles outside the jet and the perpendicular cones (refnom)
//    probe      : f_sub (perpendicular-cone subtracted), unless stated otherwise
//    acceptance : R <= 0.4 with |eta_trk| < 0.9 (acc0), R = 0.8 with
//                 |eta_trk| < 2.5 (acc1)
//
//  Figures (PDF and PNG):
//    fig_pp       v0(pT) in pp: JEWEL vacuum and PYTHIA 8 with MPI on/off,
//                 R = 0.4 and 0.8, jets 20-40 GeV/c
//    fig_dv0_R    <Delta v0> vs R for soft (1-6 GeV/c) and hard (12-100 GeV/c)
//                 probes, jets 40-60 GeV/c: recoils - vacuum and no recoils -
//                 vacuum (both with O+O MB UE embedded), and the UE-only change
//                 of JEWEL vacuum and of PYTHIA 8
//    fig_methods  significance Z of recoils - vacuum vs R for v0 with f, f_sub,
//                 f_cs and for the covariance term only (f_sub), jets 40-60
//    fig_dv0_pt   Delta v0(pT) at R = 0.4 and 0.8, jets 40-60, same curves as
//                 fig_dv0_R
//
//  Input: <resultsDir>/refnom/acc<i>/v0pt_<file>.root (one Phase2 file per
//  comparison), with TH1D objects (x = probe pT)
//    v0<probe>_<tag>_c1_R<iR>_<jetClass>    v0(pT), bootstrap errors
//    cov<probe>_<tag>_c1_R<iR>_<jetClass>   covariance term of v0 only
//  <probe> = "" (raw f), "_fsub", "_fcs"; iR = 0..3 for R = 0.1, 0.2, 0.4, 0.8;
//  jetClass e.g. jet40_60 (UE-corrected jet pT). Differences between two
//  samples of the same file use the quadrature sum of their errors.
//
//  Run: root -l -b -q 'PaperFigs.C+("/data/alice/sbhatta/JewelV0pT/results_coarse","paper_figs")'
// ============================================================================
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <algorithm>
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"

static TString gResultsDir;
// radius values and the acceptance index used for each radius (see header)
static const int kAcceptanceOfRadius[4] = {0, 0, 0, 1};
static const double kRadius[4] = {0.1, 0.2, 0.4, 0.8};
static int gUniqueId = 0;   // suffix for unique object names

// probe-pT ranges of the range averages: soft and hard probes
static const double kProbeRange[2][2] = {{1, 6}, {12, 100}};
static const char* kProbeRangeLabel[2] = {"1 < #it{p}_{T}^{probe} < 6 GeV/#it{c}", "12 < #it{p}_{T}^{probe} < 100 GeV/#it{c}"};

// one histogram from <resultsDir>/refnom/acc<acc>/v0pt_<file>.root (a detached copy), or 0
static TH1D* LoadResult(int acc, const char* file, const char* name)
{
  TFile* f = TFile::Open(Form("%s/refnom/acc%d/v0pt_%s.root", gResultsDir.Data(), acc, file));
  if (!f || f->IsZombie()) return 0;
  TH1D* h = (TH1D*)f->Get(name);
  TH1D* copy = h ? (TH1D*)h->Clone(Form("h%d", gUniqueId++)) : 0;
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

// Difference of quantity q ("v0" or "cov") with probe variant (e.g. "_fsub")
// between two samples, tagPlus - tagMinus, stored in the same Phase2 file
static TH1D* DifferenceOfSamples(int acc, const char* file, const char* quantity, const char* probe,
                                 const char* tagPlus, const char* tagMinus, int iR, const char* jetClass)
{
  return Difference(LoadResult(acc, file, Form("%s%s_%s_c1_R%d_%s", quantity, probe, tagPlus, iR, jetClass)),
                    LoadResult(acc, file, Form("%s%s_%s_c1_R%d_%s", quantity, probe, tagMinus, iR, jetClass)));
}

// Range average over probe pT: inverse-variance weighted mean of the bins whose
// lower edge lies in [ptLow, ptHigh); bins without error or with error > 0.5 are
// skipped. Returns false if no bin contributes.
static bool RangeAverage(TH1D* h, double ptLow, double ptHigh, double& mean, double& error)
{
  if (!h) return false;
  double sumW = 0, sumWX = 0;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double lowEdge = h->GetXaxis()->GetBinLowEdge(i), err = h->GetBinError(i);
    if (lowEdge < ptLow - 1e-6 || lowEdge >= ptHigh - 1e-6 || err <= 0 || err > 0.5) continue;
    sumW += 1 / (err * err); sumWX += h->GetBinContent(i) / (err * err);
  }
  if (sumW <= 0) return false;
  mean = sumWX / sumW; error = 1 / sqrt(sumW); return true;
}

static void SetPadMargins(double left = 0.15, double bottom = 0.14)
{ gPad->SetTicks(1, 1); gPad->SetLeftMargin(left); gPad->SetBottomMargin(bottom); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.04); }

static void StyleCurve(TAttLine* line, TAttMarker* marker, int color, int markerStyle, int lineStyle = 1)
{
  line->SetLineColor(color); line->SetLineWidth(2); line->SetLineStyle(lineStyle);
  marker->SetMarkerColor(color); marker->SetMarkerStyle(markerStyle); marker->SetMarkerSize(1.2);
}

// y range that follows the drawn points: the data fill the pad between a bottom
// and a top fraction reserved for labels and legends; zero is always included.
// growExtent widens [lo, hi] to the drawn points (hidden bins are ignored).
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

// empty axis frame with titles and a dashed line at y = 0
static TH1D* DrawFrame(double xMin, double xMax, double yMin, double yMax, const char* xTitle, const char* yTitle, bool logX)
{
  TH1D* f = new TH1D(Form("f%d", gUniqueId++), "", 10, xMin, xMax); f->SetDirectory(0);
  f->GetYaxis()->SetRangeUser(yMin, yMax); f->GetXaxis()->SetTitle(xTitle); f->GetYaxis()->SetTitle(yTitle);
  f->GetXaxis()->SetTitleSize(0.055); f->GetYaxis()->SetTitleSize(0.055); f->GetXaxis()->SetLabelSize(0.05); f->GetYaxis()->SetLabelSize(0.05);
  f->GetYaxis()->SetTitleOffset(1.25);
  if (logX) { f->GetXaxis()->SetMoreLogLabels(); f->GetXaxis()->SetNoExponent(); }
  f->Draw(); if (logX) gPad->SetLogx();
  TLine* zero = new TLine(xMin, 0, xMax, 0); zero->SetLineStyle(2); zero->SetLineColor(kGray + 1); zero->Draw();
  return f;
}

// hide bins with no error or error > maxErr (content -999, below every frame)
static TH1D* HideUnreliableBins(TH1D* h, double maxErr = 0.3)
{
  if (!h) return 0;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinError(i) <= 0 || h->GetBinError(i) > maxErr) { h->SetBinContent(i, -999); h->SetBinError(i, 0); }
  return h;
}

static TLegend* MakeLegend(double x0, double y0, double x1, double y1, double textSize)
{
  TLegend* l = new TLegend(x0, y0, x1, y1); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(textSize);
  return l;
}

// one Delta v0 curve: tagPlus - tagMinus, both stored in the Phase2 file `file`
struct DeltaCurve { const char* file; const char* tagPlus; const char* tagMinus; const char* label; int color, marker; };

void PaperFigs(const char* resultsDir = "/data/alice/sbhatta/JewelV0pT/results_coarse", const char* outDir = "paper_figs")
{
  gResultsDir = resultsDir; gSystem->mkdir(outDir, true);
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetEndErrorSize(0); gStyle->SetLegendFont(42); gStyle->SetTextFont(42);
  TLatex text; text.SetNDC(); text.SetTextFont(42);
  const int radiusIndexOfPanel[2] = {2, 3};   // panels at R = 0.4 and R = 0.8

  // ---------- fig_pp: v0(pT) in pp, JEWEL vacuum and PYTHIA 8 ----------
  {
    TCanvas c("cpp", "", 1200, 520); c.Divide(2, 1, 0.001, 0.001);
    for (int panel = 0; panel < 2; panel++) {
      c.cd(panel + 1); SetPadMargins();
      int iR = radiusIndexOfPanel[panel], acc = kAcceptanceOfRadius[iR];
      TH1D* frame = DrawFrame(0.2, 100, -0.3, 1.0, "#it{p}_{T}^{probe} (GeV/#it{c})", "#it{v}_{0}(#it{p}_{T})^{jet-UE}", true);
      double yLow = 1e9, yHigh = -1e9;
      const char* tags[3] = {"pp_vac_nom", "py_pp_noMPI_nom", "py_pp_nom"};
      const char* labels[3] = {"JEWEL (vacuum)", "PYTHIA 8, MPI off", "PYTHIA 8, MPI on"};
      int colors[3] = {kBlack, kAzure + 2, kRed + 1}, markers[3] = {20, 25, 21};
      TLegend* legend = MakeLegend(0.17, 0.50, 0.55, 0.73, 0.045);
      for (int s = 0; s < 3; s++) {
        TH1D* h = HideUnreliableBins(LoadResult(acc, "pp_main", Form("v0_fsub_%s_c1_R%d_jet20_40", tags[s], iR)));
        if (!h) continue;
        StyleCurve(h, h, colors[s], markers[s]); h->Draw("E1 X0 SAME"); legend->AddEntry(h, labels[s], "lp"); growExtent(h, yLow, yHigh);
      }
      fitYRange(frame, yLow, yHigh, 0.24, 0.04);   // labels on top; legend in the empty left middle
      legend->Draw(); text.SetTextSize(0.047);
      text.DrawLatex(0.18, 0.88, Form("pp #sqrt{#it{s}} = 5.36 TeV, #it{R} = %s", iR == 2 ? "0.4" : "0.8"));
      text.DrawLatex(0.18, 0.82, Form("20 < #it{p}_{T}^{jet,corr} < 40 GeV/#it{c}, |#it{#eta}_{trk}| < %s", acc ? "2.5" : "0.9"));
    }
    c.Print(Form("%s/fig_pp.pdf", outDir)); c.Print(Form("%s/fig_pp.png", outDir));
  }

  // ---------- fig_dv0_R: range-averaged Delta v0 vs R ----------
  {
    TCanvas c("cR", "", 1200, 520); c.Divide(2, 1, 0.001, 0.001);
    for (int panel = 0; panel < 2; panel++) {   // panel 0: soft probes, 1: hard probes
      c.cd(panel + 1); SetPadMargins();
      DrawFrame(0, 0.9, panel == 0 ? -0.06 : -0.35, panel == 0 ? 0.22 : 0.35, "#it{R}", "#LT#Delta#it{v}_{0}#GT", false);
      // signals at the same embedded UE (O+O MB), and UE-only changes (embedded - no UE)
      DeltaCurve curves[4] = {{"OO_embMB", "OO_RC_cs10_embMB", "OO_vac_embMB", "O+O: recoils #minus vacuum", kRed + 1, 20},
                              {"OO_embMB", "OO_NR_embMB", "OO_vac_embMB", "O+O: no recoils #minus vacuum", kAzure + 2, 21},
                              {"OO_vac_closure", "OO_vac_embMB", "OO_vac_nom", "UE only, JEWEL vacuum", kGray + 2, 24},
                              {"pp_py_closure", "py_pp_embMB", "py_pp_nom", "UE only, PYTHIA 8", kGreen + 2, 25}};
      TLegend* legend = MakeLegend(0.18, 0.72, 0.95, 0.93, 0.038); legend->SetNColumns(2);
      for (int k = 0; k < 4; k++) {
        TGraphErrors* g = new TGraphErrors();
        for (int iR = 0; iR < 4; iR++) {
          double mean, error;
          TH1D* delta = DifferenceOfSamples(kAcceptanceOfRadius[iR], curves[k].file, "v0", "_fsub", curves[k].tagPlus, curves[k].tagMinus, iR, "jet40_60");
          if (!RangeAverage(delta, kProbeRange[panel][0], kProbeRange[panel][1], mean, error)) continue;
          int n = g->GetN(); g->SetPoint(n, kRadius[iR] + 0.012 * (k - 1.5), mean); g->SetPointError(n, 0, error);   // small x offset per curve
        }
        StyleCurve(g, g, curves[k].color, curves[k].marker, k < 2 ? 1 : 3); g->Draw("PL SAME"); legend->AddEntry(g, curves[k].label, "lp");
      }
      legend->Draw(); text.SetTextSize(0.045);
      text.DrawLatex(0.18, 0.17, kProbeRangeLabel[panel]); text.DrawLatex(0.18, 0.23, "40 < #it{p}_{T}^{jet,corr} < 60 GeV/#it{c}, UE: O+O MB");
    }
    c.Print(Form("%s/fig_dv0_R.pdf", outDir)); c.Print(Form("%s/fig_dv0_R.png", outDir));
  }

  // ---------- fig_methods: significance Z of recoils - vacuum vs R ----------
  // Z = S / sqrt(sigma_S^2 + sigma_F^2 + ((F_J - F_P)/2)^2) with
  //   S       = range average of recoils - vacuum (same O+O MB UE),
  //   F_J/F_P = UE-only change of JEWEL vacuum / PYTHIA 8 (embedded - no UE),
  //   sigma_F^2 = (sigma_FJ^2 + sigma_FP^2)/4 (error of their mean); half their
  //   difference enters as a generator uncertainty. S itself is not corrected.
  //   All terms are range averages over the probe-pT range of the panel.
  {
    TCanvas c("cm", "", 1200, 520); c.Divide(2, 1, 0.001, 0.001);
    struct Method { const char* quantity; const char* probe; const char* label; int color, marker, lineStyle; };
    Method methods[4] = {{"v0", "", "#it{v}_{0}, raw #it{f}", kBlack, 20, 1}, {"v0", "_fsub", "#it{v}_{0}, perp. cones #it{f}_{sub}", kRed + 1, 21, 1},
                         {"v0", "_fcs", "#it{v}_{0}, constituent sub. #it{f}_{cs}", kAzure + 2, 22, 1}, {"cov", "_fsub", "covariance only, #it{f}_{sub}", kGray + 2, 24, 3}};
    for (int panel = 0; panel < 2; panel++) {   // panel 0: soft probes, 1: hard probes
      c.cd(panel + 1); SetPadMargins();
      DrawFrame(0, 0.9, -6, 11, "#it{R}", "significance #it{Z} of recoils #minus vacuum", false);
      for (double y : {-2., 2.}) { TLine* l = new TLine(0, y, 0.9, y); l->SetLineStyle(3); l->SetLineColor(kGray + 1); l->Draw(); }   // |Z| = 2
      TLegend* legend = MakeLegend(0.18, 0.74, 0.95, 0.94, 0.036); legend->SetNColumns(2);
      const double ptLow = kProbeRange[panel][0], ptHigh = kProbeRange[panel][1];
      for (int k = 0; k < 4; k++) {
        const Method& m = methods[k];
        TGraph* g = new TGraph();
        for (int iR = 0; iR < 4; iR++) {
          int acc = kAcceptanceOfRadius[iR];
          double signal, errSignal, changeJewel, errJewel, changePythia, errPythia;
          bool ok = RangeAverage(DifferenceOfSamples(acc, "OO_embMB", m.quantity, m.probe, "OO_RC_cs10_embMB", "OO_vac_embMB", iR, "jet40_60"), ptLow, ptHigh, signal, errSignal)
                 && RangeAverage(DifferenceOfSamples(acc, "OO_vac_closure", m.quantity, m.probe, "OO_vac_embMB", "OO_vac_nom", iR, "jet40_60"), ptLow, ptHigh, changeJewel, errJewel)
                 && RangeAverage(DifferenceOfSamples(acc, "pp_py_closure", m.quantity, m.probe, "py_pp_embMB", "py_pp_nom", iR, "jet40_60"), ptLow, ptHigh, changePythia, errPythia);
          if (!ok) continue;
          double varBackground = 0.25 * (errJewel * errJewel + errPythia * errPythia) + 0.25 * (changeJewel - changePythia) * (changeJewel - changePythia);
          g->SetPoint(g->GetN(), kRadius[iR] + 0.01 * (k - 1.5), signal / sqrt(errSignal * errSignal + varBackground));
        }
        StyleCurve(g, g, m.color, m.marker, m.lineStyle); g->Draw("PL SAME"); legend->AddEntry(g, m.label, "lp");
      }
      legend->Draw(); text.SetTextSize(0.045);
      text.DrawLatex(0.18, 0.17, kProbeRangeLabel[panel]); text.DrawLatex(0.18, 0.23, "40 < #it{p}_{T}^{jet,corr} < 60 GeV/#it{c}, UE: O+O MB");
    }
    c.Print(Form("%s/fig_methods.pdf", outDir)); c.Print(Form("%s/fig_methods.png", outDir));
  }

  // ---------- fig_dv0_pt: Delta v0(pT) at R = 0.4 and 0.8 ----------
  {
    TCanvas c("cp", "", 1200, 520); c.Divide(2, 1, 0.001, 0.001);
    for (int panel = 0; panel < 2; panel++) {
      c.cd(panel + 1); SetPadMargins();
      int iR = radiusIndexOfPanel[panel], acc = kAcceptanceOfRadius[iR];
      // x axis starts at 0.2 GeV/c so that every probe bin (first bin 0.2-0.5) lies inside the frame
      TH1D* frame = DrawFrame(0.2, 100, -0.4, 0.4, "#it{p}_{T}^{probe} (GeV/#it{c})", "#Delta#it{v}_{0}(#it{p}_{T})", true);
      double yLow = 1e9, yHigh = -1e9;
      DeltaCurve curves[4] = {{"OO_embMB", "OO_RC_cs10_embMB", "OO_vac_embMB", "recoils #minus vacuum", kRed + 1, 20},
                              {"OO_embMB", "OO_NR_embMB", "OO_vac_embMB", "no recoils #minus vacuum", kAzure + 2, 21},
                              {"OO_vac_closure", "OO_vac_embMB", "OO_vac_nom", "UE only, JEWEL vacuum", kGray + 2, 24},
                              {"pp_py_closure", "py_pp_embMB", "py_pp_nom", "UE only, PYTHIA 8", kGreen + 2, 25}};
      TLegend* legend = MakeLegend(0.18, 0.16, 0.60, 0.40, 0.042);
      for (int k = 0; k < 4; k++) {
        TH1D* h = HideUnreliableBins(DifferenceOfSamples(acc, curves[k].file, "v0", "_fsub", curves[k].tagPlus, curves[k].tagMinus, iR, "jet40_60"));
        if (!h) continue;
        // a TH1 cannot be shifted horizontally: draw the visible bins as a graph, offset by a
        // few percent in pT per curve so that the error bars do not overlap
        TGraphErrors* g = new TGraphErrors();
        for (int i = 1; i <= h->GetNbinsX(); i++) {
          if (h->GetBinContent(i) < -900) continue;
          double x = h->GetXaxis()->GetBinCenter(i) * (1 + 0.04 * (k - 1.5));
          int n = g->GetN(); g->SetPoint(n, x, h->GetBinContent(i)); g->SetPointError(n, 0, h->GetBinError(i));
        }
        StyleCurve(g, g, curves[k].color, curves[k].marker); g->Draw("P SAME"); legend->AddEntry(g, curves[k].label, "p"); growExtent(g, yLow, yHigh);
      }
      fitYRange(frame, yLow, yHigh, 0.22, 0.30);   // labels on top, legend below the data
      legend->Draw(); text.SetTextSize(0.045);
      text.DrawLatex(0.18, 0.88, Form("O+O 5.36 TeV, UE: O+O MB, #it{R} = %s", iR == 2 ? "0.4" : "0.8"));
      text.DrawLatex(0.18, 0.82, "40 < #it{p}_{T}^{jet,corr} < 60 GeV/#it{c}, probe #it{f}_{sub}");
    }
    c.Print(Form("%s/fig_dv0_pt.pdf", outDir)); c.Print(Form("%s/fig_dv0_pt.png", outDir));
  }
}
