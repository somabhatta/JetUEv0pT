// ============================================================================
//  PlotEventDisplay.C -- what jetv0ana does to one event, step by step.
//
//  Input: the event-display file written by
//    jetv0ana ... -display display.root [-displayEvents N] [-displayVariants A,B]
//  (tree "display": one entry per event, variant and (acceptance, R) slot;
//  tree "recoil": the recoil thermal-momentum subtraction; see the header of
//  generation/jetv0ana.cc). Every number drawn here was recorded inside the
//  production code path; nothing is recomputed except the drawing geometry.
//
//  For one event, one multi-page PDF and one PNG per page. Pages per displayed
//  variant (slot |eta_trk| < 0.9, anti-kT R = 0.4 unless noted):
//    (a) particles  eta-phi map of every charged track entering the clustering
//                   (pT > 0.15 GeV/c, |eta| < 0.9). Marker AREA is proportional
//                   to pT; colour = origin: generator, added UE, or generator
//                   particle whose pT was reduced by the recoil thermal-momentum
//                   subtraction (JEWEL recoil samples).
//    (b) jets       anti-kT jets with pT > 3 GeV/c: constituents coloured by jet,
//                   the jet's active-area ghosts drawn as cells of the ghost area
//                   (0.05), so the coloured patch is the jet area. Acceptance
//                   |eta| < 0.9 (solid) and the leading-jet fiducial band
//                   |eta_jet| < 0.9 - R (shaded); the selected leading jet (the
//                   highest-pT jet inside the band) is circled. Page (b2) is the
//                   same for R = 0.8 in |eta_trk| < 2.5.
//    (c) background left: jet cone and the two perpendicular cones at
//                   phi_jet +- pi/2 (radius R); tracks coloured by role: leading-
//                   jet constituent, in a perpendicular cone, reference track
//                   (0.2 <= pT < 3 GeV/c outside jet and perpendicular cones),
//                   reference track removed by the away-side cut |dphi - pi| < 1
//                   (refaw), other. Black rings mark the eta-gap reference
//                   (|eta - eta_jet| > R + 0.5, dashed lines). Right: per probe-
//                   pT bin the in-cone count n_b, the perpendicular-cone estimate
//                   b_b (mean of the two cones) and n_b - b_b; below, the
//                   fractions f_b = n_b / N and f_sub = (n_b - b_b) / (N - B).
//    (d) constituent subtraction (event-wise CS, dR_max 0.25, alpha 0, ghost
//                   area 0.01, rho from the R = 0.4 perpendicular cones): left,
//                   track pT after vs before CS (removed tracks on the bottom
//                   line); middle, eta vs phi - phi_jet around the leading jet
//                   with the CS ghost grid, tracks before (open) and after
//                   (filled) CS; right, jet pT raw, rho*A-corrected, CS and
//                   truth (same event without added UE).
//    (e) recoils    (JEWEL recoil samples only, dR_max = 1.0) rapidity-phi map of
//                   the final-state particles and the thermal scattering centres
//                   before (open) and after (filled: unsubtracted remainder) the
//                   thermal-momentum subtraction, and their pT spectra; the table
//                   gives the pT bookkeeping for dR_max = 1.0 and 0.5.
//  Each page has a title line (sample, event, variant, slot) and a one-line
//  caption at the bottom. Colours: Okabe-Ito (colour-blind safe).
//
//  Run (LCG_107 view), from this directory:
//    root -l -b -q 'PlotEventDisplay.C+("display.root","plots","PYTHIA 8 pp","py")'
//  Arguments: display file, output directory, sample label (title line), file
//  tag, event number (-1: the displayed event with the highest R = 0.4 leading-
//  jet pT in the first variant), comma-separated variants ("": all in the file).
// ============================================================================
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
//  Colours (Okabe-Ito) and roles
// ---------------------------------------------------------------------------
int Color(const char* hex) { return TColor::GetColor(hex); }
const char* kBlack = "#000000";  const char* kOrange = "#E69F00"; const char* kSkyBlue = "#56B4E9";
const char* kGreen = "#009E73";  const char* kYellow = "#F0E442"; const char* kBlue = "#0072B2";
const char* kVermillion = "#D55E00"; const char* kPurple = "#CC79A7"; const char* kGrey = "#999999";
// jet colours in pT order (yellow left out: too light on white)
const char* kJetColor[6] = {"#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"};
// track origin (jetv0ana trkOrigin)
const char* kOriginColor[3] = {"#0072B2", "#E69F00", "#009E73"};
const char* kOriginLabel[3] = {"generator", "added UE", "recoil thermal-subtracted"};
// track flags (jetv0ana kDisplayFlag*)
enum { kLeadJet = 1, kJetCone = 2, kPerpCone = 4, kRef = 8, kRefAway = 16, kGap = 32, kCSJet = 64 };

// marker AREA proportional to pT (ROOT marker size scales the diameter)
double SizeForPt(double pt) { return std::min(5.0, std::max(0.25, 0.45 * std::sqrt(pt))); }

double WrapPhi(double phi) {
  while (phi > TMath::Pi()) phi -= 2 * TMath::Pi();
  while (phi <= -TMath::Pi()) phi += 2 * TMath::Pi();
  return phi;
}

std::string VariantDescription(const std::string& v) {
  std::string base = v, ue;
  if (v.rfind("cs10_", 0) == 0) { base = "cs10"; ue = v.substr(5); }
  else if (v.rfind("emb", 0) == 0) { base = "nom"; ue = v; }
  std::map<std::string, std::string> baseText = {
    {"nom", "generator particles"}, {"nosub", "recoils, thermal momenta NOT subtracted"},
    {"cs10", "recoils, thermal momenta subtracted (#DeltaR_{max} = 1.0)"},
    {"cs05", "recoils, thermal momenta subtracted (#DeltaR_{max} = 0.5)"}};
  std::map<std::string, std::string> ueText = {
    {"embMB", "thermal UE dN/d#eta = 42"}, {"embC", "thermal UE dN/d#eta = 135"},
    {"embNp", "thermal UE N_{part}(c)-scaled"}, {"embFlow", "flow-correlated thermal UE"},
    {"embAng", "Angantyr UE"}, {"embAngC", "Angantyr UE (central)"}, {"embHij", "HIJING UE"}};
  std::string text = baseText.count(base) ? baseText[base] : base;
  if (!ue.empty()) text += " + " + (ueText.count(ue) ? ueText[ue] : ue);
  return text;
}

// ---------------------------------------------------------------------------
//  Display records
// ---------------------------------------------------------------------------
struct SlotEntry {
  int event = 0, iacc = 0, iR = 0, leadJet = -1;
  std::string variant;
  float R = 0, trackEtaMax = 0, jetEtaMax = 0, weight = 0, centr = 0, ueMult = 0;
  float leadPt = 0, leadPtCorr = 0, leadEta = 0, leadPhi = 0, leadArea = 0, rho = 0, rhoCS = 0;
  float csJetPt = 0, csJetEta = 0, csJetPhi = 0, jetPtTruth = 0, refN = 0, refS1 = 0;
  float jetGhostArea = 0, csGhostArea = 0, csMaxDR = 0, perpPhi[2] = {0, 0};
  std::vector<float> probeEdges, nIn, nPerp, nSub, fRaw, fSub;
  std::vector<float> trkPt, trkEta, trkPhi, trkPtCS, jetPt, jetEta, jetPhi, jetArea;
  std::vector<float> ghostEta, ghostPhi, csGhostEta, csGhostPhi;
  std::vector<int> trkOrigin, trkJet, trkFlags, jetNConst, jetFiducial, ghostJet;
};

struct RecoilEntry {
  int event = 0; float dRmax = 0; Long64_t nPairs = 0;
  std::vector<float> finPt, finY, finPhi, finPtAfter, thPt, thY, thPhi, thPtAfter;
  std::vector<int> finCharged;
};

template <class T> void CopyVector(std::vector<T>* from, std::vector<T>& to) { if (from) to = *from; }

std::vector<SlotEntry> ReadSlots(TFile& f) {
  std::vector<SlotEntry> out;
  TTree* t = (TTree*)f.Get("display");
  if (!t) return out;
  SlotEntry e;
  std::string* variant = nullptr;
  float probeEdges[26], nIn[25], nPerp[25], nSub[25], fRaw[25], fSub[25];
  std::map<std::string, std::vector<float>*> fv;
  std::map<std::string, std::vector<int>*> iv;
  for (auto n : {"trkPt", "trkEta", "trkPhi", "trkPtCS", "jetPt", "jetEta", "jetPhi", "jetArea", "ghostEta",
                 "ghostPhi", "csGhostEta", "csGhostPhi"}) { fv[n] = nullptr; t->SetBranchAddress(n, &fv[n]); }
  for (auto n : {"trkOrigin", "trkJet", "trkFlags", "jetNConst", "jetFiducial", "ghostJet"})
    { iv[n] = nullptr; t->SetBranchAddress(n, &iv[n]); }
  t->SetBranchAddress("variant", &variant);
  t->SetBranchAddress("event", &e.event); t->SetBranchAddress("iacc", &e.iacc); t->SetBranchAddress("iR", &e.iR);
  t->SetBranchAddress("leadJet", &e.leadJet);
  std::map<std::string, float*> fs = {
    {"R", &e.R}, {"trackEtaMax", &e.trackEtaMax}, {"jetEtaMax", &e.jetEtaMax}, {"weight", &e.weight},
    {"centr", &e.centr}, {"ueMult", &e.ueMult}, {"leadPt", &e.leadPt}, {"leadPtCorr", &e.leadPtCorr},
    {"leadEta", &e.leadEta}, {"leadPhi", &e.leadPhi}, {"leadArea", &e.leadArea}, {"rho", &e.rho},
    {"rhoCS", &e.rhoCS}, {"csJetPt", &e.csJetPt}, {"csJetEta", &e.csJetEta}, {"csJetPhi", &e.csJetPhi},
    {"jetPtTruth", &e.jetPtTruth}, {"refN", &e.refN}, {"refS1", &e.refS1}, {"jetGhostArea", &e.jetGhostArea},
    {"csGhostArea", &e.csGhostArea}, {"csMaxDR", &e.csMaxDR}};
  for (auto& b : fs) t->SetBranchAddress(b.first.c_str(), b.second);
  t->SetBranchAddress("perpPhi", e.perpPhi);
  t->SetBranchAddress("probeEdges", probeEdges); t->SetBranchAddress("nIn", nIn); t->SetBranchAddress("nPerp", nPerp);
  t->SetBranchAddress("nSub", nSub); t->SetBranchAddress("fRaw", fRaw); t->SetBranchAddress("fSub", fSub);
  for (Long64_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    SlotEntry s = e;
    s.variant = *variant;
    s.probeEdges.assign(probeEdges, probeEdges + 26); s.nIn.assign(nIn, nIn + 25); s.nPerp.assign(nPerp, nPerp + 25);
    s.nSub.assign(nSub, nSub + 25); s.fRaw.assign(fRaw, fRaw + 25); s.fSub.assign(fSub, fSub + 25);
    CopyVector(fv["trkPt"], s.trkPt); CopyVector(fv["trkEta"], s.trkEta); CopyVector(fv["trkPhi"], s.trkPhi);
    CopyVector(fv["trkPtCS"], s.trkPtCS); CopyVector(fv["jetPt"], s.jetPt); CopyVector(fv["jetEta"], s.jetEta);
    CopyVector(fv["jetPhi"], s.jetPhi); CopyVector(fv["jetArea"], s.jetArea); CopyVector(fv["ghostEta"], s.ghostEta);
    CopyVector(fv["ghostPhi"], s.ghostPhi); CopyVector(fv["csGhostEta"], s.csGhostEta);
    CopyVector(fv["csGhostPhi"], s.csGhostPhi);
    CopyVector(iv["trkOrigin"], s.trkOrigin); CopyVector(iv["trkJet"], s.trkJet); CopyVector(iv["trkFlags"], s.trkFlags);
    CopyVector(iv["jetNConst"], s.jetNConst); CopyVector(iv["jetFiducial"], s.jetFiducial);
    CopyVector(iv["ghostJet"], s.ghostJet);
    out.push_back(s);
  }
  return out;
}

std::vector<RecoilEntry> ReadRecoils(TFile& f) {
  std::vector<RecoilEntry> out;
  TTree* t = (TTree*)f.Get("recoil");
  if (!t) return out;
  RecoilEntry e;
  std::map<std::string, std::vector<float>*> fv;
  for (auto n : {"finPt", "finY", "finPhi", "finPtAfter", "thPt", "thY", "thPhi", "thPtAfter"})
    { fv[n] = nullptr; t->SetBranchAddress(n, &fv[n]); }
  std::vector<int>* charged = nullptr;
  t->SetBranchAddress("finCharged", &charged);
  t->SetBranchAddress("event", &e.event); t->SetBranchAddress("dRmax", &e.dRmax); t->SetBranchAddress("nPairs", &e.nPairs);
  for (Long64_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    RecoilEntry r = e;
    CopyVector(fv["finPt"], r.finPt); CopyVector(fv["finY"], r.finY); CopyVector(fv["finPhi"], r.finPhi);
    CopyVector(fv["finPtAfter"], r.finPtAfter); CopyVector(fv["thPt"], r.thPt); CopyVector(fv["thY"], r.thY);
    CopyVector(fv["thPhi"], r.thPhi); CopyVector(fv["thPtAfter"], r.thPtAfter); CopyVector(charged, r.finCharged);
    out.push_back(r);
  }
  return out;
}

// ---------------------------------------------------------------------------
//  Page layout: title strip, body, caption strip
// ---------------------------------------------------------------------------
struct Page {
  TCanvas* canvas = nullptr;
  TPad* body = nullptr;
};

int gPageCounter = 0;

// shrink a one-line text until it is at most maxWidth (pad fraction) wide; the
// width is measured at a quarter of the size (the measurement saturates at the
// pad width) and scaled, with a 3% margin for the rendering
void FitWidth(TLatex* text, double maxWidth) {
  const double size = text->GetTextSize();
  text->SetTextSize(0.25 * size);
  double w = 4 * text->GetXsize();
  text->SetTextSize(w > maxWidth ? 0.97 * size * maxWidth / w : size);
}

Page NewPage(const std::string& title, const std::string& caption, int width = 1600, int height = 900) {
  Page p;
  p.canvas = new TCanvas(Form("page%d", gPageCounter++), "", width, height);
  const double titleH = 0.065, captionH = 0.055;
  TPad* top = new TPad(Form("title%d", gPageCounter), "", 0, 1 - titleH, 1, 1);
  TPad* bottom = new TPad(Form("caption%d", gPageCounter), "", 0, 0, 1, captionH);
  p.body = new TPad(Form("body%d", gPageCounter), "", 0, captionH, 1, 1 - titleH);
  top->Draw(); bottom->Draw(); p.body->Draw();
  top->cd();
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.42); t.SetTextAlign(12);
  FitWidth(t.DrawLatex(0.01, 0.5, title.c_str()), 0.98);
  bottom->cd();
  TLatex c; c.SetNDC(); c.SetTextFont(42); c.SetTextSize(0.42); c.SetTextAlign(12); c.SetTextColor(Color("#333333"));
  FitWidth(c.DrawLatex(0.01, 0.5, caption.c_str()), 0.98);
  p.body->cd();
  return p;
}

TPad* SubPad(TPad* body, const char* name, double x1, double y1, double x2, double y2) {
  body->cd();
  TPad* p = new TPad(Form("%s_%d", name, gPageCounter), "", x1, y1, x2, y2);
  p->Draw(); p->cd();
  return p;
}

// eta-phi frame: x = phi, y = eta
TH2F* EtaPhiFrame(double etaRange, const char* xTitle = "#varphi (rad)", double xMin = -TMath::Pi(),
                  double xMax = TMath::Pi(), const char* yTitle = "#eta", double yTitleOffset = 0.75) {
  static int nFrames = 0;
  TH2F* h = new TH2F(Form("frame%d", nFrames++), "", 10, xMin, xMax, 10,
                     -etaRange, etaRange);
  h->SetDirectory(nullptr);
  h->GetXaxis()->SetTitle(xTitle); h->GetYaxis()->SetTitle(yTitle);
  h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetTitleOffset(yTitleOffset); h->GetYaxis()->CenterTitle(); h->GetXaxis()->CenterTitle();
  h->SetStats(0);
  h->Draw("AXIS");
  return h;
}

void HLine(double x1, double x2, double y, int color, int style, int width) {
  TLine* l = new TLine(x1, y, x2, y); l->SetLineColor(color); l->SetLineStyle(style); l->SetLineWidth(width); l->Draw();
}

// acceptance (solid) and fiducial jet band (shaded, dashed edges)
void DrawAcceptance(double trackEtaMax, double jetEtaMax, bool fiducialBand) {
  const double pi = TMath::Pi();
  if (fiducialBand) {
    TBox* band = new TBox(-pi, -jetEtaMax, pi, jetEtaMax);
    band->SetFillColorAlpha(Color(kYellow), 0.18); band->SetLineWidth(0); band->Draw();
    HLine(-pi, pi, jetEtaMax, Color(kGrey), 2, 2); HLine(-pi, pi, -jetEtaMax, Color(kGrey), 2, 2);
  }
  HLine(-pi, pi, trackEtaMax, Color(kBlack), 1, 2); HLine(-pi, pi, -trackEtaMax, Color(kBlack), 1, 2);
}

// circle of radius R in (phi, eta), wrapped in phi and clipped to |eta| < etaClip
void DrawCone(double eta, double phi, double R, int color, int style, int width, double etaClip) {
  std::vector<double> xs, ys;
  auto flush = [&]() {
    if (xs.size() > 1) {
      TGraph* g = new TGraph(xs.size(), xs.data(), ys.data());
      g->SetLineColor(color); g->SetLineStyle(style); g->SetLineWidth(width); g->Draw("L same");
    }
    xs.clear(); ys.clear();
  };
  for (int i = 0; i <= 240; i++) {
    double a = 2 * TMath::Pi() * i / 240, x = WrapPhi(phi + R * std::cos(a)), y = eta + R * std::sin(a);
    if (std::fabs(y) > etaClip || (!xs.empty() && std::fabs(x - xs.back()) > TMath::Pi())) flush();
    if (std::fabs(y) <= etaClip) { xs.push_back(x); ys.push_back(y); }
  }
  flush();
}

TMarker* Marker(double x, double y, int style, double size, int color, double alpha = 1.0) {
  TMarker* m = new TMarker(x, y, style);
  m->SetMarkerSize(size);
  if (alpha < 1) m->SetMarkerColorAlpha(color, alpha); else m->SetMarkerColor(color);
  m->Draw();
  return m;
}

// legend pad helpers: a marker + text row at NDC height y
void LegendRow(double y, int style, double size, int color, const std::string& text, double textSize = 0.055) {
  TMarker* m = new TMarker(0.10, y, style); m->SetNDC(); m->SetMarkerSize(size); m->SetMarkerColor(color); m->Draw();
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(textSize); t.SetTextAlign(12);
  t.DrawLatex(0.22, y, text.c_str());
}
void LegendLine(double y, int color, int style, int width, const std::string& text, double textSize = 0.055) {
  TLine* l = new TLine(0.04, y, 0.17, y); l->SetNDC(); l->SetLineColor(color); l->SetLineStyle(style);
  l->SetLineWidth(width); l->Draw();
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(textSize); t.SetTextAlign(12);
  t.DrawLatex(0.22, y, text.c_str());
}
void LegendBox(double y, int color, double alpha, const std::string& text, double textSize = 0.055) {
  TBox* b = new TBox(0.05, y - 0.012, 0.16, y + 0.012); b->SetFillColorAlpha(color, alpha); b->Draw();
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(textSize); t.SetTextAlign(12);
  t.DrawLatex(0.22, y, text.c_str());
}
void LegendText(double y, const std::string& text, double textSize = 0.055, int font = 42) {
  TLatex t; t.SetNDC(); t.SetTextFont(font); t.SetTextSize(textSize); t.SetTextAlign(12);
  t.DrawLatex(0.04, y, text.c_str());
}
// pT key for marker area
double DrawPtKey(double y, double textSize = 0.05) {
  LegendText(y, "marker area #propto p_{T}:", textSize, 62);
  y -= 0.045;
  for (double pt : {0.5, 2.0, 10.0, 40.0}) {
    LegendRow(y, 24, SizeForPt(pt), Color(kBlack), Form("%g GeV/#it{c}", pt), textSize);
    y -= 0.022 + 0.010 * SizeForPt(pt);
  }
  return y;
}

std::string SlotText(const SlotEntry& s) {
  return Form("|#eta_{trk}| < %.1f, anti-k_{T} R = %.1f", s.trackEtaMax, s.R);
}

// ---------------------------------------------------------------------------
//  (a) particles
// ---------------------------------------------------------------------------
void PageParticles(const SlotEntry& s, const std::string& head, const std::string& pdf, const std::string& png) {
  Page p = NewPage(head + ", " + SlotText(s),
                   "(a) Charged tracks entering the jet clustering (p_{T} > 0.15 GeV/#it{c}); marker area #propto p_{T}, "
                   "colour = origin of the particle.");
  TPad* map = SubPad(p.body, "map", 0.0, 0.0, 0.78, 1.0);
  map->SetLeftMargin(0.085); map->SetRightMargin(0.02); map->SetBottomMargin(0.12); map->SetTopMargin(0.04);
  EtaPhiFrame(s.trackEtaMax + 0.12);
  DrawAcceptance(s.trackEtaMax, s.jetEtaMax, false);
  int nOrigin[3] = {0, 0, 0}; double ptOrigin[3] = {0, 0, 0};
  std::vector<size_t> order(s.trkPt.size());
  for (size_t i = 0; i < order.size(); i++) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return s.trkPt[a] > s.trkPt[b]; });   // small on top
  for (size_t i : order) {
    int o = std::min(2, std::max(0, s.trkOrigin[i]));
    nOrigin[o]++; ptOrigin[o] += s.trkPt[i];
    Marker(s.trkPhi[i], s.trkEta[i], 20, SizeForPt(s.trkPt[i]), Color(kOriginColor[o]), 0.75);
  }
  TPad* leg = SubPad(p.body, "leg", 0.78, 0.0, 1.0, 1.0);
  double y = 0.94;
  LegendText(y, "origin (N, #Sigma p_{T})", 0.06, 62); y -= 0.06;
  for (int o = 0; o < 3; o++) {
    if (nOrigin[o] == 0) continue;
    LegendRow(y, 20, 1.6, Color(kOriginColor[o]), kOriginLabel[o]); y -= 0.045;
    LegendText(y, Form("     N = %d, #Sigma p_{T} = %.1f GeV/#it{c}", nOrigin[o], ptOrigin[o]), 0.05); y -= 0.06;
  }
  LegendLine(y, Color(kBlack), 1, 2, Form("|#eta| = %.1f acceptance", s.trackEtaMax)); y -= 0.08;
  y = DrawPtKey(y);
  y -= 0.04;
  LegendText(y, Form("event weight %.3g", s.weight), 0.05); y -= 0.045;
  if (s.centr >= 0) { LegendText(y, Form("centrality %.1f%%", s.centr), 0.05); y -= 0.045; }
  if (s.ueMult > 0) { LegendText(y, Form("UE dN/d#eta (or N_{ch}) = %.0f", s.ueMult), 0.05); y -= 0.045; }
  p.canvas->cd(); p.canvas->Print(pdf.c_str()); p.canvas->Print((png + "_a_particles.png").c_str());
}

// ---------------------------------------------------------------------------
//  (b) jets and active areas
// ---------------------------------------------------------------------------
void PageJets(const SlotEntry& s, const std::string& head, const std::string& pdf, const std::string& png,
              const std::string& pageTag) {
  Page p = NewPage(head + ", " + SlotText(s),
                   Form("(b) Anti-k_{T} R = %.1f jets with p_{T} > 3 GeV/#it{c}: constituents and active-area ghosts "
                        "(cells of area %.2f) coloured by jet; selected leading jet (|#eta_{jet}| < %.1f, shaded) circled.",
                        s.R, s.jetGhostArea, s.jetEtaMax));
  TPad* map = SubPad(p.body, "map", 0.0, 0.0, 0.74, 1.0);
  map->SetLeftMargin(0.085); map->SetRightMargin(0.02); map->SetBottomMargin(0.12); map->SetTopMargin(0.04);
  EtaPhiFrame(s.trackEtaMax + 0.12);
  DrawAcceptance(s.trackEtaMax, s.jetEtaMax, true);
  const int nJets = s.jetPt.size();
  auto jetColor = [&](int j) { return j < 6 ? Color(kJetColor[j]) : Color(kGrey); };
  const double half = 0.5 * std::sqrt(s.jetGhostArea);
  for (size_t g = 0; g < s.ghostEta.size(); g++) {
    double x1 = std::max(-TMath::Pi(), (double)s.ghostPhi[g] - half), x2 = std::min(TMath::Pi(), (double)s.ghostPhi[g] + half);
    TBox* b = new TBox(x1, s.ghostEta[g] - half, x2, s.ghostEta[g] + half);
    b->SetFillColorAlpha(jetColor(s.ghostJet[g]), s.ghostJet[g] < 6 ? 0.28 : 0.12); b->SetLineWidth(0); b->Draw();
  }
  std::vector<size_t> order(s.trkPt.size());
  for (size_t i = 0; i < order.size(); i++) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return s.trkPt[a] > s.trkPt[b]; });
  for (size_t i : order) {
    int j = s.trkJet[i];
    if (j < 0) Marker(s.trkPhi[i], s.trkEta[i], 24, SizeForPt(s.trkPt[i]), Color(kGrey), 0.6);
    else       Marker(s.trkPhi[i], s.trkEta[i], 20, SizeForPt(s.trkPt[i]), jetColor(j), 0.9);
  }
  for (int j = 0; j < nJets; j++) Marker(s.jetPhi[j], s.jetEta[j], 5, 2.2, jetColor(j));
  if (s.leadJet >= 0) {
    DrawCone(s.jetEta[s.leadJet], s.jetPhi[s.leadJet], s.R, Color(kBlack), 1, 3, s.trackEtaMax + 0.12);
    Marker(s.jetPhi[s.leadJet], s.jetEta[s.leadJet], 29, 3.0, Color(kBlack));
  }
  TPad* leg = SubPad(p.body, "leg", 0.74, 0.0, 1.0, 1.0);
  double y = 0.95;
  LegendText(y, "jets (p_{T}, area, fiducial)", 0.05, 62); y -= 0.05;
  for (int j = 0; j < nJets && j < 9; j++) {
    LegendBox(y, jetColor(j), 0.5, Form("%s%.1f GeV/#it{c}, A = %.2f, %s", j == s.leadJet ? "#bf{lead} " : "", s.jetPt[j],
                                        s.jetArea[j], s.jetFiducial[j] ? "in band" : "outside"), 0.045);
    y -= 0.045;
  }
  if (nJets > 9) { LegendText(y, Form("... %d more jets", nJets - 9), 0.045); y -= 0.045; }
  if (nJets == 0) { LegendText(y, "no jet above 3 GeV/#it{c}", 0.045); y -= 0.045; }
  y -= 0.02;
  LegendRow(y, 20, 1.4, jetColor(0), "constituent (colour of its jet)", 0.045); y -= 0.045;
  LegendRow(y, 24, 1.4, Color(kGrey), "track in no jet above 3 GeV/#it{c}", 0.045); y -= 0.045;
  LegendRow(y, 5, 1.6, Color(kBlack), "jet axis", 0.045); y -= 0.045;
  LegendRow(y, 29, 2.0, Color(kBlack), Form("selected leading jet, circle R = %.1f", s.R), 0.045); y -= 0.045;
  LegendBox(y, Color(kYellow), 0.35, Form("fiducial band |#eta_{jet}| < %.1f", s.jetEtaMax), 0.045); y -= 0.045;
  LegendLine(y, Color(kBlack), 1, 2, Form("track acceptance |#eta| < %.1f", s.trackEtaMax), 0.045); y -= 0.06;
  y = DrawPtKey(y, 0.045);
  p.canvas->cd(); p.canvas->Print(pdf.c_str()); p.canvas->Print((png + "_" + pageTag + ".png").c_str());
}

// ---------------------------------------------------------------------------
//  (c) perpendicular cones, reference region, probe spectra
// ---------------------------------------------------------------------------
void PageBackground(const SlotEntry& s, const std::string& head, const std::string& pdf, const std::string& png) {
  Page p = NewPage(head + ", " + SlotText(s),
                   "(c) Jet cone and perpendicular cones (#varphi_{jet} #pm #pi/2), track roles; right: in-cone counts n_{b}, "
                   "perpendicular-cone estimate b_{b}, n_{b} - b_{b} and the fractions f_{b}, f_{sub} per probe-p_{T} bin.");
  TPad* map = SubPad(p.body, "map", 0.0, 0.0, 0.47, 1.0);
  map->SetLeftMargin(0.13); map->SetRightMargin(0.02); map->SetBottomMargin(0.12); map->SetTopMargin(0.04);
  const double etaRange = s.trackEtaMax + 0.12;
  EtaPhiFrame(etaRange, "#varphi (rad)", -TMath::Pi(), TMath::Pi(), "#eta", 1.05);
  DrawAcceptance(s.trackEtaMax, s.jetEtaMax, false);
  const bool haveJet = s.leadPt > 0;
  if (haveJet) {
    // away-side band removed in refaw, eta-gap boundaries
    double awayCentre = WrapPhi(s.leadPhi + TMath::Pi());
    for (double shift : {-2 * TMath::Pi(), 0.0, 2 * TMath::Pi()}) {
      double x1 = std::max(-TMath::Pi(), awayCentre - 1 + shift), x2 = std::min(TMath::Pi(), awayCentre + 1 + shift);
      if (x2 <= x1) continue;
      TBox* b = new TBox(x1, -s.trackEtaMax, x2, s.trackEtaMax); b->SetFillColorAlpha(Color(kOrange), 0.12); b->SetLineWidth(0); b->Draw();
    }
    for (double sign : {-1.0, 1.0}) {
      double yGap = s.leadEta + sign * (s.R + 0.5);
      if (std::fabs(yGap) < s.trackEtaMax) HLine(-TMath::Pi(), TMath::Pi(), yGap, Color(kBlack), 2, 1);
    }
  }
  std::vector<size_t> order(s.trkPt.size());
  for (size_t i = 0; i < order.size(); i++) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return s.trkPt[a] > s.trkPt[b]; });
  int nRole[5] = {0, 0, 0, 0, 0};
  for (size_t i : order) {
    int f = s.trkFlags[i], role;
    if (f & kLeadJet) role = 0;
    else if (f & kPerpCone) role = 1;
    else if ((f & kRef) && (f & kRefAway)) role = 2;
    else if (f & kRef) role = 3;
    else role = 4;
    nRole[role]++;
    const char* col[5] = {kVermillion, kSkyBlue, kGreen, kOrange, kGrey};
    Marker(s.trkPhi[i], s.trkEta[i], role == 4 ? 24 : 20, SizeForPt(s.trkPt[i]), Color(col[role]), role == 4 ? 0.6 : 0.85);
    if (f & kGap) Marker(s.trkPhi[i], s.trkEta[i], 24, SizeForPt(s.trkPt[i]) + 0.5, Color(kBlack), 0.8);
  }
  if (haveJet) {
    DrawCone(s.leadEta, s.leadPhi, s.R, Color(kVermillion), 1, 3, etaRange);
    for (int c = 0; c < 2; c++) DrawCone(s.leadEta, s.perpPhi[c], s.R, Color(kBlue), 1, 3, etaRange);
  }
  // legend
  TPad* leg = SubPad(p.body, "leg", 0.47, 0.0, 0.64, 1.0);
  double y = 0.95;
  LegendText(y, "track role (N)", 0.07, 62); y -= 0.055;
  LegendRow(y, 20, 1.5, Color(kVermillion), Form("leading-jet constituent (%d)", nRole[0]), 0.058); y -= 0.05;
  LegendRow(y, 20, 1.5, Color(kSkyBlue), Form("perpendicular cone (%d)", nRole[1]), 0.058); y -= 0.05;
  LegendRow(y, 20, 1.5, Color(kGreen), Form("reference (%d)", nRole[2]), 0.058); y -= 0.05;
  LegendRow(y, 20, 1.5, Color(kOrange), Form("reference, not in refaw (%d)", nRole[3]), 0.058); y -= 0.05;
  LegendRow(y, 24, 1.5, Color(kGrey), Form("other (p_{T} or cone) (%d)", nRole[4]), 0.058); y -= 0.05;
  LegendRow(y, 24, 2.0, Color(kBlack), "ring: #eta-gap reference", 0.058); y -= 0.06;
  LegendLine(y, Color(kVermillion), 1, 3, "jet cone", 0.058); y -= 0.05;
  LegendLine(y, Color(kBlue), 1, 3, "perpendicular cones", 0.058); y -= 0.05;
  LegendLine(y, Color(kBlack), 2, 1, "|#eta - #eta_{jet}| = R + 0.5", 0.058); y -= 0.05;
  LegendBox(y, Color(kOrange), 0.3, "away side |#Delta#varphi - #pi| < 1", 0.058); y -= 0.07;
  LegendText(y, "reference: 0.2 #leq p_{T} < 3 GeV/#it{c}", 0.055); y -= 0.045;
  LegendText(y, "outside jet and perp. cones", 0.055); y -= 0.06;
  if (haveJet) {
    LegendText(y, Form("#rho = %.2f GeV/#it{c}", s.rho), 0.058); y -= 0.045;
    LegendText(y, Form("N_{ref} = %.0f, #LTp_{T}#GT_{ref} = %.3f", s.refN, s.refN > 0 ? s.refS1 / s.refN : 0.), 0.058); y -= 0.045;
  } else {
    LegendText(y, "no leading jet in this slot", 0.058);
  }
  // spectra
  const int nb = 25;
  TH1F* hN = new TH1F(Form("hN%d", gPageCounter), "", nb, s.probeEdges.data());
  TH1F* hB = new TH1F(Form("hB%d", gPageCounter), "", nb, s.probeEdges.data());
  TH1F* hS = new TH1F(Form("hS%d", gPageCounter), "", nb, s.probeEdges.data());
  TH1F* hF = new TH1F(Form("hF%d", gPageCounter), "", nb, s.probeEdges.data());
  TH1F* hFS = new TH1F(Form("hFS%d", gPageCounter), "", nb, s.probeEdges.data());
  double sumN = 0, sumB = 0, maxN = 0;
  for (int b = 0; b < nb; b++) {
    hN->SetBinContent(b + 1, s.nIn[b]); hB->SetBinContent(b + 1, s.nPerp[b]); hS->SetBinContent(b + 1, s.nSub[b]);
    hF->SetBinContent(b + 1, s.fRaw[b]); hFS->SetBinContent(b + 1, s.fSub[b]);
    sumN += s.nIn[b]; sumB += s.nPerp[b]; maxN = std::max(maxN, (double)s.nIn[b]);
  }
  double xMax = 1;
  for (int b = 0; b < nb; b++) if (s.nIn[b] > 0 || s.nPerp[b] > 0) xMax = s.probeEdges[b + 1];
  xMax = std::max(5.0, xMax * 1.5);
  for (TH1F* h : {hN, hB, hS, hF, hFS}) { h->SetDirectory(nullptr); h->SetStats(0); h->GetXaxis()->SetRangeUser(s.probeEdges[0] * 1.001, xMax); }
  auto style = [](TH1F* h, const char* col, int lineStyle, int width) {
    h->SetLineColor(Color(col)); h->SetLineStyle(lineStyle); h->SetLineWidth(width);
  };
  style(hN, kVermillion, 1, 3); style(hB, kBlue, 1, 3); style(hS, kBlack, 2, 2);
  style(hF, kVermillion, 1, 3); style(hFS, kBlack, 2, 2);
  TPad* top = SubPad(p.body, "spec", 0.64, 0.5, 1.0, 1.0);
  top->SetLogx(); top->SetLeftMargin(0.15); top->SetRightMargin(0.03); top->SetTopMargin(0.25); top->SetBottomMargin(0.2);
  hN->SetMinimum(std::min(-0.6, hS->GetMinimum() - 0.3)); hN->SetMaximum(maxN + 1);
  hN->GetXaxis()->SetTitle("probe p_{T} (GeV/#it{c})"); hN->GetYaxis()->SetTitle("tracks per bin");
  for (TH1F* h : {hN, hF}) {
    h->GetXaxis()->SetTitleSize(0.06); h->GetYaxis()->SetTitleSize(0.06); h->GetXaxis()->SetLabelSize(0.055);
    h->GetYaxis()->SetLabelSize(0.055); h->GetXaxis()->SetTitleOffset(1.3); h->GetYaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetMoreLogLabels();
  }
  hN->Draw("HIST"); hB->Draw("HIST same"); hS->Draw("HIST same");
  HLine(0.15, xMax, 0, Color(kGrey), 3, 1);
  TLegend* l1 = new TLegend(0.02, 0.77, 0.99, 0.99);
  l1->SetNColumns(2); l1->SetBorderSize(0); l1->SetFillStyle(0); l1->SetTextSize(0.055);
  l1->AddEntry(hN, Form("n_{b} in cone (N = %.0f)", sumN), "l");
  l1->AddEntry(hB, Form("b_{b} perp. cones (B = %.1f)", sumB), "l");
  l1->AddEntry(hS, "n_{b} - b_{b}", "l");
  l1->Draw();
  TPad* bot = SubPad(p.body, "frac", 0.64, 0.0, 1.0, 0.5);
  bot->SetLogx(); bot->SetLeftMargin(0.15); bot->SetRightMargin(0.03); bot->SetTopMargin(0.2); bot->SetBottomMargin(0.2);
  double fMax = std::max(hF->GetMaximum(), hFS->GetMaximum()), fMin = std::min(0.0, hFS->GetMinimum());
  hF->SetMinimum(fMin - 0.03); hF->SetMaximum(fMax * 1.15 + 0.01);
  hF->GetXaxis()->SetTitle("probe p_{T} (GeV/#it{c})"); hF->GetYaxis()->SetTitle("fraction per bin");
  hF->Draw("HIST"); hFS->Draw("HIST same");
  HLine(0.15, xMax, 0, Color(kGrey), 3, 1);
  TLegend* l2 = new TLegend(0.02, 0.82, 0.99, 0.99);
  l2->SetNColumns(2); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.055);
  l2->AddEntry(hF, "f_{b} = n_{b} / N", "l");
  l2->AddEntry(hFS, "f_{sub} = (n_{b} - b_{b}) / (N - B)", "l");
  l2->Draw();
  p.canvas->cd(); p.canvas->Print(pdf.c_str()); p.canvas->Print((png + "_c_background.png").c_str());
}

// ---------------------------------------------------------------------------
//  (d) constituent subtraction
// ---------------------------------------------------------------------------
void PageConstituentSubtraction(const SlotEntry& s, const std::string& head, const std::string& pdf,
                                const std::string& png) {
  Page p = NewPage(head + ", " + SlotText(s),
                   Form("(d) Event-wise constituent subtraction (#DeltaR_{max} = %.2f, #alpha = 0, ghost area %.2f, "
                        "#rho_{CS} = %.2f GeV/#it{c}): track p_{T} before/after, CS ghost grid near the jet, jet p_{T} definitions.",
                        s.csMaxDR, s.csGhostArea, s.rhoCS));
  // left: pT after vs before
  TPad* sc = SubPad(p.body, "scatter", 0.0, 0.0, 0.34, 1.0);
  sc->SetLogx(); sc->SetLogy(); sc->SetLeftMargin(0.16); sc->SetRightMargin(0.03); sc->SetBottomMargin(0.12); sc->SetTopMargin(0.13);
  const double removedY = 1.5e-3;
  TH2F* fr = new TH2F(Form("csframe%d", gPageCounter), "", 10, 0.12, 120, 10, 1e-3, 120);
  fr->SetDirectory(nullptr); fr->SetStats(0);
  fr->GetXaxis()->SetTitle("p_{T} before CS (GeV/#it{c})"); fr->GetYaxis()->SetTitle("p_{T} after CS (GeV/#it{c})");
  fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleOffset(1.4);
  fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetXaxis()->SetMoreLogLabels();
  fr->Draw("AXIS");
  TLine* diag = new TLine(0.12, 0.12, 120, 120); diag->SetLineColor(Color(kGrey)); diag->SetLineStyle(2); diag->Draw();
  HLine(0.12, 120, 2.2e-3, Color(kGrey), 3, 1);
  TLatex rl; rl.SetTextFont(42); rl.SetTextSize(0.04); rl.SetTextColor(Color("#555555"));
  rl.DrawLatex(8, 2.6e-3, "removed (p_{T}^{CS} = 0)");
  int nRemoved = 0; double ptRemoved = 0, ptBefore = 0, ptAfter = 0;
  for (size_t i = 0; i < s.trkPt.size(); i++) {
    int o = std::min(2, std::max(0, s.trkOrigin[i]));
    bool removed = s.trkPtCS[i] <= 0;
    ptBefore += s.trkPt[i]; ptAfter += s.trkPtCS[i];
    if (removed) { nRemoved++; ptRemoved += s.trkPt[i]; }
    bool inCSJet = s.trkFlags[i] & kCSJet;
    Marker(s.trkPt[i], removed ? removedY : s.trkPtCS[i], removed ? 5 : (inCSJet ? 20 : 24), removed ? 0.8 : 0.9,
           Color(kOriginColor[o]), 0.8);
  }
  TLegend* l = new TLegend(0.02, 0.875, 0.99, 0.995);
  l->SetNColumns(2); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.04);
  TMarker* mGen = new TMarker(0, 0, 24); mGen->SetMarkerColor(Color(kOriginColor[0]));
  TMarker* mUE = new TMarker(0, 0, 24); mUE->SetMarkerColor(Color(kOriginColor[1]));
  TMarker* mJet = new TMarker(0, 0, 20); mJet->SetMarkerColor(Color(kBlack));
  TMarker* mRem = new TMarker(0, 0, 5); mRem->SetMarkerColor(Color(kBlack));
  l->AddEntry(mGen, "generator", "p"); l->AddEntry(mUE, "added UE", "p");
  l->AddEntry(mJet, "filled: in CS leading jet", "p"); l->AddEntry(mRem, "removed", "p");
  l->Draw();
  // middle: zoom on the leading jet
  TPad* zm = SubPad(p.body, "zoom", 0.34, 0.0, 0.70, 1.0);
  zm->SetLeftMargin(0.17); zm->SetRightMargin(0.03); zm->SetBottomMargin(0.12); zm->SetTopMargin(0.13);
  const double dPhiRange = 1.2;
  EtaPhiFrame(s.trackEtaMax + 0.05, "#varphi - #varphi_{jet} (rad)", -dPhiRange, dPhiRange, "#eta", 1.5);
  HLine(-dPhiRange, dPhiRange, s.trackEtaMax, Color(kBlack), 1, 2); HLine(-dPhiRange, dPhiRange, -s.trackEtaMax, Color(kBlack), 1, 2);
  auto relPhi = [&](double phi) { return WrapPhi(phi - s.leadPhi); };
  for (size_t g = 0; g < s.csGhostEta.size(); g++) {
    double x = relPhi(s.csGhostPhi[g]);
    if (std::fabs(x) < dPhiRange) Marker(x, s.csGhostEta[g], 7, 0.5, Color("#BBBBBB"));
  }
  if (s.leadPt > 0) DrawCone(s.leadEta, 0.0, s.R, Color(kVermillion), 1, 2, s.trackEtaMax + 0.05);
  if (s.csJetPt > 0) DrawCone(s.csJetEta, relPhi(s.csJetPhi), s.R, Color(kBlue), 2, 2, s.trackEtaMax + 0.05);
  for (size_t i = 0; i < s.trkPt.size(); i++) {
    double x = relPhi(s.trkPhi[i]);
    if (std::fabs(x) > dPhiRange) continue;
    bool removed = s.trkPtCS[i] <= 0, inCSJet = s.trkFlags[i] & kCSJet;
    Marker(x, s.trkEta[i], 24, SizeForPt(s.trkPt[i]), Color(kBlack), 0.7);
    if (removed) Marker(x, s.trkEta[i], 5, 0.9, Color(kGrey));
    else Marker(x, s.trkEta[i], 20, SizeForPt(s.trkPtCS[i]), Color(inCSJet ? kVermillion : kSkyBlue), 0.85);
  }
  TLegend* lz = new TLegend(0.02, 0.875, 0.99, 0.995);
  lz->SetNColumns(3); lz->SetBorderSize(0); lz->SetFillStyle(0); lz->SetTextSize(0.035);
  TMarker* zb = new TMarker(0, 0, 24); zb->SetMarkerColor(Color(kBlack));
  TMarker* zj = new TMarker(0, 0, 20); zj->SetMarkerColor(Color(kVermillion));
  TMarker* zo = new TMarker(0, 0, 20); zo->SetMarkerColor(Color(kSkyBlue));
  TMarker* zr = new TMarker(0, 0, 5); zr->SetMarkerColor(Color(kGrey));
  TMarker* zg = new TMarker(0, 0, 7); zg->SetMarkerColor(Color("#999999"));
  TLine* zc = new TLine(); zc->SetLineColor(Color(kVermillion)); zc->SetLineWidth(2);
  TLine* zcs = new TLine(); zcs->SetLineColor(Color(kBlue)); zcs->SetLineWidth(2); zcs->SetLineStyle(2);
  lz->AddEntry(zb, "before CS", "p"); lz->AddEntry(zj, "after, in CS jet", "p"); lz->AddEntry(zo, "after, other", "p");
  lz->AddEntry(zr, "removed", "p"); lz->AddEntry(zg, "CS ghost", "p"); lz->AddEntry(zc, "raw jet R", "l");
  lz->AddEntry(zcs, "CS jet R", "l");
  lz->Draw();
  // right: table
  TPad* tb = SubPad(p.body, "table", 0.70, 0.0, 1.0, 1.0);
  TPaveText* pt = new TPaveText(0.04, 0.22, 0.98, 0.95, "NDC");
  pt->SetFillColor(0); pt->SetBorderSize(1); pt->SetTextFont(42); pt->SetTextSize(0.05); pt->SetTextAlign(12);
  pt->AddText(Form("#bf{leading jet, R = %.1f}", s.R));
  pt->AddText(Form("raw p_{T} = %.2f GeV/#it{c}", s.leadPt));
  pt->AddText(Form("area A = %.3f,  #rho = %.2f GeV/#it{c}", s.leadArea, s.rho));
  pt->AddText(Form("#rhoA-corrected p_{T} = %.2f GeV/#it{c}", s.leadPtCorr));
  pt->AddText(Form("CS jet p_{T} = %.2f GeV/#it{c}", s.csJetPt));
  pt->AddText(Form("truth (no added UE) p_{T} = %.2f GeV/#it{c}", s.jetPtTruth));
  pt->AddText(" ");
  pt->AddText("#bf{tracks in the acceptance}");
  pt->AddText(Form("#Sigma p_{T} before CS = %.1f GeV/#it{c}", ptBefore));
  pt->AddText(Form("#Sigma p_{T} after CS = %.1f GeV/#it{c}", ptAfter));
  pt->AddText(Form("removed: %d of %zu tracks (%.1f GeV/#it{c})", nRemoved, s.trkPt.size(), ptRemoved));
  pt->AddText(Form("subtracted #Sigma p_{T} = %.1f GeV/#it{c}", ptBefore - ptAfter));
  pt->AddText(Form("(#rho_{CS} #times acceptance = %.1f GeV/#it{c};", s.rhoCS * 2 * s.trackEtaMax * 2 * TMath::Pi()));
  pt->AddText(Form(" ghosts farther than #DeltaR_{max} keep their p_{T})"));
  pt->Draw();
  p.canvas->cd(); p.canvas->Print(pdf.c_str()); p.canvas->Print((png + "_d_cs.png").c_str());
}

// ---------------------------------------------------------------------------
//  (e) recoil thermal-momentum subtraction
// ---------------------------------------------------------------------------
void PageRecoil(const RecoilEntry& r, const RecoilEntry* other, const std::string& head, const std::string& pdf,
                const std::string& png) {
  Page p = NewPage(head + Form(", recoil thermal-momentum subtraction, #DeltaR_{max} = %.1f", r.dRmax),
                   "(e) Final-state particles and thermal scattering centres (status 3) before (open) and after (filled: "
                   "remainder) the subtraction; the smaller p_{T} of each pair (in order of #DeltaR) is removed from both.");
  const double yRange = 3.2;
  TPad* map = SubPad(p.body, "map", 0.0, 0.0, 0.56, 1.0);
  map->SetLeftMargin(0.09); map->SetRightMargin(0.02); map->SetBottomMargin(0.12); map->SetTopMargin(0.04);
  EtaPhiFrame(yRange, "#varphi (rad)", -TMath::Pi(), TMath::Pi(), "rapidity y");
  HLine(-TMath::Pi(), TMath::Pi(), 0.9, Color(kBlack), 3, 1); HLine(-TMath::Pi(), TMath::Pi(), -0.9, Color(kBlack), 3, 1);
  int nUnchanged = 0, nReduced = 0, nRemoved = 0;
  for (size_t i = 0; i < r.finPt.size(); i++) {
    if (std::fabs(r.finY[i]) > yRange) continue;
    double before = r.finPt[i], after = r.finPtAfter[i];
    if (after <= 0) { nRemoved++; Marker(r.finPhi[i], r.finY[i], 5, SizeForPt(before) + 0.3, Color(kPurple), 0.9); }
    else if (after != before) {
      nReduced++;
      Marker(r.finPhi[i], r.finY[i], 24, SizeForPt(before), Color(kGreen), 0.9);
      Marker(r.finPhi[i], r.finY[i], 20, SizeForPt(after), Color(kGreen), 0.9);
    } else { nUnchanged++; Marker(r.finPhi[i], r.finY[i], 20, SizeForPt(before), Color(kBlue), 0.55); }
  }
  int nThermal = 0, nResidual = 0;
  for (size_t i = 0; i < r.thPt.size(); i++) {
    if (std::fabs(r.thY[i]) > yRange) continue;
    nThermal++;
    Marker(r.thPhi[i], r.thY[i], 25, SizeForPt(r.thPt[i]), Color(kOrange), 0.95);
    if (r.thPtAfter[i] > 0) { nResidual++; Marker(r.thPhi[i], r.thY[i], 21, SizeForPt(r.thPtAfter[i]), Color(kOrange), 0.95); }
  }
  // spectra
  TPad* sp = SubPad(p.body, "spectra", 0.56, 0.0, 0.80, 1.0);
  sp->SetLogx(); sp->SetLogy(); sp->SetLeftMargin(0.2); sp->SetRightMargin(0.03); sp->SetBottomMargin(0.12); sp->SetTopMargin(0.17);
  double edges[31]; for (int i = 0; i <= 30; i++) edges[i] = 0.01 * std::pow(10000., i / 30.);
  TH1F* hThB = new TH1F(Form("thb%d", gPageCounter), "", 30, edges);
  TH1F* hThA = new TH1F(Form("tha%d", gPageCounter), "", 30, edges);
  TH1F* hFiB = new TH1F(Form("fib%d", gPageCounter), "", 30, edges);
  TH1F* hFiA = new TH1F(Form("fia%d", gPageCounter), "", 30, edges);
  for (size_t i = 0; i < r.thPt.size(); i++) { hThB->Fill(r.thPt[i]); if (r.thPtAfter[i] > 0) hThA->Fill(r.thPtAfter[i]); }
  for (size_t i = 0; i < r.finPt.size(); i++) { hFiB->Fill(r.finPt[i]); if (r.finPtAfter[i] > 0) hFiA->Fill(r.finPtAfter[i]); }
  auto style = [](TH1F* h, const char* col, int ls) {
    h->SetDirectory(nullptr); h->SetStats(0); h->SetLineColor(Color(col)); h->SetLineStyle(ls); h->SetLineWidth(3);
  };
  style(hFiB, kBlue, 1); style(hFiA, kGreen, 2); style(hThB, kOrange, 1); style(hThA, kVermillion, 2);
  double ymax = std::max(hFiB->GetMaximum(), hThB->GetMaximum());
  hFiB->SetMaximum(ymax * 3); hFiB->SetMinimum(0.5);
  hFiB->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})"); hFiB->GetYaxis()->SetTitle("particles per bin");
  hFiB->GetXaxis()->SetTitleSize(0.05); hFiB->GetYaxis()->SetTitleSize(0.05); hFiB->GetYaxis()->SetTitleOffset(1.7);
  hFiB->GetXaxis()->SetLabelSize(0.045); hFiB->GetYaxis()->SetLabelSize(0.045);
  hFiB->Draw("HIST"); hFiA->Draw("HIST same"); hThB->Draw("HIST same"); hThA->Draw("HIST same");
  TLegend* l = new TLegend(0.02, 0.83, 0.99, 0.995);
  l->SetNColumns(2); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.042);
  l->AddEntry(hFiB, "final, before", "l"); l->AddEntry(hFiA, "final, after", "l");
  l->AddEntry(hThB, "thermal, before", "l"); l->AddEntry(hThA, "thermal remainder", "l");
  l->Draw();
  // legend + table
  TPad* tb = SubPad(p.body, "table", 0.80, 0.0, 1.0, 1.0);
  double y = 0.95;
  LegendText(y, Form("map (|y| < %.1f)", yRange), 0.06, 62); y -= 0.05;
  LegendRow(y, 20, 1.4, Color(kBlue), Form("final, unchanged (%d)", nUnchanged), 0.05); y -= 0.045;
  LegendRow(y, 20, 1.4, Color(kGreen), Form("final, reduced (%d)", nReduced), 0.05); y -= 0.045;
  LegendRow(y, 5, 1.4, Color(kPurple), Form("final, removed (%d)", nRemoved), 0.05); y -= 0.045;
  LegendRow(y, 25, 1.4, Color(kOrange), Form("thermal, before (%d)", nThermal), 0.05); y -= 0.045;
  LegendRow(y, 21, 1.4, Color(kOrange), Form("thermal remainder (%d)", nResidual), 0.05); y -= 0.045;
  LegendLine(y, Color(kBlack), 3, 1, "|y| = 0.9", 0.05); y -= 0.06;
  auto bookkeeping = [&](const RecoilEntry& e, double y0) {
    double fb = 0, fa = 0, tbf = 0, ta = 0;
    for (size_t i = 0; i < e.finPt.size(); i++) { fb += e.finPt[i]; fa += e.finPtAfter[i]; }
    for (size_t i = 0; i < e.thPt.size(); i++) { tbf += e.thPt[i]; ta += e.thPtAfter[i]; }
    LegendText(y0, Form("#DeltaR_{max} = %.1f (%lld pairs)", e.dRmax, e.nPairs), 0.05, 62); y0 -= 0.042;
    LegendText(y0, Form("#Sigma p_{T} final: %.1f #rightarrow %.1f", fb, fa), 0.048); y0 -= 0.04;
    LegendText(y0, Form("#Sigma p_{T} thermal: %.1f #rightarrow %.1f", tbf, ta), 0.048); y0 -= 0.04;
    LegendText(y0, Form("removed from each: %.1f / %.1f", fb - fa, tbf - ta), 0.048); y0 -= 0.055;
    return y0;
  };
  LegendText(y, "#bf{p_{T} bookkeeping (GeV/#it{c})}", 0.05); y -= 0.05;
  y = bookkeeping(r, y);
  if (other) y = bookkeeping(*other, y);
  y = DrawPtKey(y - 0.01, 0.045);
  p.canvas->cd(); p.canvas->Print(pdf.c_str()); p.canvas->Print((png + "_e_recoil.png").c_str());
}

std::vector<std::string> SplitList(const std::string& s) {
  std::vector<std::string> out; std::stringstream ss(s); std::string x;
  while (std::getline(ss, x, ',')) if (!x.empty()) out.push_back(x);
  return out;
}

}  // namespace

void PlotEventDisplay(const char* displayFile = "display.root", const char* outDir = "plots_display",
                      const char* sampleLabel = "", const char* fileTag = "display", int event = -1,
                      const char* variantList = "", bool withR08 = true) {
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.5);
  gSystem->mkdir(outDir, true);
  TFile f(displayFile);
  if (f.IsZombie()) { printf("cannot open %s\n", displayFile); return; }
  std::vector<SlotEntry> slots = ReadSlots(f);
  std::vector<RecoilEntry> recoils = ReadRecoils(f);
  if (slots.empty()) { printf("no display entries in %s\n", displayFile); return; }

  // variants in file order, or the requested ones
  std::vector<std::string> variants = SplitList(variantList);
  if (variants.empty())
    for (auto& s : slots) if (std::find(variants.begin(), variants.end(), s.variant) == variants.end()) variants.push_back(s.variant);
  // event: requested, or highest R = 0.4 leading-jet pT in the first variant
  if (event < 0) {
    double best = -1;
    for (auto& s : slots)
      if (s.variant == variants[0] && s.iacc == 0 && s.leadPt > best) { best = s.leadPt; event = s.event; }
  }
  printf("event %d, variants:", event);
  for (auto& v : variants) printf(" %s", v.c_str());
  printf("\n");

  std::string base = Form("%s/%s_ev%d", outDir, fileTag, event), pdf = base + ".pdf";
  TCanvas* opener = new TCanvas("opener", "", 1600, 900);
  opener->Print((pdf + "[").c_str());
  for (const std::string& v : variants) {
    const SlotEntry *s04 = nullptr, *s08 = nullptr;
    for (auto& s : slots) if (s.event == event && s.variant == v) (s.iacc == 0 ? s04 : s08) = &s;
    if (!s04) { printf("  no R = 0.4 entry for %s\n", v.c_str()); continue; }
    std::string head = Form("%s%sevent %d, %s", sampleLabel, sampleLabel[0] ? ", " : "", event, VariantDescription(v).c_str());
    std::string png = base + "_" + v;
    PageParticles(*s04, head, pdf, png);
    PageJets(*s04, head, pdf, png, "b_jets");
    if (withR08 && s08) PageJets(*s08, head, pdf, png, "b2_jets_R08");
    PageBackground(*s04, head, pdf, png);
    PageConstituentSubtraction(*s04, head, pdf, png);
  }
  const RecoilEntry *r10 = nullptr, *r05 = nullptr;
  for (auto& r : recoils) if (r.event == event) (r.dRmax > 0.75 ? r10 : r05) = &r;
  if (r10) PageRecoil(*r10, r05, Form("%s%sevent %d", sampleLabel, sampleLabel[0] ? ", " : "", event), pdf, base);
  opener->Print((pdf + "]").c_str());
  printf("wrote %s and the PNGs %s_*.png\n", pdf.c_str(), base.c_str());
}
