// ============================================================================
//  JetUEResults.h -- reading Phase2 results and forming differences.
//
//  Phase2 writes one file per comparison triple,
//    <resultsDir>/<cfg>/acc<iacc>/v0pt_<triple>.root ,
//  where cfg = <reference>_c<lo>-<hi>|call_jt<n> (see pipeline/cfg.sh) and
//  iacc = 0 (|eta_trk| < 0.9) or 1 (|eta_trk| < 2.5). Every sample tag in the
//  triple (e.g. OO_RC_cs10_embFlow) has, per radius index iR and jet class:
//    v0<probe>_<tag>_c1_R<iR>_<jetClass>       TH1D, v0 vs probe pT (coarse
//                                              grid), bootstrap error
//    v0boot<probe>_<tag>_c1_R<iR>_<jetClass>   TH2D, x = probe-pT bin,
//                                              y = replica k (bin k+1; k=0 is
//                                              the nominal value)
//  <probe> is "" (raw fraction f), "_fsub" (perpendicular-cone subtracted) or
//  "_fcs" (constituent subtracted). "c1" is the inclusive event class.
//
//  Replica k of every sample uses the same per-event Poisson weights (seeded by
//  the global event index). For two samples built from the SAME events (e.g.
//  JEWEL vacuum with and without embedded UE) the replica-wise difference
//  therefore carries their correlation; for independent samples it reduces to
//  the usual quadrature error. All differences below use replica-wise errors.
// ============================================================================
#ifndef JetUEResults_h
#define JetUEResults_h
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace jetue {

// radius list and the acceptance index used for each radius in the paper:
// R <= 0.4 with |eta_trk| < 0.9 (acc0), R = 0.8 with |eta_trk| < 2.5 (acc1)
static const int    kNumRadii = 4;
static const double kRadius[kNumRadii]          = {0.1, 0.2, 0.4, 0.8};
static const char*  kRadiusLabel[kNumRadii]     = {"0.1", "0.2", "0.4", "0.8"};
static const int    kNominalAcceptance[kNumRadii] = {0, 0, 0, 1};

// a v0(pT) curve: central values with bootstrap errors, and its replicas
struct V0Curve {
  TH1D* central = nullptr;   // central value and error per probe-pT bin
  TH2D* replicas = nullptr;  // [pT bin, replica]; replica 0 = nominal
  bool ok() const { return central != nullptr; }
};

// Same trimmed spread as Phase2 (calErr): drop non-finite or exactly-zero
// replicas (zero = empty in that replica), take the population standard
// deviation, drop > 3 sigma outliers and recompute. Returns -1 for < 2 values.
inline double TrimmedSpread(const std::vector<double>& values) {
  int n = 0; double mean = 0;
  for (double v : values) if (TMath::Finite(v) && v != 0.0) { mean += v; n++; }
  if (n < 2) return -1;
  mean /= n;
  double var = 0;
  for (double v : values) if (TMath::Finite(v) && v != 0.0) var += (v - mean) * (v - mean);
  double sd = std::sqrt(var / n);
  if (sd <= 0) return 0;
  int n2 = 0; double mean2 = 0;
  for (double v : values) if (TMath::Finite(v) && v != 0.0 && std::fabs(v - mean) < 3 * sd) { mean2 += v; n2++; }
  if (n2 < 2) return sd;
  mean2 /= n2;
  double var2 = 0;
  for (double v : values) if (TMath::Finite(v) && v != 0.0 && std::fabs(v - mean) < 3 * sd) var2 += (v - mean2) * (v - mean2);
  return std::sqrt(var2 / n2);
}

// All Phase2 result files of one configuration and acceptance. A sample tag
// may appear in several triples; the first file that contains it is used.
class ResultSet {
 public:
  ResultSet(const std::string& resultsDir, const std::string& cfg, int acceptance) {
    TString dir = Form("%s/%s/acc%d", resultsDir.c_str(), cfg.c_str(), acceptance);
    void* d = gSystem->OpenDirectory(dir);
    if (!d) return;
    while (const char* entry = gSystem->GetDirEntry(d)) {
      TString name(entry);
      if (!name.BeginsWith("v0pt_") || !name.EndsWith(".root")) continue;
      TFile* f = TFile::Open(dir + "/" + name);
      if (f && !f->IsZombie()) files_.push_back(f);
    }
    gSystem->FreeDirectory(d);
  }
  bool empty() const { return files_.empty(); }

  // probe: "" (f), "_fsub", "_fcs"; radiusIndex 0..3; jetClass e.g. "jet40_60"
  V0Curve Get(const std::string& tag, int radiusIndex, const std::string& jetClass,
              const std::string& probe = "_fsub") const {
    TString suffix = Form("%s_%s_c1_R%d_%s", probe.c_str(), tag.c_str(), radiusIndex, jetClass.c_str());
    V0Curve c;
    for (TFile* f : files_) {
      TH1D* h = (TH1D*)f->Get("v0" + suffix);
      if (!h) continue;
      c.central = h;
      c.replicas = (TH2D*)f->Get("v0boot" + suffix);
      break;
    }
    return c;
  }

 private:
  std::vector<TFile*> files_;
};

// a - b per probe-pT bin. Error: replica-wise spread when both have replicas,
// otherwise the quadrature sum of the two bootstrap errors.
inline V0Curve Difference(const V0Curve& a, const V0Curve& b, const char* name) {
  V0Curve d;
  if (!a.ok() || !b.ok()) return d;
  d.central = (TH1D*)a.central->Clone(name);
  d.central->SetDirectory(nullptr);
  const int nBins = a.central->GetNbinsX();
  const bool haveReplicas = a.replicas && b.replicas;
  if (haveReplicas) {
    const int nRep = std::min(a.replicas->GetNbinsY(), b.replicas->GetNbinsY());
    d.replicas = new TH2D(Form("%s_replicas", name), "", nBins, 0, nBins, nRep, -0.5, nRep - 0.5);
    d.replicas->SetDirectory(nullptr);
    for (int q = 1; q <= nBins; q++)
      for (int k = 1; k <= nRep; k++) {
        double x = a.replicas->GetBinContent(q, k), y = b.replicas->GetBinContent(q, k);
        if (x != 0.0 && y != 0.0) d.replicas->SetBinContent(q, k, x - y);
      }
  }
  for (int q = 1; q <= nBins; q++) {
    d.central->SetBinContent(q, a.central->GetBinContent(q) - b.central->GetBinContent(q));
    double err = -1;
    if (haveReplicas) {
      std::vector<double> diffs;
      for (int k = 2; k <= d.replicas->GetNbinsY(); k++) diffs.push_back(d.replicas->GetBinContent(q, k));
      err = TrimmedSpread(diffs);
    }
    if (err < 0) err = std::hypot(a.central->GetBinError(q), b.central->GetBinError(q));
    d.central->SetBinError(q, err);
  }
  return d;
}

// Range summary over probe pT in [ptLow, ptHigh]: inverse-variance weighted
// mean of the bins fully inside the range. The weights are fixed from the
// central errors and applied to every replica, so the summary error is the
// replica spread of the weighted mean (correlations between bins included).
struct Summary { double value = 0, error = -1; bool ok() const { return error >= 0; } };

inline Summary RangeSummary(const V0Curve& c, double ptLow, double ptHigh) {
  Summary s;
  if (!c.ok()) return s;
  std::vector<int> bins; std::vector<double> weights;
  double sumW = 0, sumWV = 0;
  for (int q = 1; q <= c.central->GetNbinsX(); q++) {
    double lo = c.central->GetBinLowEdge(q), hi = lo + c.central->GetBinWidth(q);
    double err = c.central->GetBinError(q);
    if (lo < ptLow - 1e-6 || hi > ptHigh + 1e-6 || err <= 0) continue;
    double w = 1.0 / (err * err);
    bins.push_back(q); weights.push_back(w);
    sumW += w; sumWV += w * c.central->GetBinContent(q);
  }
  if (sumW <= 0) return s;
  s.value = sumWV / sumW;
  if (c.replicas) {
    std::vector<double> perReplica;
    for (int k = 2; k <= c.replicas->GetNbinsY(); k++) {
      double num = 0; bool complete = true;
      for (size_t i = 0; i < bins.size(); i++) {
        double v = c.replicas->GetBinContent(bins[i], k);
        if (v == 0.0) { complete = false; break; }
        num += weights[i] * v;
      }
      if (complete) perReplica.push_back(num / sumW);
    }
    s.error = TrimmedSpread(perReplica);
  }
  if (s.error < 0) s.error = 1.0 / std::sqrt(sumW);   // independent-bin fallback
  return s;
}

}  // namespace jetue
#endif
