// ============================================================================
//  Phase1_corr.C -- stage "corr": second moments vs fine Nch, per R and jet
//  class: <probe_b [pT]_ref> for the three probes and the two-particle C2 of
//  [pT]_ref. Needs the MERGED mean of the same sample and configuration.
//
//  Mean input (mu(Nch) for C2 and the closure check), first found of:
//    $V0_MEANFILE                              (batch: merged mean copied locally)
//    OUTPUT/<REFN[IREF]>/<Tag>_mean.root       (per-reference local layout)
//    OUTPUT/<Tag>_mean.root                    (flat local layout)
//
//  Event gate: identical to Phase1_mean.C (N_ref >= 2, S2 present, [pT]_ref >= 0).
//
//  Two-particle C2 estimator (dynamical variance of [pT]_ref, self-pairs removed):
//      C2 = < (pT_i - mu)(pT_j - mu) >_{i != j}
//         = [ (S1 - N mu)^2 - (S2 - 2 mu S1 + N mu^2) ] / [ N (N - 1) ]
//  with N = N_ref, S1 = sum pT, S2 = sum pT^2 of the reference particles and
//  mu = mu(Nch) read from the nominal mean profile at this event's Nch bin
//  (if the profile is missing, the event's own [pT]_ref is used). Two versions:
//    pC2_Nch    mu from pRef_Nch      (used with f)
//    pC2Sub_Nch mu from pRefFsub_Nch  (used with f_sub and f_cs in Phase2)
//  Both are filled for every gated event of the class.
//
//  Products (same event subsets as the mean stage):
//    f     : f_b [pT]_ref               always
//    f_sub : f_sub,b [pT]_ref           if sum (n - p) > 0
//    f_cs  : f_cs,b [pT]_ref            if sum cs > 0
//
//  Bootstrap: as in the mean stage (same global-event seed => replica k is
//  coherent across mean and corr).
//
//  Output <Tag>_corr_<IJob>.root (k = 0 has no "_boot<k>"):
//    pNR_NchPt<VARSUF>[_boot<k>]_R<i>_<SELN>   TProfile2D <probe_b [pT]_ref> vs (Nch, pT)
//    pC2_Nch[_boot<k>]_R<i>_<SELN>             TProfile   <c2> vs Nch (mu from pRef_Nch)
//    pC2Sub_Nch[_boot<k>]_R<i>_<SELN>          TProfile   <c2> vs Nch (mu from pRefFsub_Nch)
//    hCov_NchPt_R<i>_<SELN>  TH2D  nominal Cov(f_b,[pT]_ref) = <f_b ref> - <f_b> mu
//    Nf_NchPt_R<i>_<SELN>    copy of the nominal <f_b> (mean)
//    mu_Nch_R<i>_<SELN>      copy of the nominal mu(Nch) (mean)
//  Closure check printed at the end: max over (R, class, Nch bin) of
//  |sum_b Cov| for f, f_sub and f_cs; must be ~0 (threshold 1e-4).
// ============================================================================
void Phase1::ReadHists_corr(){
  TString meanPath = Form("OUTPUT/%s/%s_mean.root", REFN[IREF], Tag.c_str());
  { ifstream test(meanPath.Data()); if (!test.good()) meanPath = Form("OUTPUT/%s_mean.root", Tag.c_str()); }
  if (getenv("V0_MEANFILE")) meanPath = getenv("V0_MEANFILE");
  fMean = TFile::Open(meanPath, "read");
  if (!fMean || fMean->IsZombie()){
    cout << "  ERROR: cannot open " << meanPath << " (hadd mean first)\n";
    for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++){ muRefPt[iR][c] = 0; muRefPtFsub[iR][c] = 0; }
    return;
  }
  for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++){
    muRefPt[iR][c]     = (TProfile*)fMean->Get(Form("pRef_Nch_R%d_%s",     iR, SELN[c]));
    muRefPtFsub[iR][c] = (TProfile*)fMean->Get(Form("pRefFsub_Nch_R%d_%s", iR, SELN[c]));
  }
  cout << "  corr: read mu(Nch) from " << meanPath << endl;
}

void Phase1::InitHists_corr(){
  for (int k = 0; k <= NBOOT; k++) for (int iv = 0; iv < NVAR; iv++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++)
    pFracRefPt[k][iv][iR][c] = new TProfile2D(BootName(Form("pNR_NchPt%s", VARSUF[iv]), k, iR, c),
                                              ";N_{ch};p_{T}^{probe};<probe*ref>", NCHBINS, 0, NCHMAX, nB, ptEdges);
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++)
    pC2[k][iR][c] = new TProfile(BootName("pC2_Nch", k, iR, c), ";N_{ch};<c2>", NCHBINS, 0, NCHMAX);
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++)
    pC2Fsub[k][iR][c] = new TProfile(BootName("pC2Sub_Nch", k, iR, c), ";N_{ch};<c2>", NCHBINS, 0, NCHMAX);
}

// two-particle C2 of [pT]_ref for N reference particles with sums S1, S2 about mu
static inline double TwoParticleC2(double N, double S1, double S2, double mu){
  double P1 = S1 - N*mu, P2 = S2 - 2.0*mu*S1 + N*mu*mu;
  return (P1*P1 - P2)/(N*(N - 1.0));
}

void Phase1::Fill_corr(Event* evt, Long64_t globalEvt){
  double nch = (evt->nch >= 0) ? evt->nch : 0.0;
  double bootWeight[NBOOT+1]; bootWeight[0] = evt->weight;
  for (int k = 1; k <= NBOOT; k++) bootWeight[k] = evt->weight*PoissonBootWeight(globalEvt, k);

  for (int iR = 0; iR < nR; iR++){
    double refMeanPt = evt->refMeanPt[iR];
    bool passGate = (refMeanPt >= 0 && evt->refCount[iR] >= 2 && evt->refSumPt2[iR] >= 0);   // == Phase1_mean.C
    if (!passGate) continue;
    double N = evt->refCount[iR], S1 = evt->refSumPt[iR], S2 = evt->refSumPt2[iR];
    for (int c = 0; c < NSEL; c++){
      if (!evt->inJetClass[iR][c]) continue;
      TProfile* mu  = muRefPt[iR][c];
      TProfile* muS = muRefPtFsub[iR][c];
      double muF    = mu  ? mu ->GetBinContent(mu ->FindBin(nch)) : refMeanPt;
      double muFsub = muS ? muS->GetBinContent(muS->FindBin(nch)) : refMeanPt;
      double c2F    = TwoParticleC2(N, S1, S2, muF);
      double c2Fsub = TwoParticleC2(N, S1, S2, muFsub);

      const double* n  = evt->nInCone[iR][c];
      const double* p  = evt->nPerpCone[iR][c];
      const double* cs = evt->nCsJet[iR][c];
      double sumRaw = 0.0, sumSub = 0.0, sumCs = 0.0;
      for (int ip = 0; ip < nB; ip++){ sumRaw += n[ip]; sumSub += (n[ip] - p[ip]); sumCs += cs[ip]; }
      if (sumRaw <= 0) continue;
      bool haveFsub = (sumSub > 0), haveFcs = (sumCs > 0);
      for (int k = 0; k <= NBOOT; k++){
        double wk = bootWeight[k]; if (wk <= 0) continue;
        pC2    [k][iR][c]->Fill(nch, c2F,    wk);
        pC2Fsub[k][iR][c]->Fill(nch, c2Fsub, wk);
        for (int ip = 0; ip < nB; ip++){
          double pt = ptCenter[ip];
          pFracRefPt[k][IVAR_F][iR][c]->Fill(nch, pt, (n[ip]/sumRaw)*refMeanPt, wk);
          if (haveFsub) pFracRefPt[k][IVAR_FSUB][iR][c]->Fill(nch, pt, ((n[ip] - p[ip])/sumSub)*refMeanPt, wk);
          if (haveFcs)  pFracRefPt[k][IVAR_FCS ][iR][c]->Fill(nch, pt, (cs[ip]/sumCs)*refMeanPt, wk);
        }
      }
    }
  }
}

void Phase1::Run_corr(){
  Long64_t nEntries = fChain->GetEntries(), first, last; jobRange(nEntries, IJob, NEv, first, last);
  cout << "  corr: events [" << first << "," << last << ") of " << nEntries << ", " << NBOOT
       << " bootstrap replicas/event, " << NVAR << " probe variants" << endl;
  for (Long64_t e = first; e < last; e++){
    fChain->GetEntry(e); Tool::progress(e-first, last-first);
    Event ev; Fill_Evt(&ev); Fill_corr(&ev, e + evtOffset());
  }
  cout << "\r  corr: done (" << last-first << ")        " << endl;
}

// max over (R, class, Nch bin with entries) of |sum_b (<probe_b ref> - <probe_b> mu)|
static double WorstCovSum(TFile* fMean, TProfile2D* (*probeRef)[NSEL], const char* probeStem, const char* muStem,
                          int nR, TString& where, bool& found){
  double worst = 0; found = false;
  for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++){
    TProfile2D* probe = (TProfile2D*)(fMean ? fMean->Get(Form("%s_R%d_%s", probeStem, iR, SELN[c])) : 0);
    TProfile*   mu    = (TProfile*)  (fMean ? fMean->Get(Form("%s_R%d_%s", muStem,    iR, SELN[c])) : 0);
    if (!probe || !mu) continue;
    found = true;
    TProfile2D* pr = probeRef[iR][c];
    int nx = pr->GetNbinsX(), ny = pr->GetNbinsY();
    for (int jx = 1; jx <= nx; jx++){
      if (mu->GetBinEntries(jx) <= 0) continue;
      double m = mu->GetBinContent(jx), s = 0;
      for (int ip = 1; ip <= ny; ip++) s += pr->GetBinContent(jx, ip) - probe->GetBinContent(jx, ip)*m;
      if (TMath::Abs(s) > TMath::Abs(worst)){ worst = s; where = Form("R%d_%s Nchbin%d", iR, SELN[c], jx); }
    }
  }
  return worst;
}

void Phase1::Final_corr(){
  fOut->cd();
  for (int k = 0; k <= NBOOT; k++) for (int iv = 0; iv < NVAR; iv++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++)
    pFracRefPt[k][iv][iR][c]->Write();
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++) pC2[k][iR][c]->Write();
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++) pC2Fsub[k][iR][c]->Write();

  // nominal covariance of f and its ingredients (Phase2 does not read these)
  bool warned = false;
  for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++){
    TProfile2D* fracMean = (TProfile2D*)(fMean ? fMean->Get(Form("pN_NchPt_R%d_%s", iR, SELN[c])) : 0);
    TProfile*   mu       = muRefPt[iR][c];
    if (!fracMean || !mu){
      if (!warned) cout << "  corr CHECK: mean pN/pRef missing -> cannot store cov or self-check\n";
      warned = true; continue;
    }
    TProfile2D* fracRef = pFracRefPt[0][IVAR_F][iR][c];
    int nx = fracRef->GetNbinsX(), ny = fracRef->GetNbinsY();
    TH2D* hCov = new TH2D(Form("hCov_NchPt_R%d_%s", iR, SELN[c]),
                          ";N_{ch};p_{T}^{probe};cov(f,[p_{T}]_{ref})", nx, 0, NCHMAX, ny, ptEdges);
    for (int jx = 1; jx <= nx; jx++){
      double m = mu->GetBinContent(jx);
      for (int ip = 1; ip <= ny; ip++)
        hCov->SetBinContent(jx, ip, fracRef->GetBinContent(jx, ip) - fracMean->GetBinContent(jx, ip)*m);
    }
    hCov->Write();
    fracMean->Clone(Form("Nf_NchPt_R%d_%s", iR, SELN[c]))->Write();
    mu->Clone(Form("mu_Nch_R%d_%s", iR, SELN[c]))->Write();
  }

  // closure check sum_b Cov(probe_b,[pT]_ref) = 0 for each probe with its matched mean
  TString whereF, whereS, whereC; bool haveF, haveS, haveC;
  double worstF = WorstCovSum(fMean, pFracRefPt[0][IVAR_F],    "pN_NchPt",      "pRef_Nch",     nR, whereF, haveF);
  double worstS = WorstCovSum(fMean, pFracRefPt[0][IVAR_FSUB], "pN_NchPt_fsub", "pRefFsub_Nch", nR, whereS, haveS);
  double worstC = WorstCovSum(fMean, pFracRefPt[0][IVAR_FCS],  "pN_NchPt_fcs",  "pRefFcs_Nch",  nR, whereC, haveC);
  if (TMath::Abs(worstC) > 1e-4) cout << "  corr CHECK: FAIL  f_cs max |Sum_pT cov| = " << worstC << "  (" << whereC << ")\n";
  else                           cout << "  corr CHECK: OK  f_cs max |Sum_pT cov| = " << worstC << "\n";
  if (TMath::Abs(worstF) > 1e-4)
    cout << "  corr CHECK: FAIL  max |Sum_pT cov| = " << worstF << "  (" << whereF << ")\n"
         << "             mean and corr are INCONSISTENT (stale mean? gate mismatch?): regenerate\n"
         << "             mean.root with the current config, then rerun corr + Phase2.\n";
  else
    cout << "  corr CHECK: OK  max |Sum_pT cov| = " << worstF << "  (~0: mean/corr consistent, f variant)\n";
  if (!haveS)
    cout << "  corr CHECK: f_sub -> pRefFsub_Nch absent in mean; regenerate mean.\n";
  else if (TMath::Abs(worstS) > 1e-4)
    cout << "  corr CHECK: FAIL  f_sub max |Sum_pT cov| = " << worstS << "  (" << whereS << ")\n";
  else
    cout << "  corr CHECK: OK  f_sub max |Sum_pT cov| = " << worstS << "  (~0: matched mean, f_sub exact)\n";

  fOut->Close(); if (fMean) fMean->Close();
}
