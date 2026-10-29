// ============================================================================
//  Phase1_mean.C -- stage "mean": first moments vs fine Nch, per R and jet
//  class, for the three probe definitions and their matched reference means.
//
//  Event gate (identical in Phase1_corr.C; change both or neither):
//      [pT]_ref >= 0  AND  N_ref >= 2  AND  S2 present (>= 0)
//  N_ref >= 2 is needed by the two-particle C2 (division by N(N-1)). Because
//  mean and corr average over the identical event set, sum_b Cov(f_b,[pT]_ref)
//  = 0 holds exactly (sum_b f_b = 1); a mismatched gate or a stale mean shows up
//  as a positive pT-independent offset of the covariance.
//
//  Per event in jet class c (sums over all probe bins b):
//    f      : f_b     = n_b / sum n                 always (sum n > 0 by the class)
//    f_sub  : f_sub,b = (n_b - p_b) / sum (n - p)    only if sum (n - p) > 0
//    f_cs   : f_cs,b  = cs_b / sum cs                only if sum cs > 0
//  The reference mean <[pT]_ref> is filled once per event set:
//    pRef_Nch     all gated events of the class  (used with f)
//    pRefFsub_Nch the subset with sum (n - p) > 0 (used with f_sub)
//    pRefFcs_Nch  the subset with sum cs > 0      (used with f_cs)
//  so each probe has a reference mean over exactly its own events.
//
//  Bootstrap: every profile is filled for k = 0 (weight w) and k = 1..NBOOT
//  (weight w * Poisson(1), seeded by the global event index, Rule.h). The same
//  seed is used in the corr stage, so replica k is coherent across both files.
//  Replicas with weight 0 (~37% per event) are skipped.
//
//  Nch axis: events without a valid Nch (< 0) are filled at Nch = 0, so the
//  inclusive combination (all Nch bins) stays complete.
//
//  Output <Tag>_mean_<IJob>.root (k = 0 has no "_boot<k>"):
//    pN_NchPt<VARSUF>[_boot<k>]_R<i>_<SELN>   TProfile2D <probe_b> vs (Nch, probe pT)
//    pRef_Nch[_boot<k>]_R<i>_<SELN>           TProfile   <[pT]_ref> vs Nch (f events)
//    pRefFsub_Nch[_boot<k>]_R<i>_<SELN>       TProfile   same over f_sub events
//    pRefFcs_Nch[_boot<k>]_R<i>_<SELN>        TProfile   same over f_cs events
//    REFWIN                                   TNamed     reference stamp (Rule.h)
// ============================================================================
static TString BootName(const char* stem, int k, int iR, int c){
  return (k == 0) ? Form("%s_R%d_%s", stem, iR, SELN[c])
                  : Form("%s_boot%d_R%d_%s", stem, k, iR, SELN[c]);
}

void Phase1::InitHists_mean(){
  for (int k = 0; k <= NBOOT; k++) for (int iv = 0; iv < NVAR; iv++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++)
    pFrac[k][iv][iR][c] = new TProfile2D(BootName(Form("pN_NchPt%s", VARSUF[iv]), k, iR, c),
                                         ";N_{ch};p_{T}^{probe};<probe>", NCHBINS, 0, NCHMAX, nB, ptEdges);
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++)
    pRefPt[k][iR][c] = new TProfile(BootName("pRef_Nch", k, iR, c), ";N_{ch};<[p_{T}]_{ref}>", NCHBINS, 0, NCHMAX);
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++)
    pRefPtFsub[k][iR][c] = new TProfile(BootName("pRefFsub_Nch", k, iR, c), ";N_{ch};<[p_{T}]_{ref}> (f_sub set E')", NCHBINS, 0, NCHMAX);
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++)
    pRefPtFcs[k][iR][c] = new TProfile(BootName("pRefFcs_Nch", k, iR, c), ";N_{ch};<[p_{T}]_{ref}> (f_cs set E'')", NCHBINS, 0, NCHMAX);
}

void Phase1::Fill_mean(Event* evt, Long64_t globalEvt){
  double nch = (evt->nch >= 0) ? evt->nch : 0.0;
  double bootWeight[NBOOT+1]; bootWeight[0] = evt->weight;
  for (int k = 1; k <= NBOOT; k++) bootWeight[k] = evt->weight*PoissonBootWeight(globalEvt, k);

  for (int iR = 0; iR < nR; iR++){
    double refMeanPt = evt->refMeanPt[iR];
    bool passGate = (refMeanPt >= 0 && evt->refCount[iR] >= 2 && evt->refSumPt2[iR] >= 0);
    if (!passGate) continue;
    for (int c = 0; c < NSEL; c++){
      if (!evt->inJetClass[iR][c]) continue;
      const double* n  = evt->nInCone[iR][c];
      const double* p  = evt->nPerpCone[iR][c];
      const double* cs = evt->nCsJet[iR][c];
      double sumRaw = 0.0, sumSub = 0.0, sumCs = 0.0;   // normalisations of f, f_sub, f_cs
      for (int ip = 0; ip < nB; ip++){ sumRaw += n[ip]; sumSub += (n[ip] - p[ip]); sumCs += cs[ip]; }
      if (sumRaw <= 0) continue;                         // already excluded by inJetClass
      bool haveFsub = (sumSub > 0);
      bool haveFcs  = (sumCs  > 0);
      for (int k = 0; k <= NBOOT; k++){
        double wk = bootWeight[k]; if (wk <= 0) continue;
        pRefPt[k][iR][c]->Fill(nch, refMeanPt, wk);
        if (haveFsub) pRefPtFsub[k][iR][c]->Fill(nch, refMeanPt, wk);
        if (haveFcs)  pRefPtFcs [k][iR][c]->Fill(nch, refMeanPt, wk);
        for (int ip = 0; ip < nB; ip++){
          double pt = ptCenter[ip];
          pFrac[k][IVAR_F][iR][c]->Fill(nch, pt, n[ip]/sumRaw, wk);
          if (haveFsub) pFrac[k][IVAR_FSUB][iR][c]->Fill(nch, pt, (n[ip] - p[ip])/sumSub, wk);
          if (haveFcs)  pFrac[k][IVAR_FCS ][iR][c]->Fill(nch, pt, cs[ip]/sumCs, wk);
        }
      }
    }
  }
}

void Phase1::Run_mean(){
  Long64_t nEntries = fChain->GetEntries(), first, last; jobRange(nEntries, IJob, NEv, first, last);
  cout << "  mean: events [" << first << "," << last << ") of " << nEntries << ", " << NBOOT
       << " bootstrap replicas/event, " << NVAR << " probe variants" << endl;
  for (Long64_t e = first; e < last; e++){
    fChain->GetEntry(e); Tool::progress(e-first, last-first);
    Event ev; Fill_Evt(&ev); Fill_mean(&ev, e + evtOffset());   // global index -> same seed as corr
  }
  cout << "\r  mean: done (" << last-first << ")        " << endl;
}

void Phase1::Final_mean(){
  fOut->cd();
  WriteRefWin();
  for (int k = 0; k <= NBOOT; k++) for (int iv = 0; iv < NVAR; iv++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++)
    pFrac[k][iv][iR][c]->Write();
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++) pRefPt[k][iR][c]->Write();
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++) pRefPtFsub[k][iR][c]->Write();
  for (int k = 0; k <= NBOOT; k++) for (int iR = 0; iR < nR; iR++) for (int c = 0; c < NSEL; c++) pRefPtFcs[k][iR][c]->Write();
  fOut->Close();
}
