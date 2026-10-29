// ============================================================================
//  Phase1_evt.C -- Fill_Evt: unpack the current tree entry for acceptance IACC.
//
//  Per R (slot k = IACC*TREE_NR + iR):
//   * reference particles of reference IREF: N_ref, S1 = sum pT, S2 = sum pT^2,
//     [pT]_ref = S1/N_ref (-1 if N_ref = 0);
//   * jet-class selection: the leading jet (anti-kT, raw jetPt > 0 = jet found)
//     is put in jet class c if its pT (definition V0_JETPT: jetPtCorr, csJetPt
//     or jetPtTruth) is in [SELLO[c], SELHI[c]) AND the event passes the
//     centrality selection. The classes do not overlap, so an event is in at
//     most one class per R;
//   * for that class the spectra n_b (nIn), p_b (nPerp) and cs_b (nCS) are
//     copied (zero for every other class). The class is kept (inJetClass) only
//     if sum_b n_b > 0, which f_b needs.
//  Centrality selection: V0_CENTLO <= centr < V0_CENTHI; events with centr < 0
//  (no added UE: vacuum/PYTHIA trees without embedding) always pass.
//  The event gate on the reference (N_ref >= 2, S2 present) is applied in the
//  mean and corr fills, not here.
// ============================================================================
void Phase1::Fill_Evt(Event* evt){
  evt->weight = br_weight;
  evt->nch    = (double)br_nchMid;
  bool passCentrality = (br_centr < 0) || (br_centr >= CENTLO && br_centr < CENTHI);
  for (int iR = 0; iR < nR; iR++){
    int k = slot(iR);
    double refN = br_refN[k], refS1 = br_refS1[k], refS2 = br_refS2[k];
    evt->refCount [iR] = (int)llround(refN);
    evt->refSumPt [iR] = refS1;
    evt->refSumPt2[iR] = refS2;
    evt->refMeanPt[iR] = (refN > 0) ? refS1/refN : -1.0;

    bool   hasJet = (br_jetPt[k] > 0);
    double jetPt  = JETPTDEF == 1 ? br_csJetPt[k] : JETPTDEF == 2 ? br_jetPtTruth[k] : br_jetPtCorr[k];
    for (int c = 0; c < NSEL; c++){
      bool inClass = passCentrality && hasJet && jetPt >= SELLO[c] && jetPt < SELHI[c];
      double sumInCone = 0.0;
      for (int ip = 0; ip < nB; ip++){
        int i = k*TREE_NB + ip;
        double n  = inClass ? br_nIn[i]   : 0.;
        double p  = inClass ? br_nPerp[i] : 0.;
        double cs = inClass ? br_nCS[i]   : 0.;
        evt->nInCone[iR][c][ip] = n; evt->nPerpCone[iR][c][ip] = p; evt->nCsJet[iR][c][ip] = cs;
        sumInCone += n;
      }
      evt->inJetClass[iR][c] = inClass && (sumInCone > 0.0);
    }
  }
}
