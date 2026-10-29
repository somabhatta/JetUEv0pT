// ============================================================================
//  Phase1_ini.C -- stage "ini": QA distributions of the per-event inputs
//  (event weighted, no event gate). Used to check acceptance, binning and the
//  Nch range: hNch must not pile up in the last bin (raise V0_NCHMAX if it does).
//
//  Output <Tag>_ini_<IJob>.root:
//    hNch               Nch (nchMid)
//    hLeadHi            tree centr (percent; -1 = none)
//    hNref_R<i>         N_ref of reference IREF
//    hRef_R<i>          [pT]_ref (events with N_ref > 0)
//    hNcone_R<i>        leading-jet jetPtCorr (events with a jet)
//    hNin_R<i>          sum_b n_b (in-cone yield, all events)
//    hJetFlag_R<i>      jet class index of the event, -1 = none
//    hRefVsNch_R<i>     [pT]_ref vs Nch
// ============================================================================
void Phase1::InitHists_ini(){
  hNch        = new TH1D("hNch", ";N_{ch} (|#eta|<0.8, 0.2<p_{T}<5 GeV);events", NCHBINS, 0, NCHMAX);
  hCentrality = new TH1D("hLeadHi", ";centrality (JEWEL, %);events", 101, -1, 100);
  for (int iR = 0; iR < nR; iR++){
    hRefCount[iR]    = new TH1D(Form("hNref_R%d", iR),  Form("R=%.1f;N_{ref} (reference-region count);events", Rval[iR]), 400, -0.5, 399.5);
    hRefMeanPt[iR]   = new TH1D(Form("hRef_R%d", iR),   Form("R=%.1f;[p_{T}]_{ref} (GeV/c);events", Rval[iR]), 120, 0, 3);
    hLeadJetPt[iR]   = new TH1D(Form("hNcone_R%d", iR), Form("R=%.1f;UE-corrected leading-jet p_{T} (GeV/c);events", Rval[iR]), 150, -30, 270);
    hInConeYield[iR] = new TH1D(Form("hNin_R%d", iR),   Form("R=%.1f;in-cone charged yield;events", Rval[iR]), 200, -0.5, 199.5);
    hJetClass[iR]    = new TH1D(Form("hJetFlag_R%d", iR), Form("R=%.1f;jet class (-1 none);events", Rval[iR]), NSEL+1, -1.5, (double)NSEL - 0.5);
    hRefMeanPtVsNch[iR] = new TH2D(Form("hRefVsNch_R%d", iR), Form("R=%.1f;N_{ch};[p_{T}]_{ref}", Rval[iR]),
                                   NCHBINS, 0, NCHMAX, 120, 0, 3);
  }
}

void Phase1::Fill_ini(Event* evt){
  double w = evt->weight;
  hNch->Fill(evt->nch, w);
  hCentrality->Fill(br_centr, w);
  for (int iR = 0; iR < nR; iR++){
    int k = slot(iR);
    hRefCount[iR]->Fill(evt->refCount[iR], w);
    if (evt->refMeanPt[iR] >= 0){
      hRefMeanPt[iR]->Fill(evt->refMeanPt[iR], w);
      hRefMeanPtVsNch[iR]->Fill(evt->nch, evt->refMeanPt[iR], w);
    }
    if (br_jetPt[k] > 0) hLeadJetPt[iR]->Fill(br_jetPtCorr[k], w);
    double inConeYield = 0; for (int ip = 0; ip < nB; ip++) inConeYield += br_nIn[k*TREE_NB + ip];
    hInConeYield[iR]->Fill(inConeYield, w);
    int jetClass = -1; for (int c = 0; c < NSEL; c++) if (evt->inJetClass[iR][c]) jetClass = c;
    hJetClass[iR]->Fill(jetClass, w);
  }
}

void Phase1::Run_ini(){
  Long64_t nEntries = fChain->GetEntries(), first, last; jobRange(nEntries, IJob, NEv, first, last);
  cout << "  ini: events [" << first << "," << last << ") of " << nEntries << endl;
  for (Long64_t e = first; e < last; e++){
    fChain->GetEntry(e); Tool::progress(e-first, last-first);
    Event ev; Fill_Evt(&ev); Fill_ini(&ev);
  }
  cout << "\r  ini: done (" << last-first << ")        " << endl;
}

void Phase1::Final_ini(){
  fOut->cd();
  hNch->Write(); hCentrality->Write();
  for (int iR = 0; iR < nR; iR++){
    hRefCount[iR]->Write(); hRefMeanPt[iR]->Write(); hLeadJetPt[iR]->Write();
    hInConeYield[iR]->Write(); hJetClass[iR]->Write(); hRefMeanPtVsNch[iR]->Write();
  }
  fOut->Close();
}
