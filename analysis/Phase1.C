// ============================================================================
//  Phase1.C -- driver of the staged Phase1 event loop (ini / mean / corr).
//
//  Build : root -l -b -q -e '.L Phase1.C+O'           (NBOOT 200)
//          ACLiC with -DV0_NBOOT=40 for the 40-replica production build, e.g.
//          root -l -b -q -e 'gSystem->AddIncludePath("-DV0_NBOOT=40"); .L Phase1.C+O'
//          (or set "#define V0_NBOOT 40" before including Rule.h)
//  Run   : gSystem->Load("Phase1_C.so");
//          Phase1 p("<files.list | file.root>", "<tag>", "<outdir>", <ijob>, <nev>);
//          p.Exec_mean();   // then hadd the slices, then Exec_corr() with V0_MEANFILE
//  Input : TChain of "t_<V0_TREE>"; run options from the environment (Rule.h).
//  Output: <outdir>/<tag>_<stage>_<ijob>.root (object names in Phase1_<stage>.C).
// ============================================================================
#define Phase1_cxx
#include "Phase1.h"
#include "Phase1_evt.C"
#include "Phase1_ini.C"
#include "Phase1_mean.C"
#include "Phase1_corr.C"

Phase1::Phase1(string input, string tag, string outdir, int ijob, int nevts){
  InInput = input; Tag = tag; OutDir = outdir; IJob = ijob; NEv = nevts;
  InitTree(); ReadAxis();
  cout << "Phase1 [" << Tag << "] nR=" << nR << " nProbeBins=" << nB
       << " NVAR=" << NVAR << " entries=" << fChain->GetEntries() << " NBOOT=" << NBOOT << endl;
}

// Enable and attach only the branches the stages use. Reference branches are
// <REFBR[IREF]>{N,S1,S2}; jetPtTruth only when it defines the jet classes.
void Phase1::InitTree(){
  fChain = Tool::makeChain(InInput);
  fChain->SetBranchStatus("*", 0);
  const char* used[] = {"weight", "centr", "nchMid", "jetPt", "jetPtCorr", "csJetPt", "nIn", "nPerp", "nCS"};
  for (auto b : used) fChain->SetBranchStatus(b, 1);
  fChain->SetBranchAddress("weight",    &br_weight);
  fChain->SetBranchAddress("centr",     &br_centr);
  fChain->SetBranchAddress("nchMid",    &br_nchMid);
  fChain->SetBranchAddress("jetPt",     br_jetPt);
  fChain->SetBranchAddress("jetPtCorr", br_jetPtCorr);
  fChain->SetBranchAddress("csJetPt",   br_csJetPt);
  fChain->SetBranchAddress("nIn",       br_nIn);
  fChain->SetBranchAddress("nPerp",     br_nPerp);
  fChain->SetBranchAddress("nCS",       br_nCS);
  if (JETPTDEF == 2){ fChain->SetBranchStatus("jetPtTruth", 1); fChain->SetBranchAddress("jetPtTruth", br_jetPtTruth); }
  const char* ref = REFBR[IREF];
  for (auto q : {"N", "S1", "S2"}) fChain->SetBranchStatus(Form("%s%s", ref, q), 1);
  fChain->SetBranchAddress(Form("%sN",  ref), br_refN);
  fChain->SetBranchAddress(Form("%sS1", ref), br_refS1);
  fChain->SetBranchAddress(Form("%sS2", ref), br_refS2);
  cout << "  tree " << treeName() << "  acceptance |eta|<" << TREE_ACCS[IACC]
       << "  reference IREF=" << IREF << " (" << REFN[IREF] << ", " << PTREF_LO[IREF] << "-" << PTREF_HI[IREF] << " GeV/c)"
       << "  jet pT " << JETPTN[JETPTDEF] << "  centr [" << CENTLO << "," << CENTHI << ")\n";
}

void Phase1::ReadAxis(){
  nR = TREE_NR; nB = TREE_NB;
  for (int iR = 0; iR < nR; iR++) Rval[iR] = TREE_RVALS[iR];
  for (int ip = 0; ip < nB; ip++){
    ptLo[ip] = TREE_PTEDGE[ip]; ptHi[ip] = TREE_PTEDGE[ip+1];
    ptCenter[ip] = 0.5*(ptLo[ip] + ptHi[ip]); ptEdges[ip] = ptLo[ip];
  }
  ptEdges[nB] = ptHi[nB-1];
}

void Phase1::Exec_ini(){
  Phase = "ini"; TStopwatch timer; timer.Start();
  fOut = new TFile(Form("%s/%s_ini_%d.root", OutDir.c_str(), Tag.c_str(), IJob), "RECREATE");
  InitHists_ini(); Run_ini(); Final_ini();
  timer.Stop(); cout << "[ini " << Tag << "] "; timer.Print();
}
void Phase1::Exec_mean(){
  Phase = "mean"; TStopwatch timer; timer.Start();
  fOut = new TFile(Form("%s/%s_mean_%d.root", OutDir.c_str(), Tag.c_str(), IJob), "RECREATE");
  InitHists_mean(); Run_mean(); Final_mean();
  timer.Stop(); cout << "[mean " << Tag << "] "; timer.Print();
}
void Phase1::Exec_corr(){
  Phase = "corr"; TStopwatch timer; timer.Start();
  ReadHists_corr();
  fOut = new TFile(Form("%s/%s_corr_%d.root", OutDir.c_str(), Tag.c_str(), IJob), "RECREATE");
  InitHists_corr(); Run_corr(); Final_corr();
  timer.Stop(); cout << "[corr " << Tag << ", " << NBOOT << " bootstrap replicas] "; timer.Print();
}
