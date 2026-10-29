// compare_analysis_outputs.C: bin-by-bin comparison of every object (highest cycle) in two
// ROOT files (Phase1 mean/corr or Phase2 outputs). Used to check that a code change leaves
// the analysis output unchanged.
// Run: root -l -b -q 'compare_analysis_outputs.C+("old.root","new.root")'
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TNamed.h"
#include "TClass.h"
#include <set>
#include <cstdio>
#include <cmath>
int compare_analysis_outputs(const char* fa, const char* fb, int verbose=5){
  TFile* A=TFile::Open(fa); TFile* B=TFile::Open(fb);
  if(!A||!B||A->IsZombie()||B->IsZombie()){ printf("OPENFAIL %s %s\n",fa,fb); return 1; }
  std::set<std::string> na, nb;
  for(auto k:*A->GetListOfKeys()) na.insert(k->GetName());
  for(auto k:*B->GetListOfKeys()) nb.insert(k->GetName());
  int onlyA=0, onlyB=0; for(auto&s:na) if(!nb.count(s)){ onlyA++; if(onlyA<=verbose) printf("  only in A: %s\n",s.c_str()); }
  for(auto&s:nb) if(!na.count(s)){ onlyB++; if(onlyB<=verbose) printf("  only in B: %s\n",s.c_str()); }
  long nObj=0,nHist=0,nBins=0,nDiffObj=0,nSkip=0; double maxAbs=0, maxRel=0; std::string worst;
  for(auto& s:na){ if(!nb.count(s)) continue; nObj++;
    TObject* a=A->Get(s.c_str()); TObject* b=B->Get(s.c_str());
    if(strcmp(a->ClassName(),b->ClassName())){ printf("  class differs %s\n",s.c_str()); nDiffObj++; continue; }
    if(a->InheritsFrom("TH1")){ TH1* ha=(TH1*)a; TH1* hb=(TH1*)b; nHist++;
      int n=ha->GetNcells(); bool diff=false;
      if(n!=hb->GetNcells()){ printf("  ncells differ %s\n",s.c_str()); nDiffObj++; continue; }
      if(strcmp(ha->GetTitle(),hb->GetTitle())){ printf("  title differs %s: '%s' vs '%s'\n",s.c_str(),ha->GetTitle(),hb->GetTitle()); diff=true; }
      if(ha->GetEntries()!=hb->GetEntries()) diff=true;
      bool prof2=a->InheritsFrom("TProfile2D"), prof=a->InheritsFrom("TProfile");
      for(int i=0;i<n;i++){ nBins++;
        double va[4]={ha->GetBinContent(i),ha->GetBinError(i),0,0}, vb[4]={hb->GetBinContent(i),hb->GetBinError(i),0,0};
        if(prof){ va[2]=((TProfile*)ha)->GetBinEntries(i); vb[2]=((TProfile*)hb)->GetBinEntries(i); }
        if(prof2){ va[2]=((TProfile2D*)ha)->GetBinEntries(i); vb[2]=((TProfile2D*)hb)->GetBinEntries(i); }
        if(ha->GetSumw2N()&&hb->GetSumw2N()){ va[3]=ha->GetSumw2()->At(i); vb[3]=hb->GetSumw2()->At(i); }
        for(int q=0;q<4;q++){ bool na_=std::isnan(va[q]), nb_=std::isnan(vb[q]);
          if(na_&&nb_) continue;
          if(va[q]!=vb[q]){ diff=true; double d=fabs(va[q]-vb[q]); double r=d/std::max(fabs(va[q]),fabs(vb[q]));
            if(d>maxAbs){maxAbs=d; worst=s;} if(r>maxRel) maxRel=r; } } }
      if(diff){ nDiffObj++; if(nDiffObj<=verbose) printf("  DIFF %s\n",s.c_str()); }
    } else if(a->InheritsFrom("TNamed") && !a->InheritsFrom("TPad")){
      if(strcmp(((TNamed*)a)->GetTitle(),((TNamed*)b)->GetTitle())){ nDiffObj++; printf("  TNamed title differs %s\n",s.c_str()); }
    } else nSkip++;
  }
  printf("RESULT %s : keysA=%zu keysB=%zu onlyA=%d onlyB=%d common=%ld hist=%ld bins=%ld diffObj=%ld skipped(non-hist)=%ld maxAbs=%g maxRel=%g %s\n",
         fb, na.size(), nb.size(), onlyA, onlyB, nObj, nHist, nBins, nDiffObj, nSkip, maxAbs, maxRel, worst.c_str());
  return (onlyA||onlyB||nDiffObj)?1:0;
}
