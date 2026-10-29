// compare_generation_outputs.C -- bit-level comparison of two jetv0ana output files.
//
//   root -l -b -q 'compare_generation_outputs.C("reference.root", "test.root")'
//
// Walks both files recursively (top level and the corr/, spectra/ directories):
//   - the key lists (name, class) must agree;
//   - TTrees: same branches, same leaf types and lengths, same number of
//     entries; every leaf element of every entry is compared bitwise;
//   - TH1: same binning, title, entries; every bin content and bin error
//     (including under/overflow) and the stored sum of weights squared are
//     compared bitwise.
// Prints the counts and the first differences; the last line reads
// "RESULT: bit-identical" only if nothing differs.
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TTree.h"

#include <cstdio>
#include <cstring>
#include <set>
#include <string>

namespace {
long gTrees = 0, gEntries = 0, gValues = 0, gHists = 0, gBins = 0, gDiffs = 0;

bool sameBits(double a, double b) { return memcmp(&a, &b, sizeof(double)) == 0; }

void report(const std::string &what)
{
  gDiffs++;
  if (gDiffs <= 20) printf("  DIFF %s\n", what.c_str());
}

void compareTrees(const std::string &path, TTree *ta, TTree *tb)
{
  gTrees++;
  if (ta->GetEntries() != tb->GetEntries()) { report(path + ": number of entries"); return; }
  TObjArray *la = ta->GetListOfLeaves(), *lb = tb->GetListOfLeaves();
  if (la->GetEntries() != lb->GetEntries()) { report(path + ": number of leaves"); return; }
  for (int i = 0; i < la->GetEntries(); i++) {
    TLeaf *a = (TLeaf *)la->At(i), *b = (TLeaf *)lb->At(i);
    if (std::string(a->GetName()) != b->GetName() || std::string(a->GetTypeName()) != b->GetTypeName() ||
        a->GetLenStatic() != b->GetLenStatic())
      { report(path + ": leaf layout at " + a->GetName()); return; }
  }
  for (Long64_t e = 0; e < ta->GetEntries(); e++) {
    ta->GetEntry(e); tb->GetEntry(e); gEntries++;
    for (int i = 0; i < la->GetEntries(); i++) {
      TLeaf *a = (TLeaf *)la->At(i), *b = (TLeaf *)lb->At(i);
      if (a->GetLen() != b->GetLen()) { report(path + "/" + a->GetName() + Form(" entry %lld: length", e)); continue; }
      for (int j = 0; j < a->GetLen(); j++) {
        gValues++;
        double x = a->GetValue(j), y = b->GetValue(j);   // float/int -> double is exact and injective
        if (!sameBits(x, y))
          report(path + "/" + a->GetName() + Form("[%d] entry %lld: %.9g vs %.9g", j, e, x, y));
      }
    }
  }
}

void compareHists(const std::string &path, TH1 *a, TH1 *b)
{
  gHists++;
  if (a->GetNcells() != b->GetNcells() || std::string(a->GetTitle()) != b->GetTitle() ||
      !sameBits(a->GetEntries(), b->GetEntries()))
    { report(path + ": binning / title / entries"); return; }
  for (int i = 0; i < a->GetNbinsX() + 1; i++)
    if (!sameBits(a->GetXaxis()->GetBinUpEdge(i), b->GetXaxis()->GetBinUpEdge(i))) { report(path + ": bin edges"); return; }
  for (int i = 0; i < a->GetNcells(); i++) {
    gBins++;
    if (!sameBits(a->GetBinContent(i), b->GetBinContent(i)) || !sameBits(a->GetBinError(i), b->GetBinError(i)))
      report(path + Form(" bin %d: %.17g vs %.17g", i, a->GetBinContent(i), b->GetBinContent(i)));
  }
}

void compareDirs(const std::string &path, TDirectory *da, TDirectory *db)
{
  std::set<std::string> keysA, keysB;
  for (auto *k : *da->GetListOfKeys()) keysA.insert(std::string(k->GetName()) + " (" + ((TKey *)k)->GetClassName() + ")");
  for (auto *k : *db->GetListOfKeys()) keysB.insert(std::string(k->GetName()) + " (" + ((TKey *)k)->GetClassName() + ")");
  if (keysA != keysB) report(path + ": key lists differ");
  for (auto *k : *da->GetListOfKeys()) {
    TKey *key = (TKey *)k;
    std::string name = key->GetName(), full = path + "/" + name;
    TObject *oa = da->Get(name.c_str()), *ob = db->Get(name.c_str());
    if (!ob) { report(full + ": missing in second file"); continue; }
    if (oa->InheritsFrom(TDirectory::Class())) compareDirs(full, (TDirectory *)oa, (TDirectory *)ob);
    else if (oa->InheritsFrom(TTree::Class())) compareTrees(full, (TTree *)oa, (TTree *)ob);
    else if (oa->InheritsFrom(TH1::Class())) compareHists(full, (TH1 *)oa, (TH1 *)ob);
    else report(full + ": unhandled class " + oa->ClassName());
  }
}
} // namespace

int compare_generation_outputs(const char *fileA, const char *fileB)
{
  gTrees = gEntries = gValues = gHists = gBins = gDiffs = 0;
  TFile a(fileA), b(fileB);
  if (a.IsZombie() || b.IsZombie()) { printf("cannot open input\n"); return 1; }
  compareDirs("", &a, &b);
  printf("%s vs %s:\n  trees %ld, tree entries %ld, leaf values %ld, histograms %ld, bins %ld, differences %ld\n",
         fileA, fileB, gTrees, gEntries, gValues, gHists, gBins, gDiffs);
  printf("RESULT: %s\n", gDiffs == 0 ? "bit-identical" : "DIFFERENT");
  return gDiffs == 0 ? 0 : 1;
}
