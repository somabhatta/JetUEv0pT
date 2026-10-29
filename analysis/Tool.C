// ============================================================================
//  Tool.C -- small helpers shared by Phase1 and Phase2.
//    makeChain : TChain of the tree "t_<V0_TREE>" from one .root file or a file
//                list (one path per line, '#' = comment; if a line has three
//                columns and the third ends in .root, that column is used)
//    progress  : percent counter on stdout
//    readCuts / writeCuts : wide-centrality Nch edges (INPUT/cutNch_<sample>.txt,
//                optional input of Phase2)
// ============================================================================
#ifndef Tool_C_
#define Tool_C_
#include "Rule.h"
class Tool {
public:
  Tool(){}
  static TChain* makeChain(const string& input){
    TChain* chain = new TChain(treeName().c_str());
    if (input.size() > 5 && input.substr(input.size()-5) == ".root"){ chain->Add(input.c_str()); }
    else {
      ifstream list(input.c_str()); string line; int nFiles = 0;
      while (getline(list, line)){
        if (line.empty() || line[0] == '#') continue;
        istringstream cols(line); string c1, c2, c3; cols >> c1 >> c2 >> c3;
        string path = (!c3.empty() && c3.substr(c3.size()-5) == ".root") ? c3 : c1;
        chain->Add(path.c_str()); nFiles++;
      }
      cout << "Tool: " << nFiles << " files from " << input << endl;
    }
    return chain;
  }
  static void progress(Long64_t i, Long64_t n){
    static int last = -1; int pct = n > 0 ? int(100.*i/n) : 0;
    if (pct != last){ cout << "\r  " << pct << "%  " << flush; last = pct; }
  }
  static void writeCuts(const string& fn, const double* edges, int nEdges){
    ofstream o(fn.c_str()); for (int i = 0; i < nEdges; i++) o << edges[i] << "\n";
  }
  static bool readCuts(const string& fn, double* edges, int nEdges){
    ifstream in(fn.c_str()); if (!in) return false;
    for (int i = 0; i < nEdges; i++){ if (!(in >> edges[i])) return false; }
    return true;
  }
};
#endif
