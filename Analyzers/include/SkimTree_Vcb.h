#ifndef __SkimTree_Vcb_h__
#define __SkimTree_Vcb_h__

#include "AnalyzerCore.h"

#include "Vcb_Def.h"

using namespace std;

class SkimTree_Vcb : public AnalyzerCore
{
 public:
  SkimTree_Vcb();
  ~SkimTree_Vcb();

  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter(AnalyzerParameter param);
  void WriteHist();

 protected:
  TTree* newtree;

  vector<TString> vec_single_muon_trigger;
  vector<TString> vec_single_electron_trigger;

  float mu_trig_safe_pt_cut;
  float el_trig_safe_pt_cut;

  vector<Muon> vec_muon;
  vector<Electron> vec_electron;
  vector<Jet> vec_jet;
};

#endif /* SkimTree_Vcb_h */
