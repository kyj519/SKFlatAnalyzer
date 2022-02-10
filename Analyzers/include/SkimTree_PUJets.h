#ifndef SkimTree_PUJets_h
#define SkimTree_PUJets_h

#include "AnalyzerCore.h"


class SkimTree_PUJets : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();
  void WriteHist();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  TString IsoMuTriggerName;
  double TriggerSafePtCut;
  double weight_Prefire;

  vector<TString> MuonIDs, MuonIDSFKeys, MuonISOSFKeys,MuonTrigSFKeys, vec_muon_trigger, ElectronIDs;
  vector<Muon> vec_muon;
  vector<Jet> vec_jet;
  vector<Electron> vec_electron;




  SkimTree_PUJets();
  ~SkimTree_PUJets();

protected:
 TTree* newtree;
};



#endif