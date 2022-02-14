#ifndef PUJets_MCTruthMatcher_h
#define PUJets_MCTruthMatcher_h

#include "AnalyzerCore.h"


class PUJets_MCTruthMatcher : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();
  int GetGenMatchedJet(const Jet& jet, const std::vector<Gen>& gens);

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
  vector<Gen> vec_gens;




  PUJets_MCTruthMatcher();
  ~PUJets_MCTruthMatcher();

};



#endif