#ifndef PUJets_h
#define PUJets_h

#include "AnalyzerCore.h"


class PUJets : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

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


  unsigned int getBinID(const double eta, const double Pt);
  vector<TString> getHistKey(bool isDATA, bool pt_balanced, unsigned int binID, int PUID);

  PUJets();
  ~PUJets();

};



#endif