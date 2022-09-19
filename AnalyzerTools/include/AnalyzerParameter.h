#ifndef AnalyzerParameter_h
#define AnalyzerParameter_h

#include "TString.h"
#include <iostream>

using namespace std;

class AnalyzerParameter{

public:

  TString Name;

  bool MCCorrrectionIgnoreNoHist;

  TString Electron_Tight_ID, Electron_Loose_ID, Electron_Veto_ID;
  TString Electron_ID_SF_Key, Electron_Trigger_SF_Key;
  TString Electron_FR_ID, Electron_FR_Key;
  TString Electron_CF_ID, Electron_CF_Key;
  double Electron_Tight_RelIso, Electron_Loose_RelIso, Electron_Veto_RelIso;
  bool Electron_UseMini, Electron_UsePtCone;
  double Electron_MinPt;
  
  TString Muon_Tight_ID, Muon_Loose_ID, Muon_Veto_ID;
  TString Muon_RECO_SF_Key, Muon_ID_SF_Key, Muon_ISO_SF_Key, Muon_Trigger_SF_Key;
  TString Muon_FR_ID, Muon_FR_Key;
  TString Muon_CF_ID, Muon_CF_Key;
  double Muon_Tight_RelIso, Muon_Loose_RelIso, Muon_Veto_RelIso;
  bool Muon_UseMini, Muon_UsePtCone, Muon_UseTuneP;
  double Muon_MinPt;

  TString Jet_ID, FatJet_ID;

  enum Syst{
    Central,
    JetResUp, JetResDown, 
    JetEnUp, JetEnDown,
    /* MuonEnUp, MuonEnDown,  */
    /* ElectronResUp, ElectronResDown, */
    /* ElectronEnUp, ElectronEnDown, */
    PileupReweightUp, PileupReweightDown,
    PrefireReweightUp, PrefireReweightDown,
    MuTrigUp, MuTrigDown,
    MuIDUp, MuIDDown,
    MuISOUp, MuISODown,
    PileupJetVetoUp, PileupJetVetoDown,
    Up_HF, Down_HF, 
    Up_JES, Down_JES, 
    Up_LFStats1, Down_LFStats1,
    Up_LFStats2, Down_LFStats2,
    Up_CFERR1, Down_CFERR1,
    Up_CFERR2, Down_CFERR2,
    Up_HFStats1, Down_HFStats1,  
    Up_HFStats2, Down_HFStats2,
    Up_LF, Down_LF,
    NSyst
  };
  Syst syst_;
  TString GetSystType();

  void Clear();

  AnalyzerParameter();
  ~AnalyzerParameter();

};

#endif
