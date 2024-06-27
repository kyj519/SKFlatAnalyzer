#ifndef AnalyzerParameter_h
#define AnalyzerParameter_h

#include "TString.h"
#include <iostream>

using namespace std;

class AnalyzerParameter
{
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
  
  TString Jet_ID, FatJet_ID, PUJet_Veto_ID;
  
  //syst which affect selection
  //JES Breakdown List:       JECSources = {"Absolute", "BBEC1", "EC2", "FlavorQCD", "HF", "RelativeBal"}
  //JES Year Dependent List:  JECSources = {"Absolute", "BBEC1", "EC2", "HF", "RelativeSample"};
  enum Syst
  {
    Central,
    JetResUp, JetResDown, 
    JetEnUp, JetEnDown,
    /* MuonEnUp, MuonEnDown,  */
    ElectronResUp, ElectronResDown, 
    ElectronEnUp, ElectronEnDown, 
    UnclusteredEnergyUp, UnclusteredEnergyDown,
    NSyst,
    JetEnAbsoluteUp, JetEnAbsoluteDown,
    JetEnBBEC1Up, JetEnBBEC1Down,
    JetEnEC2Up, JetEnEC2Down,
    JetEnFlavorQCDUp, JetEnFlavorQCDDown,
    JetEnHFUp, JetEnHFDown,
    JetEnRelativeBalUp, JetEnRelativeBalDown,
    JetEnAbsolute2018Up, JetEnAbsolute2018Down,
    JetEnBBEC12018Up, JetEnBBEC12018Down,
    JetEnEC22018Up, JetEnEC22018Down,
    JetEnHF2018Up, JetEnHF2018Down,
    JetEnRelativeSample2018Up, JetEnRelativeSample2018Down,
    JetEnAbsolute2017Up, JetEnAbsolute2017Down,
    JetEnBBEC12017Up, JetEnBBEC12017Down,
    JetEnEC22017Up, JetEnEC22017Down,
    JetEnHF2017Up, JetEnHF2017Down,
    JetEnRelativeSample2017Up, JetEnRelativeSample2017Down,
    JetEnAbsolute2016Up, JetEnAbsolute2016Down,
    JetEnBBEC12016Up, JetEnBBEC12016Down,
    JetEnEC22016Up, JetEnEC22016Down,
    JetEnHF2016Up, JetEnHF2016Down,
    JetEnRelativeSample2016Up, JetEnRelativeSample2016Down
  };
  Syst syst_;
  
  TString GetSystType();

  void Clear();

  AnalyzerParameter();
  ~AnalyzerParameter();

};

#endif
