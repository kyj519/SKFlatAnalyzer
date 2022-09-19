#include "AnalyzerParameter.h"

void AnalyzerParameter::Clear(){

  Name = "";

  MCCorrrectionIgnoreNoHist = false;

  Electron_Tight_ID = "";
  Electron_Loose_ID = "";
  Electron_Veto_ID = "";
  Electron_ID_SF_Key = "";
  Electron_FR_ID = "";
  Electron_FR_Key = "";
  Electron_CF_ID = "";
  Electron_CF_Key = "";
  Electron_Tight_RelIso = 999.;
  Electron_Loose_RelIso = 999.;
  Electron_Veto_RelIso = 999.;
  Electron_UseMini = false;
  Electron_UsePtCone = false;
  Electron_MinPt = 10.;

  Muon_Tight_ID = "";

  Muon_Loose_ID = "";
  Muon_Veto_ID = "";
  Muon_RECO_SF_Key = "";
  Muon_ID_SF_Key = "";
  Muon_ISO_SF_Key = "";
  Muon_Trigger_SF_Key = "";
  Muon_FR_ID = "";
  Muon_FR_Key = "";
  Muon_CF_ID = "";
  Muon_CF_Key = "";
  Muon_Tight_RelIso = 999.;
  Muon_Loose_RelIso = 999.;
  Muon_Veto_RelIso = 999.;
  Muon_UseMini = false;
  Muon_UsePtCone = false;
  Muon_UseTuneP = false;
  Muon_MinPt = 10.;

  Jet_ID = "";
  FatJet_ID = "";

  syst_ = Central;

}

AnalyzerParameter::AnalyzerParameter(){

  Name = "Default";

  MCCorrrectionIgnoreNoHist = false;

  Electron_Tight_ID = "passTightID";
  Electron_Loose_ID = "passLooseID";
  Electron_Veto_ID = "passVetoID";
  Electron_ID_SF_Key = "passTightID";

  Muon_Tight_ID = "POGTightWithTightIso";
  Muon_Loose_ID = "POGLoose";
  Muon_Veto_ID = "POGLoose";
  Muon_RECO_SF_Key = "";
  Muon_ID_SF_Key = "NUM_TightID_DEN_genTracks";
  Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";
  Muon_Trigger_SF_Key = "POGTight";

  Jet_ID = "HN";
  FatJet_ID = "HN";

  syst_ = Central;

}

TString AnalyzerParameter::GetSystType(){

  if(syst_==Syst::Central) return "Central";
  else if(syst_==Syst::JetResUp) return "JetResUp";
  else if(syst_==Syst::JetResDown) return "JetResDown";
  else if(syst_==Syst::JetEnUp) return "JetEnUp";
  else if(syst_==Syst::JetEnDown) return "JetEnDown";
  // else if(syst_==Syst::MuonEnUp) return "MuonEnUp";
  // else if(syst_==Syst::MuonEnDown) return "MuonEnDown";
  // else if(syst_==Syst::ElectronResUp){
  //   return "ElectronResUp";
  // }
  // else if(syst_==Syst::ElectronResDown){
  //   return "ElectronResDown";
  // }
  // else if(syst_==Syst::ElectronEnUp){
  //   return "ElectronEnUp";
  // }
  // else if(syst_==Syst::ElectronEnDown){
  //   return "ElectronEnDown";
  // }
  else if(syst_==Syst::PileupReweightUp) return "PileupReweightUp";
  else if(syst_==Syst::PileupReweightDown) return "PileupReweightDown";
  else if(syst_==Syst::PrefireReweightUp) return "PrefireReweightUp";
  else if(syst_==Syst::PrefireReweightDown) return "PrefireReweightDown";
  else if(syst_==Syst::MuTrigUp) return "MuTrigUp";
  else if(syst_==Syst::MuTrigDown) return "MuTrigDown";
  else if(syst_==Syst::MuIDUp) return "MuIDUp";
  else if(syst_==Syst::MuIDDown) return "MuIDDown";
  else if(syst_==Syst::MuISOUp) return "MuISOUp";
  else if(syst_==Syst::MuISODown) return "MuISODown";
  else if(syst_==Syst::PileupJetVetoUp) return "PileupJetVetoUp";
  else if(syst_==Syst::PileupJetVetoDown) return "PileupJetVetoDown";
  else if(syst_==Syst::Up_HF) return "Up_HF";
  else if(syst_==Syst::Down_HF) return "Down_HF";
  else if(syst_==Syst::Up_JES) return "Up_JES";
  else if(syst_==Syst::Down_JES) return "Down_JES";
  else if(syst_==Syst::Up_LFStats1) return "Up_LFStats1";
  else if(syst_==Syst::Down_LFStats1) return "Down_LFStats1";
  else if(syst_==Syst::Up_LFStats2) return "Up_LFStats2";
  else if(syst_==Syst::Down_LFStats2) return "Down_LFStats2";
  else if(syst_==Syst::Up_CFERR1) return "Up_CFERR1";
  else if(syst_==Syst::Down_CFERR1) return "Down_CFERR1";
  else if(syst_==Syst::Up_CFERR2) return "Up_CFERR2";
  else if(syst_==Syst::Down_CFERR2) return "Down_CFERR2";
  else if(syst_==Syst::Up_HFStats1) return "Up_HFStats1";
  else if(syst_==Syst::Down_HFStats1) return "Down_HFStats1";
  else if(syst_==Syst::Up_HFStats2) return "Up_HFStats2";
  else if(syst_==Syst::Down_HFStats2) return "Down_HFStats2";
  else if(syst_==Syst::Up_LF) return "Up_LF";
  else if(syst_==Syst::Down_LF) return "Down_LF";

  else{
    cout << "[AnalyzerParameter::GetSystType] Wrong Syst " << syst_ << endl;
    exit(ENODATA);
    return "ERROR";
  }

}

AnalyzerParameter::~AnalyzerParameter(){
}
