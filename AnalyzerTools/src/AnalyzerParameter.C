#include "AnalyzerParameter.h"

void AnalyzerParameter::Clear()
{
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
  PUJet_Veto_ID = "";

  syst_ = Central;
} // void AnalyzerParameter::Clear()

//////////

AnalyzerParameter::AnalyzerParameter()
{
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
} // AnalyzerParameter::AnalyzerParameter()

TString AnalyzerParameter::GetSystType()
{
  if (syst_ == Syst::Central)
    return "Central";
  else if (syst_ == Syst::JetResUp)
    return "JetResUp";
  else if (syst_ == Syst::JetResDown)
    return "JetResDown";
  else if (syst_ == Syst::JetEnUp)
    return "JetEnUp";
  else if (syst_ == Syst::JetEnDown)
    return "JetEnDown";
  // else if(syst_==Syst::MuonEnUp) return "MuonEnUp";
  // else if(syst_==Syst::MuonEnDown) return "MuonEnDown";
  else if (syst_ == Syst::ElectronResUp)
    return "ElectronResUp";
  else if (syst_ == Syst::ElectronResDown)
    return "ElectronResDown";
  else if (syst_ == Syst::ElectronEnUp)
    return "ElectronEnUp";
  else if (syst_ == Syst::ElectronEnDown)
    return "ElectronEnDown";
  else if (syst_ == Syst::UnclusteredEnergyUp)
    return "UEUp";
  else if (syst_ == Syst::UnclusteredEnergyDown)
    return "UEDown";
  //JES Breakdown: Total 42 variations
  else if(syst_ == Syst::JetEnAbsoluteUp)
    return "JetEnAbsoluteUp";
  else if(syst_ == Syst::JetEnAbsoluteDown)
   return "JetEnAbsoluteDown";
  else if(syst_ == Syst::JetEnBBEC1Up)
    return "JetEnBBEC1Up";
  else if(syst_ == Syst::JetEnBBEC1Down)
    return "JetEnBBEC1Down";
  else if(syst_ == Syst::JetEnEC2Up)
    return "JetEnEC2Up";
  else if(syst_ == Syst::JetEnEC2Down)
    return "JetEnEC2Down";
  else if(syst_ == Syst::JetEnFlavorQCDUp)
    return "JetEnFlavorQCDUp";
  else if(syst_ == Syst::JetEnFlavorQCDDown)
    return "JetEnFlavorQCDDown";
  else if(syst_ == Syst::JetEnHFUp)
    return "JetEnHFUp";
  else if(syst_ == Syst::JetEnHFDown)
    return "JetEnHFDown";
  else if(syst_ == Syst::JetEnRelativeBalUp)
    return "JetEnRelativeBalUp";
  else if(syst_ == Syst::JetEnRelativeBalDown)
    return "JetEnRelativeBalDown";
  else if(syst_ == Syst::JetEnAbsolute2018Up)
    return "JetEnAbsolute2018Up";
  else if(syst_ == Syst::JetEnAbsolute2018Down)
    return "JetEnAbsolute2018Down";
  else if(syst_ == Syst::JetEnBBEC12018Up)
    return "JetEnBBEC12018Up";
  else if(syst_ == Syst::JetEnBBEC12018Down)
    return "JetEnBBEC12018Down";
  else if(syst_ == Syst::JetEnEC22018Up)
    return "JetEnEC22018Up";
  else if(syst_ == Syst::JetEnEC22018Down)
    return "JetEnEC22018Down";
  else if(syst_ == Syst::JetEnHF2018Up)
    return "JetEnHF2018Up";
  else if(syst_ == Syst::JetEnHF2018Down)
    return "JetEnHF2018Down";
  else if(syst_ == Syst::JetEnRelativeSample2018Up)
    return "JetEnRelativeSample2018Up";
  else if(syst_ == Syst::JetEnRelativeSample2018Down)
    return "JetEnRelativeSample2018Down";
    
  else if(syst_ == Syst::JetEnRelativeBalUp)
    return "JetEnRelativeBalUp";
  else if(syst_ == Syst::JetEnRelativeBalDown)
    return "JetEnRelativeBalDown";
  else if(syst_ == Syst::JetEnAbsolute2017Up)
    return "JetEnAbsolute2017Up";
  else if(syst_ == Syst::JetEnAbsolute2017Down)
    return "JetEnAbsolute2017Down";
  else if(syst_ == Syst::JetEnBBEC12017Up)
    return "JetEnBBEC12017Up";
  else if(syst_ == Syst::JetEnBBEC12017Down)
    return "JetEnBBEC12017Down";
  else if(syst_ == Syst::JetEnEC22017Up)
    return "JetEnEC22017Up";
  else if(syst_ == Syst::JetEnEC22017Down)
    return "JetEnEC22017Down";
  else if(syst_ == Syst::JetEnHF2017Up)
    return "JetEnHF2017Up";
  else if(syst_ == Syst::JetEnHF2017Down)
    return "JetEnHF2017Down";
  else if(syst_ == Syst::JetEnRelativeSample2017Up)
    return "JetEnRelativeSample2017Up";
  else if(syst_ == Syst::JetEnRelativeSample2017Down)
    return "JetEnRelativeSample2017Down";

  else if(syst_ == Syst::JetEnRelativeBalUp)
    return "JetEnRelativeBalUp";
  else if(syst_ == Syst::JetEnRelativeBalDown)
    return "JetEnRelativeBalDown";
  else if(syst_ == Syst::JetEnAbsolute2016Up)
    return "JetEnAbsolute2016Up";
  else if(syst_ == Syst::JetEnAbsolute2016Down)
    return "JetEnAbsolute2016Down";
  else if(syst_ == Syst::JetEnBBEC12016Up)
    return "JetEnBBEC12016Up";
  else if(syst_ == Syst::JetEnBBEC12016Down)
    return "JetEnBBEC12016Down";
  else if(syst_ == Syst::JetEnEC22016Up)
    return "JetEnEC22016Up";
  else if(syst_ == Syst::JetEnEC22016Down)
    return "JetEnEC22016Down";
  else if(syst_ == Syst::JetEnHF2016Up)
    return "JetEnHF2016Up";
  else if(syst_ == Syst::JetEnHF2016Down)
    return "JetEnHF2016Down";
  else if(syst_ == Syst::JetEnRelativeSample2016Up)
    return "JetEnRelativeSample2016Up";
  else if(syst_ == Syst::JetEnRelativeSample2016Down)
    return "JetEnRelativeSample2016Down";


  cout << "[AnalyzerParameter::GetSystType] Wrong Syst " << syst_ << endl;

  exit(ENODATA);
  return "ERROR";
}

AnalyzerParameter::~AnalyzerParameter()
{
}
