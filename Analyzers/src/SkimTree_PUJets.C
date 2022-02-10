#include "SkimTree_PUJets.h"

SkimTree_PUJets::SkimTree_PUJets(){
  newtree = NULL;
}

SkimTree_PUJets::~SkimTree_PUJets(){

}

void SkimTree_PUJets::initializeAnalyzer(){
  outfile->cd();
  cout << "[SkimTree_Dilepton::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

	MuonIDs = {"POGLooseWithRelIso04","POGTight"};
  ElectronIDs = {"passLooseID","passTightID"};
	MuonIDSFKeys = {"NUM_LooseID_DEN_TrackerMuons"};
  MuonISOSFKeys = {"NUM_TightRelIso_DEN_TightIDandIPCut"};
  MuonTrigSFKeys = {"IsoMu24"};
  if(DataYear==2016){
    IsoMuTriggerName = "HLT_IsoMu24_v";
    TriggerSafePtCut = 26.;
  }
  else if(DataYear==2017){
    IsoMuTriggerName = "HLT_IsoMu27_v";
    TriggerSafePtCut = 29.;
  }
	else if(DataYear==2018){
    IsoMuTriggerName = "HLT_IsoMu24_v";
    TriggerSafePtCut = 26.;
  }
}

void SkimTree_PUJets::executeEvent(){
	vec_electron = GetAllElectrons();
	vec_muon = GetAllMuons();
	vec_jet = GetAllJets();
  weight_Prefire = GetPrefireWeight(0);
  AnalyzerParameter param;

  TString MuonID = MuonIDs.at(0);
  TString MuonIDSFKey = MuonIDSFKeys.at(0);
  TString MuonTrigSFKey = MuonTrigSFKeys.at(0);
  TString MuonISOSFKey = MuonISOSFKeys.at(0);
  param.Clear();
  param.syst_ = AnalyzerParameter::Central;
  param.Name = MuonID+"_"+"Central";

  param.Muon_Loose_ID = MuonID;
  param.Muon_Tight_ID = MuonIDs.at(1);
  param.Muon_ID_SF_Key = MuonIDSFKey;
  param.Muon_Trigger_SF_Key = MuonTrigSFKey;
  param.Muon_ISO_SF_Key = MuonISOSFKey;
  param.Electron_Loose_ID = ElectronIDs.at(0);
  param.Jet_ID = "tight";

  executeEventFromParameter(param);

}

void SkimTree_PUJets::executeEventFromParameter(AnalyzerParameter param){
	Event ev = GetEvent();
	if(!(PassMETFilter())) return;
 	Particle met = ev.GetMETVector();
	if(!ev.PassTrigger(IsoMuTriggerName)) return;

	vector<Muon> AllMuons = vec_muon;
	vector<Electron> AllElectrons = vec_electron;
	vector<Jet> AllJets = vec_jet;


	vector<Muon> muons = SelectMuons(AllMuons, param.Muon_Loose_ID, 30., 2.4);
  vector<Electron> electrons = SelectElectrons(AllElectrons, param.Electron_Loose_ID ,30. ,2.5);
	vector<Jet> jets_unVeto = SelectJets(AllJets, param.Jet_ID, 15, 4.7);
  vector<Jet> jets = JetsVetoLeptonInside(jets_unVeto,electrons,muons,0.4);

	std::sort(muons.begin(),muons.end(),PtComparing);
	std::sort(jets.begin(),jets.end(),PtComparing);

	if(muons.size() != 2) return;
  if(jets.size() != 1) return;
  if( muons.at(0).Pt() <= TriggerSafePtCut ) return;

	if(newtree->Fill()<0) exit(EIO);

}

void SkimTree_PUJets::WriteHist(){
  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

  return;
}