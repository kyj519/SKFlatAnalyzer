#include "PUJets.h"

PUJets::PUJets(){

}

PUJets::~PUJets(){

}

void PUJets::initializeAnalyzer(){
	MuonIDs = {"POGLooseWithRelIso04","POGTight"};
  ElectronIDs = {"passLooseID","passTightID"};
	MuonIDSFKeys = {"NUM_LooseID_DEN_TrackerMuons"};
  MuonISOSFKeys = {"NUM_TightRelIso_DEN_TightIDandIPCut"};
  
  if(DataYear==2016){
    IsoMuTriggerName = "HLT_IsoMu24_v";
    MuonTrigSFKeys = {"IsoMu24"};
    TriggerSafePtCut = 26.;
  }
  else if(DataYear==2017){
    IsoMuTriggerName = "HLT_IsoMu27_v";
    MuonTrigSFKeys = {"IsoMu27"};
    TriggerSafePtCut = 29.;
  }
	else if(DataYear==2018){
    IsoMuTriggerName = "HLT_IsoMu24_v";
    MuonTrigSFKeys = {"IsoMu24"};
    TriggerSafePtCut = 26.;
  }
}

void PUJets::executeEvent(){
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


    if(RunSyst){

      for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){

        //==== Everything else remains same, but only change syst_ and parameter name

        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }

    } 
}

void PUJets::executeEventFromParameter(AnalyzerParameter param){
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
	Particle ZCand = muons.at(0) + muons.at(1);
  Particle ZCand_reverse = Particle(-1.*ZCand.Px(),-1.*ZCand.Py(),-1.*ZCand.Pz(),ZCand.E());
 

	double ZPt = ZCand.Pt();
  double JetPt = jets.at(0).Pt();
  double dPhi_1 = ZCand_reverse.DeltaPhi(jets.at(0));
	double weight = 1.;
  bool run_debug = true;
  
  

  if(!IsData){
      //lumi
    double lumi_weight = weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= lumi_weight;
    if(run_debug == true) cout << "lumiw" << lumi_weight <<endl;

    //MCweight +1 or -1
    double mc_weight = ev.MCweight();
    if(run_debug == true) cout << "mcw" << mc_weight <<endl;
    weight *= mc_weight;


    //pileup reweight
    double pileup_weight = mcCorr->GetPileUpWeight(nPileUp, 0);
    if(run_debug == true) cout << "puw" << pileup_weight <<endl;
    weight *= pileup_weight;

    //L1 prefire
    double prefire_weight = GetPrefireWeight(0);
    if(run_debug == true) cout << "l1W" << prefire_weight <<endl;
    weight *= prefire_weight;

    double muon_trig_sf = mcCorr->MuonTrigger_SF(param.Muon_Tight_ID,param.Muon_Trigger_SF_Key,muons,0);
    if(run_debug == true) cout << "trigW" << muon_trig_sf <<endl;
    weight *= muon_trig_sf;

    for(int i = 0; i<muons.size();i++){
      double muon_ISO_sf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(),0);
      if(run_debug == true) cout << "iso w" << muon_ISO_sf <<endl;
      weight *= muon_ISO_sf;  
    }

    
    


  }


  
  bool trigger = true;
  if(jets.at(0).Pt() < 10 || jets.at(0).Pt() >= 50) return;
  if((ZPt/JetPt)>1.5) return;
  if((ZPt/JetPt)<0.5) trigger = false;
  if(trigger == true){
    FillHist(param.Name+"/dPhi_region1/"+std::to_string(jets.at(0).GetPUID(jets.at(0).Pt(),jets.at(0).Eta())), dPhi_1, weight, 50, -3.15/4., 3.15/4.);
    FillHist(param.Name+"/dPhi_region1/ZMass/"+std::to_string(jets.at(0).GetPUID(jets.at(0).Pt(),jets.at(0).Eta())),ZCand.M(),weight,50,70,110);
  }
  else{
    FillHist(param.Name+"/dPhi_region2/"+std::to_string(jets.at(0).GetPUID(jets.at(0).Pt(),jets.at(0).Eta())), dPhi_1, weight, 50, -3.15/4., 3.15/4.);
    FillHist(param.Name+"/dPhi_region2/ZMass/"+std::to_string(jets.at(0).GetPUID(jets.at(0).Pt(),jets.at(0).Eta())),ZCand.M(),weight,50,70,110);
  }
}