#include "PUJets.h"

PUJets::PUJets(){

}

PUJets::~PUJets(){

}

void PUJets::initializeAnalyzer(){
	MuonIDs = {"POGLooseWithMediumRelIso","POGTightWithTightIso"};
  ElectronIDs = {"passLooseID","passTightID"};
	MuonIDSFKeys = {"NUM_TightID_DEN_TrackerMuons"};
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


	vector<Muon> muons = SelectMuons(AllMuons, param.Muon_Loose_ID, 15., 2.4);
  vector<Electron> electrons = SelectElectrons(AllElectrons, param.Electron_Loose_ID ,15. ,2.5);
	vector<Jet> jets_unVeto = SelectJets(AllJets, param.Jet_ID, 15, 4.7);

  //----baseline selection----
  if(!((muons.size() == 2 && electrons.size() == 0) || (muons.size() == 0 && electrons.size() == 2)))  return;

  //----Jet Cleaning----
  vector<Jet> jets = JetsVetoLeptonInside(jets_unVeto,electrons,muons,0.4);

	std::sort(muons.begin(),muons.end(),PtComparing);
	std::sort(jets.begin(),jets.end(),PtComparing);

  //----Tight Muon Selection----

  vector<Muon> muons_tight = SelectMuons(muons, param.Muon_Tight_ID, TriggerSafePtCut, 2.4);
  if(muons_tight.size() != 2) return;
  if((muons_tight.at(0).Charge()+muons_tight.at(1).Charge()) != 0) return;
	Particle ZCand = muons_tight.at(0) + muons_tight.at(1);
  Particle ZCand_reverse = Particle(-1.*ZCand.Px(),-1.*ZCand.Py(),-1.*ZCand.Pz(),ZCand.E());
  if(ZCand.M()<50) return;

 
  bool run_debug = true;
  double weight = 1.;  

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

    for(unsigned int i = 0; i<muons.size();i++){
      double muon_ISO_sf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(),0);
      double muon_ID_sf = mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muons.at(i).Eta(), muons.at(i).Pt(), 0);
      if(run_debug == true) cout << "iso w" << muon_ISO_sf <<endl;
      if(run_debug == true) cout << "ID w" << muon_ID_sf <<endl;
      weight *= muon_ISO_sf*muon_ID_sf;  
    }

    


  }


  FillHist(param.Name+"/ZMass",ZCand.M(),weight,50,70,110);
  
  if(jets.size() != 1) return;	
  double ZPt = ZCand.Pt();
  double JetPt = jets.at(0).Pt();
  double dPhi_1 = ZCand_reverse.DeltaPhi(jets.at(0));
  bool toggle = true;
  if(jets.at(0).Pt() < 10 || jets.at(0).Pt() >= 50) return;
  if((ZPt/JetPt)>1.5) return;
  if((ZPt/JetPt)<0.5) toggle = false;
  
  vector<TString> histIDs = getHistKey(IsDATA, toggle, getBinID(jets.at(0).Eta(),jets.at(0).Pt()), jets.at(0).GetPUID(jets.at(0).Pt(),jets.at(0).Eta()));
  for(unsigned int i = 0; i < histIDs.size(); i++){
    FillHist(histIDs.at(i),dPhi_1,weight,50,-3.15/4,3.15/4);
  }
}

unsigned int PUJets::getBinID(const double eta, const double Pt){
  unsigned PtNum = 0;
  unsigned etaNum = 0;
  if(10.<=Pt && Pt < 20.) PtNum = 1;
  else if(20.<=Pt && Pt < 30.) PtNum = 5;
  else if(30.<=Pt && Pt < 40.) PtNum = 9;
  else if(40.<=Pt && Pt < 50.) PtNum = 13;
  if( fabs(eta) < 2.5) etaNum = 0;
  else if(2.5<=fabs(eta) && fabs(eta) < 2.75) etaNum = 1;
  else if(2.75<=fabs(eta) && fabs(eta) < 3.0) etaNum = 2;
  else if(3.0<=fabs(eta) && fabs(eta) < 5.0) etaNum = 3; 
  unsigned int binID = PtNum + etaNum;
  return binID;
}

std::vector<TString> PUJets::getHistKey(bool isDATA, bool pt_balanced, unsigned int binID, int PUID){
  TString isDATAKey,pt_balancedKey,binIDKey,PUIDKey,passPUIDKey;
  isDATAKey = isDATA ? "isDATA" : "isMC";
  pt_balancedKey = pt_balanced ? "pt_balanced" : "pt_unbalanced";
  binIDKey = to_string(binID);
  vector<TString> result;
  result.push_back(isDATAKey + "/byPUID/" + pt_balancedKey + "/" + binIDKey);
  switch(PUID){
    case 0b111:
      result.push_back(isDATAKey + "/byPassPUID/" + pt_balancedKey + "/passTight"+ "/" + binIDKey);
      result.push_back(isDATAKey + "/byPassPUID/" + pt_balancedKey + "/passMedium"+ "/" + binIDKey);
      result.push_back(isDATAKey + "/byPassPUID/" + pt_balancedKey + "/passLoose"+ "/" + binIDKey);
      break;
    case 0b110:
      result.push_back(isDATAKey + "/byPassPUID/" + pt_balancedKey + "/passMedium"+ "/" + binIDKey);
      result.push_back(isDATAKey + "/byPassPUID/" + pt_balancedKey + "/passLoose"+ "/" + binIDKey);
      break;
    case 0b100:
      result.push_back(isDATAKey + "/byPassPUID/" + pt_balancedKey + "/passLoose"+ "/" + binIDKey);
      break;
    case 0b0:
      result.push_back(isDATAKey + "/byPassPUID/" + pt_balancedKey + "/fail"+ "/" + binIDKey);
      break;
  }
  return result;
}
