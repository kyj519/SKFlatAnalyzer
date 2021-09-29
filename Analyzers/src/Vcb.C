#include "Vcb.h"

void Vcb::initializeAnalyzer(){
  vec_muid = {"POGTight"};
  //vec_muid_sf_key = {""};
  
  if(DataYear==2017){
    iso_mu_trig_name = "HLT_IsoMu27_v";
    trig_safe_pt_cut = 29.;
  }

  vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::comb) );
  mcCorr->SetJetTaggingParameters(jtps);
  
  fitter_driver = new TKinFitterDriver(DataYear);

  return;
}//void Vcb::initializeAnalyzer()

//////////

void Vcb::executeEvent(){
  vec_muon = GetAllMuons();
  vec_jet = GetAllJets();

  AnalyzerParameter param;

  for(auto&  muid : vec_muid)
    {
      param.Clear();
      
      param.syst_ = AnalyzerParameter::Central;
      param.Name = muid + "_" + "Central";
      param.Muon_Tight_ID = muid;
      param.Jet_ID = "Tight";

      executeEventFromParameter(param);
    }

  return;
}//void Vcb::executeEvent()

//////////

void Vcb::executeEventFromParameter(AnalyzerParameter param){
  FillHist(param.Name+"/NoCut_"+param.Name, 0., 1., 1, 0., 1.);
  
  if(!PassMETFilter()) return;
  FillHist(param.Name+"/MetFilter_"+param.Name, 0., 1., 1, 0., 1.);

  Event ev = GetEvent();
  Particle met = ev.GetMETVector();
  
  if(!ev.PassTrigger(iso_mu_trig_name)) return;
  FillHist(param.Name+"/Trig_"+param.Name, 0., 1., 1, 0., 1.);
  
  vector<Muon> vec_this_muon = vec_muon;
  vector<Jet> vec_this_jet = vec_jet;

  if(param.syst_ == AnalyzerParameter::Central){}
  else
    {
      cerr << "Vcb::executeEventFromParameter: Wrong syst_" << endl;
      exit(1);
    }
  
  vector<Muon> vec_sel_muon = SelectMuons(vec_this_muon, param.Muon_Tight_ID, MUON_PT, MUON_ETA);
  vector<Jet> vec_sel_jet = SelectJets(vec_this_jet, param.Jet_ID, JET_PT, JET_ETA);

  //sort
  sort(vec_sel_muon.begin(), vec_sel_muon.end(), PtComparing);
  sort(vec_sel_jet.begin(), vec_sel_jet.end(), PtComparing);

  //coutn n of btag
  int nbtag = 0;
  vector<bool> vec_btag;
  for(auto& jet : vec_sel_jet) 
    {
      double tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score)
	{
	  nbtag++;
	  vec_btag.push_back(true);
	}
      else vec_btag.push_back(false);
    }

  Muon muon = vec_sel_muon.at(0);

  //baseline selection
  if(vec_sel_muon.size()!=1) return;
  if(vec_sel_jet.size()<4) return;
  if(met.Pt()<MET_PT) return;
  if(nbtag<2) return;
  FillHist(param.Name+"/BaselineSelection_"+param.Name, 0., 1., 1, 0., 1.);

  //kinematic fitter
  fitter_driver->Set_Objects(vec_jet, vec_btag, muon, met);

  return;
}//void Vcb::executeEventFromParameter(AnalyzerParameter param)

//////////

Vcb::Vcb(){

}//Vcb::Vcb()

Vcb::~Vcb(){
  delete fitter_driver;
}//Vcb::~Vcb()

//////////

