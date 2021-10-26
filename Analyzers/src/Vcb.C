#include "Vcb.h"

void Vcb::initializeAnalyzer()
{
  vec_muid = {"POGTight"};
  //vec_muid_sf_key = {""};
  
  if(DataYear==2017)
    {
      iso_mu_trig_name = "HLT_IsoMu27_v";
      trig_safe_pt_cut = 29.;
    }

  vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::comb) );
  mcCorr->SetJetTaggingParameters(jtps);

  //Retrieve POG JER
  std::string jetPtResolutionPath = "";
  std::string jetPtResolutionSFPath = "";
  std::string BASE_PATH = std::getenv("SKFlat_WD") + std::string("/data/Run2UltraLegacy_v2/");
  if(DataYear==2016)
    {
      jetPtResolutionPath   = BASE_PATH + "2016/JME/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
      jetPtResolutionSFPath = BASE_PATH + "2016/JME/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
    }
  else if(DataYear==2017)
    {
      jetPtResolutionPath   = BASE_PATH + "2017/JME/Fall17_V3_MC_PtResolution_AK4PFchs.txt";
      jetPtResolutionSFPath = BASE_PATH + "2017/JME/Fall17_V3_MC_SF_AK4PFchs.txt";
    }
  else if(DataYear==2018)
    {
      jetPtResolutionPath   = BASE_PATH + "2018/JME/Autumn18_V7_MC_PtResolution_AK4PFchs.txt";
      jetPtResolutionSFPath = BASE_PATH + "2018/JME/Autumn18_V7_MC_SF_AK4PFchs.txt";
    }
  else
    {
      throw std::runtime_error("no configuration for year");
    }
  jet_resolution    = JME::JetResolution(jetPtResolutionPath);
  jet_resolution_sf = JME::JetResolutionScaleFactor(jetPtResolutionSFPath);

  fitter_driver = new TKinFitterDriver(DataYear);

  return;
}//void Vcb::initializeAnalyzer()

//////////

void Vcb::executeEvent()
{
  vec_muon = GetAllMuons();
  vec_jet = GetAllJets();

  AnalyzerParameter param;

  for(auto&  muid : vec_muid)
    {
      param.Clear();
      
      param.syst_ = AnalyzerParameter::Central;
      param.Name = muid + "_" + "Central";
      param.Muon_Tight_ID = muid;
      param.Jet_ID = "tight";

      executeEventFromParameter(param);
    }

  return;
}//void Vcb::executeEvent()

//////////

void Vcb::executeEventFromParameter(AnalyzerParameter param)
{
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

  //n n of btag
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

   //baseline selection
  if(vec_sel_muon.size()!=1) return;
  if(vec_sel_jet.size()<4) return;
  if(met.Pt()<MET_PT) return;
  if(nbtag<2) return;
  FillHist(param.Name+"/BaselineSelection_"+param.Name, 0., 1., 1, 0., 1.);
  
  Muon muon = vec_sel_muon.at(0);
  
  vector<float> vec_resolution_pt;
  for(auto& jet : vec_sel_jet)
    {
      float resolution_pt = jet_resolution.getResolution({{JME::Binning::JetPt, jet.Pt()}, {JME::Binning::JetEta, jet.Eta()}, {JME::Binning::Rho, Rho}});
      float resolution_pt_sf = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jet.Pt()},{JME::Binning::JetEta, jet.Eta()}}, Variation::NOMINAL);
      
      vec_resolution_pt.push_back(resolution_pt*resolution_pt_sf);
    }

  //kinematic fitter
  fitter_driver->Set_Objects(vec_sel_jet, vec_resolution_pt, vec_btag, muon, met);
  fitter_driver->Scan();
  
  cout << "Vcb: test" << endl;

  return;
}//void Vcb::executeEventFromParameter(AnalyzerParameter param)

//////////

Vcb::Vcb(){

}//Vcb::Vcb()

Vcb::~Vcb(){
  //delete fitter_driver;
}//Vcb::~Vcb()

//////////

