#include "Vcb_El.h"

//////////

Vcb_El::Vcb_El()
{
}//Vcb_El::Vcb_El()

//////////

Vcb_El::~Vcb_El()
{
  delete fitter_driver;
}//Vcb_El::~Vcb_El()

//////////

void Vcb_El::initializeAnalyzer()
{
  vec_el_id = {"passTightID"};
  vec_el_id_sf_key = {"ID_SF_passTightID"};

  if(DataYear==2016)
    {
      vec_el_trig.push_back("HLT_Ele27_WPTight_Gsf_v");
      trig_safe_pt_cut = 30.;
    }
  else if(DataYear==2017)
    {
      vec_el_trig.push_back("HLT_Ele35_WPTight_Gsf_v");
      trig_safe_pt_cut = 37.;
    }
  else if(DataYear==2018)
    {
      vec_el_trig.push_back("HLT_Ele32_WPTight_Gsf_v");
      trig_safe_pt_cut = 35.;
    }
  else std::runtime_error("No trigger configuration for year");
  
  //for(auto& trigger_name : vec_mu_trig) cout << "[Vcb_El::initializeAnalyzer] Iso Electron Trigger Name = " << trigger_name << endl;

  //Jet Tagging Parameters
  vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
  //vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::iterativefit, JetTagging::iterativefit));
  mcCorr->SetJetTaggingParameters(vec_jet_tagging_para);

  //Retrieve POG JER
  std::string jetPtResolutionPath = "";
  std::string jetPtResolutionSFPath = "";
  std::string BASE_PATH = std::getenv("SKFlat_WD") + std::string("/data/Run2UltraLegacy_v2/");

  if(DataEra=="2016preVFP")
    {
      jetPtResolutionPath   = BASE_PATH + "2016preVFP/JME/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt";
      jetPtResolutionSFPath = BASE_PATH + "2016preVFP/JME/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt";
    }
  else if(DataEra=="2016postVFP")
    {
      jetPtResolutionPath   = BASE_PATH + "2016postVFP/JME/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt";
      jetPtResolutionSFPath = BASE_PATH + "2016postVFP/JME/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt";
    }
  else if(DataEra=="2017")
    {
      jetPtResolutionPath   = BASE_PATH + "2017/JME/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt";
      jetPtResolutionSFPath = BASE_PATH + "2017/JME/Summer19UL17_JRV2_MC_SF_AK4PFchs.txt";
    }
  else if(DataEra=="2018")
    {
      jetPtResolutionPath   = BASE_PATH + "2018/JME/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt";
      jetPtResolutionSFPath = BASE_PATH + "2018/JME/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt";
    }
  else
    {
      throw std::runtime_error("No JME configuration for year");
    }
  jet_resolution = JME::JetResolution(jetPtResolutionPath);
  jet_resolution_sf = JME::JetResolutionScaleFactor(jetPtResolutionSFPath);

  fitter_driver = new TKinFitterDriver(DataYear);

  run_syst = HasFlag("RunSyst");
  cout << "[Vcb_El::initializeAnalyzer] RunSyst = " << run_syst << endl;

  return;
}//Vcb_El::initializeAnalyzer()

//////////

void Vcb_El::executeEvent()
{
  vec_electron = GetAllElectrons();
  vec_muon = GetAllMuons();
  vec_jet = GetAllJets();

  AnalyzerParameter param;
  
  //electron ids
  for(unsigned int i=0; i<vec_el_id.size(); i++)
    {
      param.Clear();
      
      param.syst_ = AnalyzerParameter::Central;
      
      param.Name = vec_el_id.at(i) + "_" + "Central";

      param.Electron_Tight_ID = vec_el_id.at(i);
      param.Electron_ID_SF_Key = vec_el_id_sf_key.at(i);

      param.Muon_Loose_ID = "POGLoose";

      param.Jet_ID = "tight";

      executeEventFromParameter(param);

      //run syst
      if(run_syst)
	{
	  executeEventFromParameter(param);
	}

    }//loop over electron ids
  
  return;
}//void Vcb_El::executeEvent()

//////////

void Vcb_El::executeEventFromParameter(AnalyzerParameter param)
{
  //no cut
  FillHist(param.Name+"/NoCut_"+param.Name, 0., 1., 1, 0., 1.);

  //met filter
  if(!PassMETFilter()) return;
  FillHist(param.Name+"/MetFilter_"+param.Name, 0., 1., 1, 0., 1.);

  Event ev = GetEvent();
  Particle met = ev.GetMETVector();

  //trigger
  if(!ev.PassTrigger(vec_el_trig)) return;
  FillHist(param.Name+"/Trig_"+param.Name, 0., 1., 1, 0., 1.);

  vector<Muon> vec_this_muon = vec_muon;
  vector<Electron> vec_this_electron = vec_electron;
  vector<Jet> vec_this_jet = vec_jet;

  // //syst basics
  // if(param.syst_ == AnalyzerParameter::Central){}
  // else if(param.syst_ == AnalyzerParameter::JetResUp) vec_this_jet = SmearJets(vec_this_jet, +1);
  // else if(param.syst_ == AnalyzerParameter::JetResDown) vec_this_jet = SmearJets(vec_this_jet, -1);
  // else if(param.syst_ == AnalyzerParameter::JetEnUp) vec_this_jet = SmearJets(vec_this_jet, +1);
  // else if(param.syst_ == AnalyzerParameter::JetEnDown) vec_this_jet = SmearJets(vec_this_jet, -1);
  // else if(param.syst_ == AnalyzerParameter::MuonEnUp) vec_this_muon = ScaleMuons(vec_this_muon, +1);
  // else if(param.syst_ == AnalyzerParameter::MuonEnDown) vec_this_muon = ScaleMuons(vec_this_muon, -1);
  // else if(param.syst_ == AnalyzerParameter::ElectronResUp) return;
  // else if(param.syst_ == AnalyzerParameter::ElectronResDown) return;
  // else if(param.syst_ == AnalyzerParameter::ElectronEnUp) return;
  // else if(param.syst_ == AnalyzerParameter::ElectronEnDown) return;
  // else
  //   {
  //     cerr << "Vcb_Mu::executeEventFromParameter: Wrong syst_" << endl;
  //     exit(1);
  //   }

  vector<Electron> vec_sel_electron = SelectElectrons(vec_this_electron, param.Electron_Tight_ID, ELECTRON_PT, ELECTRON_ETA);
  vector<Muon> vec_sel_muon = SelectMuons(vec_this_muon, param.Muon_Loose_ID, MUON_PT, MUON_ETA);
  
  vector<Jet> vec_sel_jet = SelectJets(vec_this_jet, param.Jet_ID, JET_PT, JET_ETA);
  vec_sel_jet = JetsVetoLeptonInside(vec_sel_jet, vec_sel_electron, vec_sel_muon, DR_LEPTON_VETO);
  
  //sort
  sort(vec_sel_jet.begin(), vec_sel_jet.end(), PtComparing);

  //n of btag
  int nbtag = 0;
  vector<bool> vec_btag;
  for(auto& jet : vec_sel_jet)
    {
      float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score)
        {
          nbtag++;
          vec_btag.push_back(true);
        }
      else vec_btag.push_back(false);
    }

  if(met.Pt()<MET_PT) return;
  if(vec_sel_electron.size()!=1) return;
  if(vec_sel_electron.at(0).Pt()<=trig_safe_pt_cut) return;
  if(vec_sel_muon.size()!=0) return;
  if(vec_sel_jet.size()<4) return;
  if(nbtag<2) return;
  FillHist(param.Name+"/BaselineSelection_"+param.Name, 0., 1., 1, 0., 1.);

  Electron electron = vec_sel_electron.at(0);

  //caculate weight
  float weight = 1.;
  if(!IsData)
    {
      //lumi
      weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");

      //MCweight +1 or -1
      weight *= ev.MCweight();

      //pileup reweight
      weight *= mcCorr->GetPileUpWeight(nPileUp, 0);
      
      //L1 prefire
      weight *= GetPrefireWeight(0);
      
      //SF for electron trigger effi
      weight *= mcCorr->ElectronTrigger_SF(param.Electron_Tight_ID, param.Electron_ID_SF_Key, vec_sel_electron, 0);
      
      cout << "E Trig SF = " << mcCorr->ElectronTrigger_SF(param.Electron_Tight_ID, param.Electron_ID_SF_Key, vec_sel_electron, 0) << endl;
      
    }//if(!IsData)

  return;

  vector<float> vec_resolution_pt;
  for(auto& jet : vec_sel_jet)
    {
      float resolution_pt = jet_resolution.getResolution({{JME::Binning::JetPt, jet.Pt()}, {JME::Binning::JetEta, jet.Eta()}, {JME::Binning::Rho, Rho}});
      float resolution_pt_sf = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jet.Pt()},{JME::Binning::JetEta, jet.Eta()}}, Variation::NOMINAL);

      vec_resolution_pt.push_back(resolution_pt*resolution_pt_sf);
    }

  return;
}//void Vcb_El::executeEventFromParameter(AnalyzerParameter param)

//////////
