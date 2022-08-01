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
  run_debug = HasFlag("RunDebug");
  cout << "[Vcb_El::initializeAnalyzer] RunDebug = " << run_debug << endl;

  run_syst = HasFlag("RunSyst");
  cout << "[Vcb_El::initializeAnalyzer] RunSyst = " << run_syst << endl;

  vec_el_id = {"passTightID"};
  vec_el_id_sf_key = {"ID_SF_passTightID"};

  if(DataYear==2016)
    {
      vec_el_trig.push_back("HLT_Ele27_WPTight_Gsf_v");
      el_trig = "Ele27";
      trig_safe_pt_cut = 30.;
    }
  else if(DataYear==2017)
    {
      vec_el_trig.push_back("HLT_Ele35_WPTight_Gsf_v");
      el_trig = "Ele35";
      trig_safe_pt_cut = 37.;
    }
  else if(DataYear==2018)
    {
      vec_el_trig.push_back("HLT_Ele32_WPTight_Gsf_v");
      el_trig = "Ele32";
      trig_safe_pt_cut = 35.;
    }
  else std::runtime_error("No trigger configuration for year");
  
  for(auto& trigger_name : vec_el_trig) cout << "[Vcb_El::initializeAnalyzer] Electron Trigger Name = " << trigger_name << endl;

  //Jet Tagging Parameters
  if(run_debug) vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
  else vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::iterativefit, JetTagging::iterativefit));
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
      param.Electron_Loose_ID = "passLooseID";
      param.Electron_ID_SF_Key = vec_el_id_sf_key.at(i);

      param.Muon_Tight_ID = "POGTightWithTightIso";

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

  vector<Electron> vec_electron_veto = SelectElectrons(vec_this_electron, param.Electron_Loose_ID, ELECTRON_PT_VETO, ELECTRON_ETA);
  vector<Muon> vec_muon_veto = SelectMuons(vec_this_muon, param.Muon_Tight_ID, MUON_PT_VETO, MUON_ETA);
  
  vector<Electron> vec_sel_electron = SelectElectrons(vec_this_electron, param.Electron_Tight_ID, trig_safe_pt_cut, ELECTRON_ETA);

  vector<Jet> vec_sel_jet = SelectJets(vec_this_jet, param.Jet_ID, JET_PT, JET_ETA);
  vec_sel_jet = JetsVetoLeptonInside(vec_sel_jet, vec_electron_veto, vec_muon_veto, DR_LEPTON_VETO);
  
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
  if(vec_electron_veto.size()!=1) return;
  if(vec_muon_veto.size()!=0) return;
  if(vec_sel_electron.size()!=1) return;
  if(vec_sel_electron.at(0).Pt()<=trig_safe_pt_cut) return;
  if(vec_sel_jet.size()<4) return;
  if(nbtag<2) return;
  FillHist(param.Name+"/BaselineSelection_"+param.Name, 0., 1., 1, 0., 1.);

  Electron electron = vec_sel_electron.at(0);

  //caculate weight
  float weight = 1.;
  if(!IsData)
    {
      //lumi
      //float lumi_weight = weight_norm_1invpb*ev.GetTriggerLumi("Full");
      //weight *= lumi_weight;

      //MCweight +1 or -1
      //float mc_weight = ev.MCweight();
      //weight *= mc_weight;
      
      //pileup reweight
      float pileup_weight = mcCorr->GetPileUpWeight(nPileUp, 0);
      weight *= pileup_weight;
      
      //L1 prefire
      float prefire_weight = GetPrefireWeight(0);
      weight *= prefire_weight;

      //SF for electron trigger effi
      float sf_el_trig_effi =  mcCorr->ElectronTrigger_SF(param.Electron_Tight_ID, el_trig, vec_sel_electron, 0);
      weight *= sf_el_trig_effi;
      
      if(run_debug) cout << "E Trig SF = " << electron.scEta() << "\t" << electron.UncorrPt() << "\t" << sf_el_trig_effi << endl;

      //SF for electron ID eff
      float sf_el_id_effi = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electron.scEta(), electron.UncorrPt(), 0);
      weight *= sf_el_id_effi;
      
      if(run_debug) cout << "E ID SF = " << sf_el_id_effi << endl;
      
      //SF for electron Reco eff
      float sf_el_reco_effi = mcCorr->ElectronReco_SF(electron.scEta(), electron.UncorrPt(), 0); 
      weight *= sf_el_reco_effi;
      
      if(run_debug) cout << "E Reco SF = " << sf_el_reco_effi << endl;
      
      //SF for b-tagging
      float sf_b_tag = 1;
      if(run_debug) sf_b_tag = mcCorr->GetBTaggingReweight_1a(vec_sel_jet, vec_jet_tagging_para.at(0));
      else sf_b_tag = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0));
      weight *= sf_b_tag;
      
      if(run_debug) cout << "B Tagging SF = " << sf_b_tag << endl;

    }//if(!IsData)

  if(nbtag==2)
    {
      FillHist(param.Name+"/TwoB/Met", met.Vect().Mag(), weight, 50, 0, 300);
      FillHist(param.Name+"/TwoB/N_Vertex", nPV, weight, 100, 0, 100);
      FillHist(param.Name+"/TwoB/Electron_Pt", electron.Pt(), weight, 50, 0, 200);
      FillHist(param.Name+"/TwoB/Electron_Eta", electron.Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/TwoB/Leading_Jet_Pt", vec_sel_jet.at(0).Pt(), weight, 50, 0, 300);
      FillHist(param.Name+"/TwoB/Leading_Jet_Eta", vec_sel_jet.at(0).Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/TwoB/Subleading_Jet_Pt", vec_sel_jet.at(1).Pt(), weight, 50, 0, 300);
      FillHist(param.Name+"/TwoB/Subleading_Jet_Eta", vec_sel_jet.at(1).Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/TwoB/N_Jet", vec_sel_jet.size(), weight, 10, 0, 10);
      FillHist(param.Name+"/TwoB/N_BJet", nbtag, weight, 10, 0, 10);
    }
  else
    {
      FillHist(param.Name+"/ThreeB/Met", met.Vect().Mag(), weight, 50, 0, 300);
      FillHist(param.Name+"/ThreeB/N_Vertex", nPV, weight, 100, 0, 100);
      FillHist(param.Name+"/ThreeB/Electron_Pt", electron.Pt(), weight, 50, 0, 200);
      FillHist(param.Name+"/ThreeB/Electron_Eta", electron.Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/ThreeB/Leading_Jet_Pt", vec_sel_jet.at(0).Pt(), weight, 50, 0, 300);
      FillHist(param.Name+"/ThreeB/Leading_Jet_Eta", vec_sel_jet.at(0).Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/ThreeB/Subleading_Jet_Pt", vec_sel_jet.at(1).Pt(), weight, 50, 0, 300);
      FillHist(param.Name+"/ThreeB/Subleading_Jet_Eta", vec_sel_jet.at(1).Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/ThreeB/N_Jet", vec_sel_jet.size(), weight, 10, 0, 10);
      FillHist(param.Name+"/ThreeB/N_BJet", nbtag, weight, 10, 0, 10);
    }
  
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
