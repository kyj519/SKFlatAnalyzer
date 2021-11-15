#include "Vcb_Mu.h"

void Vcb_Mu::initializeAnalyzer()
{
  vec_mu_id = {"POGTight"};
  vec_mu_id_sf_key = {"NUM_TightID_DEN_TrackerMuons"};//should be checked

  vec_mu_iso_sf_key = {"NUM_TightRelIso_DEN_TightIDandIPCut"};//should be checked

  iso_mu_trig_name.clear();
  if(DataYear==2016)
    {
      iso_mu_trig_name.push_back("HLT_IsoMu24_v");
      //iso_mu_trig_name.push_back("HLT_IsoTkMu24_v");
      trig_safe_pt_cut = 26.;
    }  
  else if(DataYear==2017)
    {
      iso_mu_trig_name.push_back("HLT_IsoMu27_v");
      trig_safe_pt_cut = 30.;
    }
  else if(DataYear==2018)
    {
      iso_mu_trig_name.push_back("HLT_IsoMu24_v");
      trig_safe_pt_cut = 26.;
    }
  else std::runtime_error("No trigger configuration for year"); 
    
  for(auto& trigger_name : iso_mu_trig_name) cout << "[Vcb_Mu::initializeAnalyzer] Iso Muon Trigger Name = " << trigger_name << endl;
  cout << "[Vcb_Mu::initializeAnalyzer Trigger Safe Pt Cut = " << trig_safe_pt_cut << endl;
  
  //Jet Tagging Parameters
  //vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb));//should be checked
  vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::iterativefit, JetTagging::iterativefit));
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
  cout << "[Vcb_Mu::initializeAnalyzer] RunSyst = " << run_syst << endl;
  
  return;
}//void Vcb_Mu::initializeAnalyzer()

//////////

void Vcb_Mu::executeEvent()
{
  vec_muon = GetAllMuons();
  vec_electron = GetAllElectrons();
  vec_jet = GetAllJets();
  
  AnalyzerParameter param;

  //muon ids
  for(unsigned int i=0; i<vec_mu_id.size(); i++)
    {
      param.Clear();
      
      param.syst_ = AnalyzerParameter::Central;
      
      param.Name = vec_mu_id.at(i) + "_" + "Central";
      
      param.Muon_Tight_ID = vec_mu_id.at(i);
      param.Muon_ID_SF_Key = vec_mu_id_sf_key.at(i);
      param.Muon_ISO_SF_Key = vec_mu_iso_sf_key.at(i);

      param.Electron_Loose_ID = "passLooseID";

      param.Jet_ID = "tight";

      executeEventFromParameter(param);

      //run syst
      if(run_syst)
	{
	  //basic variation defined in AnalyzerParameter
	  for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++)
	    {
	      param.syst_ = AnalyzerParameter::Syst(it_syst);
	      param.Name = vec_mu_id.at(i) + "_" + "Syst_" + param.GetSystType();
	      executeEventFromParameter(param);
	    }

	  //pileup reweight
	  for(int i=0; i<2; i++)
	    {
	      
	    }

	  //prefire reweight
	  for(int i=0; i<2; i++)
	    {
	    }
	  
	}//run syst
    }//loop over muon ids

  return;
}//void Vcb_Mu::executeEvent()

//////////

void Vcb_Mu::executeEventFromParameter(AnalyzerParameter param)
{
  //no cut
  FillHist(param.Name+"/NoCut_"+param.Name, 0., 1., 1, 0., 1.);
  
  //met filter
  if(!PassMETFilter()) return;
  FillHist(param.Name+"/MetFilter_"+param.Name, 0., 1., 1, 0., 1.);
  
  Event ev = GetEvent();
  Particle met = ev.GetMETVector();
  
  //trigger
  if(!ev.PassTrigger(iso_mu_trig_name)) return;
  FillHist(param.Name+"/Trig_"+param.Name, 0., 1., 1, 0., 1.);
  
  vector<Muon> vec_this_muon = vec_muon;
  vector<Electron> vec_this_electron = vec_electron;
  vector<Jet> vec_this_jet = vec_jet;
  
  //syst basics
  if(param.syst_ == AnalyzerParameter::Central){}
  else if(param.syst_ == AnalyzerParameter::JetResUp) vec_this_jet = SmearJets(vec_this_jet, +1);
  else if(param.syst_ == AnalyzerParameter::JetResDown) vec_this_jet = SmearJets(vec_this_jet, -1);
  else if(param.syst_ == AnalyzerParameter::JetEnUp) vec_this_jet = SmearJets(vec_this_jet, +1);
  else if(param.syst_ == AnalyzerParameter::JetEnDown) vec_this_jet = SmearJets(vec_this_jet, -1);
  else if(param.syst_ == AnalyzerParameter::MuonEnUp) vec_this_muon = ScaleMuons(vec_this_muon, +1);
  else if(param.syst_ == AnalyzerParameter::MuonEnDown) vec_this_muon = ScaleMuons(vec_this_muon, -1);
  else if(param.syst_ == AnalyzerParameter::ElectronResUp) return;
  else if(param.syst_ == AnalyzerParameter::ElectronResDown) return;
  else if(param.syst_ == AnalyzerParameter::ElectronEnUp) return;
  else if(param.syst_ == AnalyzerParameter::ElectronEnDown) return;
  else
    {
      cerr << "Vcb_Mu::executeEventFromParameter: Wrong syst_" << endl;
      exit(1);
    }
  
  vector<Muon> vec_sel_muon = SelectMuons(vec_this_muon, param.Muon_Tight_ID, MUON_PT, MUON_ETA);
  vector<Electron> vec_sel_electron = SelectElectrons(vec_this_electron, param.Electron_Loose_ID, ELECTRON_PT, ELECTRON_ETA);
  
  vector<Jet> vec_sel_jet = SelectJets(vec_this_jet, param.Jet_ID, JET_PT, JET_ETA);
  vec_sel_jet = JetsVetoLeptonInside(vec_sel_jet, vec_sel_electron, vec_sel_muon, DR_LEPTON_VETO);
  
  //sort
  sort(vec_sel_muon.begin(), vec_sel_muon.end(), PtComparing);
  sort(vec_sel_jet.begin(), vec_sel_jet.end(), PtComparing);

  //n of btag
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
  if(vec_sel_muon.at(0).Pt()<=trig_safe_pt_cut) return;
  if(vec_sel_electron.size()!=0) return;
  if(vec_sel_jet.size()<4) return;
  if(met.Pt()<MET_PT) return;
  if(nbtag!=2) return;//
  FillHist(param.Name+"/BaselineSelection_"+param.Name, 0., 1., 1, 0., 1.);
  
  Muon muon = vec_sel_muon.at(0);
  
  vector<float> vec_resolution_pt;
  for(auto& jet : vec_sel_jet)
    {
      float resolution_pt = jet_resolution.getResolution({{JME::Binning::JetPt, jet.Pt()}, {JME::Binning::JetEta, jet.Eta()}, {JME::Binning::Rho, Rho}});
      float resolution_pt_sf = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jet.Pt()},{JME::Binning::JetEta, jet.Eta()}}, Variation::NOMINAL);
      
      vec_resolution_pt.push_back(resolution_pt*resolution_pt_sf);
    }

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
     
      //SF for muon id
      weight *= mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muon.Eta(), muon.MiniAODPt());
     
      //SF for muon iso
      weight *= mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muon.Eta(), muon.MiniAODPt());

      //SF for muon trigger effi
      weight *= mcCorr->MuonTrigger_SF(param.Muon_Tight_ID, param.Muon_ID_SF_Key, vec_sel_muon, 0);
     
      //SF for b-tagging
      weight *= mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0));
     
    }//if(!isData)
    
  FillHist(param.Name+"/Met", met.Vect().Mag(), weight, 50, 0, 300);
  FillHist(param.Name+"/Pileup", nPileUp, weight, 100, 0, 100);
  FillHist(param.Name+"/Muon_Pt", muon.Pt(), weight, 50, 0, 200);
  FillHist(param.Name+"/Muon_Eta", muon.Eta(), weight, 50, -5, 5);
  FillHist(param.Name+"/Leading_Jet_Pt", vec_sel_jet.at(0).Pt(), weight, 50, 0, 300);
  FillHist(param.Name+"/Leading_Jet_Eta", vec_sel_jet.at(0).Eta(), weight, 50, -5, 5);
  FillHist(param.Name+"/Subleading_Jet_Pt", vec_sel_jet.at(1).Pt(), weight, 50, 0, 300);
  FillHist(param.Name+"/Subleading_Jet_Eta", vec_sel_jet.at(1).Eta(), weight, 50, -10, 10);
  FillHist(param.Name+"/N_BJet", nbtag, weight, 10, 0, 10);

  return;

  //kinematic fitter
  fitter_driver->Set_Objects(vec_sel_jet, vec_resolution_pt, vec_btag, muon, met);
  fitter_driver->Scan();
  if(fitter_driver->Check_Status()==true)
    {
      Results_Container results_container = fitter_driver->Get_Results();
      
      float chi2 = results_container.best_chi2;

      float initial_had_t_m = results_container.best_initial_had_t_mass;
      float initial_had_w_m = results_container.best_initial_had_w_mass;
      float initial_lep_t_m = results_container.best_initial_lep_t_mass;
      float initial_lep_w_m = results_container.best_initial_lep_w_mass;
      
      float fitted_had_t_m = results_container.best_fitted_had_t_mass;
      float fitted_had_w_m = results_container.best_fitted_had_w_mass;
      float fitted_lep_t_m = results_container.best_fitted_lep_t_mass;
      float fitted_lep_w_m = results_container.best_fitted_lep_w_mass;
	
      FillHist(param.Name+"/Chi2", results_container.best_chi2, weight, 30, 0, 30);

      FillHist(param.Name+"/Had_T_M_Initial", chi2, initial_had_t_m, weight, 50, 0, 500, 50, 0, 350);
      FillHist(param.Name+"/Had_W_M_Initial", chi2, initial_had_w_m, weight, 50, 0, 500, 50, 0, 200);
      FillHist(param.Name+"/Lep_T_M_Initial", chi2, initial_lep_t_m, weight, 50, 0, 500, 50, 0, 350);
      FillHist(param.Name+"/Lep_W_M_Initial", chi2, initial_lep_w_m, weight, 50, 0, 500, 50, 0, 200);
    
      FillHist(param.Name+"/Had_T_M_Fitted", chi2, fitted_had_t_m, weight, 50, 0, 500, 50, 0, 350);
      FillHist(param.Name+"/Had_W_M_Fitted", chi2, fitted_had_w_m, weight, 50, 0, 500, 50, 0, 200);
      FillHist(param.Name+"/Lep_T_M_Fitted", chi2, fitted_lep_t_m, weight, 50, 0, 500, 50, 0, 350);
      FillHist(param.Name+"/Lep_W_M_Fitted", chi2, fitted_lep_w_m, weight, 50, 0, 500, 50, 0, 200);
    
      int index_w_u = results_container.best_index_w_u;
      Jet jet_w_u = vec_sel_jet.at(index_w_u);
      
      float bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
      float cvsb_w_u = Get_CvsB(jet_w_u);
      float cvsl_w_u = Get_CvsL(jet_w_u);
      
      int index_w_d = results_container.best_index_w_d;
      Jet jet_w_d = vec_sel_jet.at(index_w_d);

      float bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
      float cvsb_w_d = Get_CvsB(jet_w_d);
      float cvsl_w_d = Get_CvsL(jet_w_d);

      FillHist(param.Name+"/B_vs_C_W_U", bvsc_w_u, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_B_W_U", cvsb_w_u, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_L_W_U", cvsl_w_u, weight, 10, 0, 1);
    
      FillHist(param.Name+"/B_vs_C_W_D", bvsc_w_d, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_B_W_D", cvsb_w_d, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_L_W_D", cvsl_w_d, weight, 10, 0, 1);

      //n b-tagged
      float nbtag_w = 0;
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < bvsc_w_u) nbtag_w++;
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < bvsc_w_d) nbtag_w++;
      nbtag_w += 0.5;
      
      FillHist(param.Name+"/N_BTag_W_Candidate", nbtag_w, weight, 3, 0, 3);

      //n c-tagged
      
      //naive template
      FillHist(param.Name+"/Naive_Template", bvsc_w_u, cvsb_w_u, cvsl_w_u, weight, 20, 0, 1, 20, 0, 1, 20, 0, 1);
      FillHist(param.Name+"/Naive_Template", bvsc_w_d, cvsb_w_d, cvsl_w_d, weight, 20, 0, 1, 20, 0, 1, 20, 0, 1);
    }//if(fitter_driver->Check_Status()==true)
  
  return;
}//void Vcb_Mu::executeEventFromParameter(AnalyzerParameter param)

//////////

Vcb_Mu::Vcb_Mu()
{
}//Vcb_Mu::Vcb_Mu()

Vcb_Mu::~Vcb_Mu()
{
  delete fitter_driver;
}//Vcb_Mu::~Vcb_Mu()

//////////

float Vcb_Mu::Get_CvsB(const Jet& jet)
{
  float prob_c = jet.GetTaggerResult(JetTagging::DeepFlavour_c);
  float prob_b = jet.GetTaggerResult(JetTagging::DeepFlavour_b);
  float prob_bb = jet.GetTaggerResult(JetTagging::DeepFlavour_bb);
  float prob_lepb = jet.GetTaggerResult(JetTagging::DeepFlavour_lepb); 
  
  float c_vs_b = prob_c/(prob_c+ prob_b+prob_bb+prob_lepb);

  return c_vs_b;
}//float Vcb_Mu::Get_CvsB(const Jet& jet)

//////////

float Vcb_Mu::Get_CvsL(const Jet& jet)
{
  float prob_c = jet.GetTaggerResult(JetTagging::DeepFlavour_c);
  float prob_uds = jet.GetTaggerResult(JetTagging::DeepFlavour_uds);
  float prob_g = jet.GetTaggerResult(JetTagging::DeepFlavour_g);
  
  float c_vs_l = prob_c/(prob_c+prob_uds+prob_g);

  return c_vs_l;
}//float Vcb_Mu::Get_CvsL(const Jet& jet)

//////////
