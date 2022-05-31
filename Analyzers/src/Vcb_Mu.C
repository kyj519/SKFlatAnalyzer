#include "Vcb_Mu.h"

//////////

Vcb_Mu::Vcb_Mu()
{
  result_tree = NULL;
  permutation_tree_correct = NULL;
  permutation_tree_wrong = NULL;
}//Vcb_Mu::Vcb_Mu()

//////////

Vcb_Mu::~Vcb_Mu()
{
  delete fitter_driver;

  outfile->cd();
  result_tree->Write();

  if(run_permutation_tree)
    {
      outfile->cd();
      permutation_tree_correct->Write();
      permutation_tree_wrong->Write();
    }
}//Vcb_Mu::~Vcb_Mu()

//////////

void Vcb_Mu::initializeAnalyzer()
{
  run_debug = HasFlag("RunDebug");
  cout << "[Vcb_Mu::initializeAnalyzer] RunDebug = " << run_debug << endl;

  run_permutation_tree = HasFlag("RunPermutationTree");
  cout << "[Vcb_Mu::initializeAnalyzer] RunPermutationTree = " << run_permutation_tree << endl;

  run_kf_eval = HasFlag("RunKFEval");
  cout << "[Vcb_Mu::initializeAnalyzer] RunKFEval = " << run_kf_eval << endl;

  run_syst = HasFlag("RunSyst");
  cout << "[Vcb_Mu::initializeAnalyzer] RunSyst = " << run_syst << endl;

  //to check additional figure of merit using hadronic w mass constrain to suppress background  
  rm_wm_constraint = HasFlag("RM_WM");
  cout << "[Vcb_Mu::initializeAnalyzer] RM_WM = " << rm_wm_constraint << endl;

  vec_mu_id = {"POGTightWithTightIso"};
  vec_mu_id_sf_key = {"NUM_TightID_DEN_TrackerMuons"};
  vec_mu_iso_sf_key = {"NUM_TightRelIso_DEN_TightIDandIPCut"};

  vec_mu_trig.clear();
  if(DataYear==2016)
    {
      vec_mu_trig.push_back("HLT_IsoMu24_v");
      mu_trig = "IsoMu24";
      trig_safe_pt_cut = 26.;
    }  
  else if(DataYear==2017)
    {
      vec_mu_trig.push_back("HLT_IsoMu27_v");
      mu_trig = "IsoMu27";
      trig_safe_pt_cut = 30.;
    }
  else if(DataYear==2018)
    {
      vec_mu_trig.push_back("HLT_IsoMu24_v");
      mu_trig = "IsoMu24";
      trig_safe_pt_cut = 26.;
    }
  else std::runtime_error("No trigger configuration for year"); 
    
  for(auto& trigger_name : vec_mu_trig) cout << "[Vcb_Mu::initializeAnalyzer] Iso Muon Trigger Name = " << trigger_name << endl;
  cout << "[Vcb_Mu::initializeAnalyzer Trigger Safe Pt Cut = " << trig_safe_pt_cut << endl;
  
  //Jet Tagging Parameters
  if(run_debug) vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb));//should be checked
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

  fitter_driver = new TKinFitterDriver(DataYear, rm_wm_constraint);
    
  result_tree = new TTree("Result_Tree", "Result_Tree");
  result_tree->Branch("weight", &weight);
  result_tree->Branch("chi2", &chi2);
  result_tree->Branch("bvsc_w_u", &bvsc_w_u);
  result_tree->Branch("cvsb_w_u", &cvsb_w_u);
  result_tree->Branch("cvsl_w_u", &cvsl_w_u);
  result_tree->Branch("bvsc_w_d", &bvsc_w_d);
  result_tree->Branch("cvsb_w_d", &cvsb_w_d);
  result_tree->Branch("cvsl_w_d", &cvsl_w_d);

  if(run_permutation_tree)
    {
      chk_matched_jets_only = false;
      
      permutation_tree_correct = new TTree("Permutation_Correct", "Permutation_Correct");
      permutation_tree_correct->Branch("gen_neutrino_px", &gen_neutrino_px);
      permutation_tree_correct->Branch("gen_neutrino_py", &gen_neutrino_py);
      permutation_tree_correct->Branch("gen_neutrino_pz", &gen_neutrino_pz);
      permutation_tree_correct->Branch("met_px", &met_px);
      permutation_tree_correct->Branch("met_py", &met_py);
      permutation_tree_correct->Branch("met_rebalance_px", &met_rebalance_px);
      permutation_tree_correct->Branch("met_rebalance_py", &met_rebalance_py);
      permutation_tree_correct->Branch("neutrino_pz_sol", &neutrino_pz_sol);
      permutation_tree_correct->Branch("mt_gen", &mt_gen);
      permutation_tree_correct->Branch("mt_met", &mt_met);
      permutation_tree_correct->Branch("mt_met_rebalance", &mt_met_rebalance);
      permutation_tree_correct->Branch("weight", &weight);
      permutation_tree_correct->Branch("n_jets", &n_sel_jet);
      permutation_tree_correct->Branch("had_t_b_pt", &had_t_b_pt);
      permutation_tree_correct->Branch("w_u_pt", &w_u_pt);
      permutation_tree_correct->Branch("w_d_pt", &w_d_pt);
      permutation_tree_correct->Branch("lep_t_b_pt", &lep_t_b_pt);
      permutation_tree_correct->Branch("had_t_b_bscore", &had_t_b_bscore);
      permutation_tree_correct->Branch("lep_t_b_bscore", &lep_t_b_bscore);
      permutation_tree_correct->Branch("met_rebalance_px", &met_rebalance_px);
      permutation_tree_correct->Branch("met_rebalance_py", &met_rebalance_py);
      permutation_tree_correct->Branch("tt_delta_phi", &tt_delta_phi);
      permutation_tree_correct->Branch("had_t_mass", &had_t_mass);
      permutation_tree_correct->Branch("had_w_mass", &had_w_mass);
      permutation_tree_correct->Branch("lep_t_mass", &lep_t_mass);
      permutation_tree_correct->Branch("lep_t_partial_mass", &lep_t_partial_mass);
      permutation_tree_correct->Branch("chi2_jet_had_t_b", &chi2_jet_had_t_b);
      permutation_tree_correct->Branch("chi2_jet_w_u", &chi2_jet_w_u);
      permutation_tree_correct->Branch("chi2_jet_w_d", &chi2_jet_w_d);
      permutation_tree_correct->Branch("chi2_jet_lep_t_b", &chi2_jet_lep_t_b);
      permutation_tree_correct->Branch("chi2_jet_extra", &chi2_jet_extra);
      permutation_tree_correct->Branch("chi2_constraint_had_t", &chi2_constraint_had_t);
      permutation_tree_correct->Branch("chi2_constraint_had_w", &chi2_constraint_had_w);
      permutation_tree_correct->Branch("chi2_constraint_lep_t", &chi2_constraint_lep_t);
      permutation_tree_correct->Branch("chi2_constraint_lep_w", &chi2_constraint_lep_w);
      permutation_tree_correct->Branch("chi2", &chi2);

      permutation_tree_wrong = new TTree("Permutation_Wrong", "Permutation_Wrong");
      permutation_tree_wrong->Branch("weight", &weight);
      permutation_tree_wrong->Branch("gen_neutrino_px", &gen_neutrino_px);
      permutation_tree_wrong->Branch("gen_neutrino_py", &gen_neutrino_py);
      permutation_tree_wrong->Branch("gen_neutrino_pz", &gen_neutrino_pz);
      permutation_tree_wrong->Branch("met_px", &met_px);
      permutation_tree_wrong->Branch("met_py", &met_py);
      permutation_tree_wrong->Branch("met_rebalance_px", &met_rebalance_px);
      permutation_tree_wrong->Branch("met_rebalance_py", &met_rebalance_py);
      permutation_tree_wrong->Branch("neutrino_pz_sol", &neutrino_pz_sol);
      permutation_tree_wrong->Branch("mt_gen", &mt_gen);
      permutation_tree_wrong->Branch("mt_met", &mt_met);
      permutation_tree_wrong->Branch("mt_met_rebalance", &mt_met_rebalance);
      permutation_tree_wrong->Branch("n_jets", &n_sel_jet);
      permutation_tree_wrong->Branch("had_t_b_pt", &had_t_b_pt);
      permutation_tree_wrong->Branch("w_u_pt", &w_u_pt);
      permutation_tree_wrong->Branch("w_d_pt", &w_d_pt);
      permutation_tree_wrong->Branch("lep_t_b_pt", &lep_t_b_pt);
      permutation_tree_wrong->Branch("had_t_b_bscore", &had_t_b_bscore);
      permutation_tree_wrong->Branch("lep_t_b_bscore", &lep_t_b_bscore);
      permutation_tree_wrong->Branch("met_rebalance_px", &met_rebalance_px);
      permutation_tree_wrong->Branch("met_rebalance_py", &met_rebalance_py);
      permutation_tree_wrong->Branch("tt_delta_phi", &tt_delta_phi);
      permutation_tree_wrong->Branch("had_t_mass", &had_t_mass);
      permutation_tree_wrong->Branch("had_w_mass", &had_w_mass);
      permutation_tree_wrong->Branch("lep_t_mass", &lep_t_mass);
      permutation_tree_wrong->Branch("lep_t_partial_mass", &lep_t_partial_mass);
      permutation_tree_wrong->Branch("chi2_jet_had_t_b", &chi2_jet_had_t_b);
      permutation_tree_wrong->Branch("chi2_jet_w_u", &chi2_jet_w_u);
      permutation_tree_wrong->Branch("chi2_jet_w_d", &chi2_jet_w_d);
      permutation_tree_wrong->Branch("chi2_jet_lep_t_b", &chi2_jet_lep_t_b);
      permutation_tree_wrong->Branch("chi2_jet_extra", &chi2_jet_extra);
      permutation_tree_wrong->Branch("chi2_constraint_had_t", &chi2_constraint_had_t);
      permutation_tree_wrong->Branch("chi2_constraint_had_w", &chi2_constraint_had_w);
      permutation_tree_wrong->Branch("chi2_constraint_lep_t", &chi2_constraint_lep_t);
      permutation_tree_wrong->Branch("chi2_constraint_lep_w", &chi2_constraint_lep_w);
      permutation_tree_wrong->Branch("chi2", &chi2);
    }
  
  if(run_kf_eval) chk_matched_jets_only = true;

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
  Event ev = GetEvent();
  Particle met = ev.GetMETVector();

  //caculate weight for cut flow
  float weight_cut_flow = 1;
  if(!IsData)
    {
      //MCweight +1 or -1
      float mc_weight = ev.MCweight();
      weight_cut_flow *= mc_weight;
    }

  //no cut
  FillHist(param.Name+"/Cut_Flow", Cut_Flow::No_Cut, weight_cut_flow, n_cut_flow, 0, n_cut_flow);
  
  //met filter
  if(!PassMETFilter()) return;
  FillHist(param.Name+"/Cut_Flow", Cut_Flow::Met_Filter, weight_cut_flow, n_cut_flow, 0, n_cut_flow);
  
  //trigger
  if(!ev.PassTrigger(vec_mu_trig)) return;
  FillHist(param.Name+"/Cut_Flow", Cut_Flow::Trigger, weight_cut_flow, n_cut_flow, 0, n_cut_flow);

  vector<Muon> vec_this_muon = vec_muon;
  vector<Electron> vec_this_electron = vec_electron;
  vector<Jet> vec_this_jet = vec_jet;

  //syst basics
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
  vector<Muon> vec_sel_muon = SelectMuons(vec_this_muon, param.Muon_Tight_ID, trig_safe_pt_cut, MUON_ETA);

  vector<Muon> vec_muon_veto = SelectMuons(vec_this_muon, param.Muon_Tight_ID, MUON_PT_VETO, MUON_ETA);
  vector<Electron> vec_electron_veto = SelectElectrons(vec_this_electron, param.Electron_Loose_ID, ELECTRON_PT_VETO, ELECTRON_ETA);
  
  vector<Jet> vec_sel_jet = SelectJets(vec_this_jet, param.Jet_ID, JET_PT, JET_ETA);
  vec_sel_jet = JetsVetoLeptonInside(vec_sel_jet, vec_electron_veto, vec_muon_veto, DR_LEPTON_VETO);
  n_sel_jet = vec_sel_jet.size();

  vector<Jet> vec_sel_jet_match;
  if(!IsData) vec_sel_jet_match = SelectJets(vec_this_jet, param.Jet_ID, 20, JET_ETA);

  //sort jet as pt ordering
  sort(vec_sel_jet.begin(), vec_sel_jet.end(), PtComparing);
  sort(vec_sel_jet_match.begin(), vec_sel_jet_match.end(), PtComparing);
  
  //n of btag
  int n_btag = 0;
  vector<bool> vec_btag;
  for(auto& jet : vec_sel_jet)
    {
      float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score)
        {
          n_btag++;
          vec_btag.push_back(true);
        }
      else vec_btag.push_back(false);
    }

  vector<bool> vec_btag_match;
  for(auto& jet : vec_sel_jet_match)
    {
      float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score) vec_btag_match.push_back(true);
      else vec_btag_match.push_back(false);
    }

  //baseline selection
  if(met.Pt()<MET_PT) return;
  if(vec_muon_veto.size()!=1) return;
  if(vec_electron_veto.size()!=0) return;
  if(vec_sel_muon.size()!=1) return;
  if(n_sel_jet<4) return;
  if(n_btag<2) return;

  FillHist(param.Name+"/Cut_Flow", Cut_Flow::TT_Selection, weight_cut_flow, n_cut_flow, 0, n_cut_flow);

  Muon muon = vec_sel_muon.at(0);

  //caculate weight
  weight = 1.;
  if(!IsData)
    {
      //lumi
      float lumi_weight = weight_norm_1invpb*ev.GetTriggerLumi("Full");
      weight *= lumi_weight;

      //MCweight +1 or -1
      float mc_weight = ev.MCweight();
      weight *= mc_weight;

      //pileup reweight
      float pileup_weight = mcCorr->GetPileUpWeight(nPileUp, 0);
      weight *= pileup_weight;

      //L1 prefire
      float prefire_weight = GetPrefireWeight(0);
      weight *= prefire_weight;

      //SF for muon trigger effi
      float sf_mu_trig_effi = mcCorr->MuonTrigger_SF("POGTight", mu_trig, vec_sel_muon, 0);
      weight *= sf_mu_trig_effi;

      //if(run_debug) cout << "M Trig SF = " << muon.Eta() << "\t" << muon.MiniAODPt() << "\t" << sf_mu_trig_effi << endl;

      //SF for muon id
      float sf_mu_id_effi = mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muon.Eta(), muon.MiniAODPt());
      weight *= sf_mu_id_effi;

      //if(run_debug) cout << "M ID SF = " << muon.Eta() << "\t" << muon.MiniAODPt() << "\t" << sf_mu_trig_effi << endl;

      //SF for muon iso
      float sf_mu_iso_effi = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muon.Eta(), muon.MiniAODPt());
      weight *= sf_mu_iso_effi;

      //SF for b-tagging
      float sf_b_tag = 1;
      if(run_debug) sf_b_tag *= mcCorr->GetBTaggingReweight_1a(vec_sel_jet, vec_jet_tagging_para.at(0));
      else sf_b_tag *= mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0));
      weight *= sf_b_tag;
    }//if(!IsData)

  //Data and MC comparison
  if(n_btag==2)
    {
      FillHist(param.Name+"/TwoB/Met", met.Vect().Mag(), weight, 50, 0, 300);
      FillHist(param.Name+"/TwoB/N_Vertex", nPV, weight, 100, 0, 100);
      FillHist(param.Name+"/TwoB/Muon_Pt", muon.Pt(), weight, 50, 0, 200);
      FillHist(param.Name+"/TwoB/Muon_Eta", muon.Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/TwoB/Leading_Jet_Pt", vec_sel_jet.at(0).Pt(), weight, 50, 0, 300);
      FillHist(param.Name+"/TwoB/Leading_Jet_Eta", vec_sel_jet.at(0).Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/TwoB/Subleading_Jet_Pt", vec_sel_jet.at(1).Pt(), weight, 50, 0, 300);
      FillHist(param.Name+"/TwoB/Subleading_Jet_Eta", vec_sel_jet.at(1).Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/TwoB/N_Jet", vec_sel_jet.size(), weight, 10, 0, 10); 
      FillHist(param.Name+"/TwoB/N_BJet", n_btag, weight, 10, 0, 10);
    }
  else
    {
      FillHist(param.Name+"/ThreeB/Met", met.Vect().Mag(), weight, 50, 0, 300);
      FillHist(param.Name+"/ThreeB/N_Vertex", nPV, weight, 100, 0, 100);
      FillHist(param.Name+"/ThreeB/Muon_Pt", muon.Pt(), weight, 50, 0, 200);
      FillHist(param.Name+"/ThreeB/Muon_Eta", muon.Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/ThreeB/Leading_Jet_Pt", vec_sel_jet.at(0).Pt(), weight, 50, 0, 300);
      FillHist(param.Name+"/ThreeB/Leading_Jet_Eta", vec_sel_jet.at(0).Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/ThreeB/Subleading_Jet_Pt", vec_sel_jet.at(1).Pt(), weight, 50, 0, 300);
      FillHist(param.Name+"/ThreeB/Subleading_Jet_Eta", vec_sel_jet.at(1).Eta(), weight, 60, -3, 3);
      FillHist(param.Name+"/ThreeB/N_Jet", vec_sel_jet.size(), weight, 10, 0, 10);
      FillHist(param.Name+"/ThreeB/N_BJet", n_btag, weight, 10, 0, 10);
    }
  
  vector<int> vec_hf_flavour;
  vector<int> vec_hf_origin;
  int index_matched_jet_match[4];//index of sel_jet_match matched to gen
  int index_matched_jet[4];//index of sel_jet matched to gen
  int index_gen[4];
  int chk_included = -1;
  if(!IsData)
    {
      //scan GenHFHadronMatcher results and construct JER for gen match
      vector<float> vec_jer_match;
      for(auto& jet : vec_sel_jet_match)
	{
	  int jet_flavour = jet.GenHFHadronMatcherFlavour();
	  int jet_origin = jet.GenHFHadronMatcherOrigin();

	  vec_hf_flavour.push_back(jet_flavour);
	  vec_hf_origin.push_back(jet_origin);
	  
	  float resolution_pt = jet_resolution.getResolution({{JME::Binning::JetPt, jet.Pt()}, {JME::Binning::JetEta, jet.Eta()}, {JME::Binning::Rho, Rho}});
	  float resolution_pt_sf = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jet.Pt()},{JME::Binning::JetEta, jet.Eta()}}, Variation::NOMINAL);

	  vec_jer_match.push_back(resolution_pt*resolution_pt_sf);
	}

      //try to find jets from W using closed delta R matching 
      vector<Gen> vec_gen = GetGens();
      
      bool surely_matched[4];
      float matched_jet_dr[4];
      Gen_Match_TT(vec_sel_jet_match, vec_gen, vec_hf_flavour, vec_hf_origin, vec_jer_match, index_gen, index_matched_jet_match, surely_matched, matched_jet_dr);
      
      //PrintGen(vec_gen);
      
      Index_Converter(vec_sel_jet, vec_sel_jet_match, index_matched_jet_match, index_matched_jet);
      
      cout << "n_size sel_jet = " << vec_sel_jet.size() << ", Index of matched_sel_jet[4]" << endl;
      for(int i=0; i<4; i++)
	{
	  cout << index_matched_jet[i] << " ";
	}
      cout << endl;

      //For the performance evaluation of KF
      //only nicely matched matched and well selected (= every matched jets passed baseline selection)
      if(run_kf_eval)
	{
	  cout << "test run_kf_eval" << endl;

	  //reject badly selected event
          if(Chk_Included(index_matched_jet)!=0) return;

	  for(int i=0; i<4; i++)
	    {
	      cout << i << ", " << index_matched_jet_match[i] << ", " << vec_hf_flavour[index_matched_jet_match[i]] << ", " << vec_btag_match[index_matched_jet_match[i]] << endl;

	      //reject unmatched events
	      if(index_matched_jet_match[i]==-999) return;

	      //reject events with too big (>0.4) DR
	      if(surely_matched[i]==false && 0.4<matched_jet_dr[i]) return;
	
	      //require btag for b jet
	      if(vec_hf_flavour[index_matched_jet_match[i]]==5 && vec_btag_match[index_matched_jet_match[i]]==false) return;

	      //require no btag for c jet
	      if(vec_hf_flavour[index_matched_jet_match[i]]==4 && vec_btag_match[index_matched_jet_match[i]]==true) return;
	    }
	}//if(run_kf_eval)
      
      //if all four jets are successfully match
      if(index_matched_jet_match[0]!=-999 && index_matched_jet_match[1]!=-999 && index_matched_jet_match[2]!=-999 && index_matched_jet_match[3]!=-999)
	{
	  FillHist(param.Name+"/PF_Gen_Matched", 1, weight, 2, 0, 2);
 
	  TLorentzVector w_gen_matched = vec_sel_jet_match.at(index_matched_jet_match[1]) + vec_sel_jet_match.at(index_matched_jet_match[2]);
	  TLorentzVector t_gen_matched = vec_sel_jet_match.at(index_matched_jet_match[0]) + w_gen_matched;

	  FillHist(param.Name+"/W_Gen_Matched_Mass", w_gen_matched.M(), weight, 100, 0, 400);
	  FillHist(param.Name+"/T_Gen_Matched_Mass", t_gen_matched.M(), weight, 100, 0, 600);

	  //correction
	  Jet j0 = vec_sel_jet_match.at(index_matched_jet_match[0]);
	  Jet j1 = vec_sel_jet_match.at(index_matched_jet_match[1]);
          Jet j2 = vec_sel_jet_match.at(index_matched_jet_match[2]);

	  float j0_corr = j0.GetBJetRegressionNN("BCorr");
	  j0.SetPtEtaPhiM(j0_corr*j0.Pt(), j0.Eta(), j0.Phi(), j0.M());
	  
	  float j1_corr;
	  if(vec_hf_flavour[index_matched_jet_match[1]]==4) j1_corr = j1.GetBJetRegressionNN("CCorr");
	  else if(vec_hf_flavour[index_matched_jet_match[1]]==5) j1_corr = j1.GetBJetRegressionNN("BCorr");
	  
	  j1.SetPtEtaPhiM(j1_corr*j1.Pt(), j1.Eta(), j1.Phi(), j1.M());

	  float j2_corr;
	  if(vec_hf_flavour[index_matched_jet_match[2]]==4) j2_corr = j2.GetBJetRegressionNN("CCorr");
          else if(vec_hf_flavour[index_matched_jet_match[2]]==5) j2_corr = j2.GetBJetRegressionNN("BCorr");
	  
	  j2.SetPtEtaPhiM(j2_corr*j2.Pt(), j2.Eta(), j2.Phi(), j2.M());

	  TLorentzVector w_gen_matched_corr = j1 + j2;
          TLorentzVector t_gen_matched_corr = j0 + w_gen_matched_corr;

	  FillHist(param.Name+"/W_Gen_Matched_Mass_CorrNN", w_gen_matched_corr.M(), weight, 100, 0, 400);
          FillHist(param.Name+"/T_Gen_Matched_Mass_CorrNN", t_gen_matched_corr.M(), weight, 100, 0, 600);

	  for(int i=0; i<4; i++)
            {
              if(surely_matched[i]==true) FillHist(param.Name+Form("/DR_Surely_Matched_%d", i), matched_jet_dr[i], weight, 30, 0, 3);
            }

	  if(matched_jet_dr[1]<matched_jet_dr[2])
	    {
	      FillHist(param.Name+"/W_DR_Small", matched_jet_dr[1], weight, 30, 0, 3);
	      FillHist(param.Name+"/W_DR_Large", matched_jet_dr[2], weight, 30, 0, 3);
	    }
	  else
	    {
	      FillHist(param.Name+"/W_DR_Small", matched_jet_dr[2], weight, 30, 0, 3);
	      FillHist(param.Name+"/W_DR_Large", matched_jet_dr[1], weight, 30, 0, 3);
	    }
	  
	  //if all matched four jets from tt system are included in vec_jet_sel, chk_included=true
	  chk_included = Chk_Included(index_matched_jet);
	}
      else 
	{
	  FillHist(param.Name+"/PF_Gen_Matched", 0, weight, 2, 0, 2);
	}
      
      //-1 gen matching failed
      // 0 gen matching succeeded && four jets are in selection
      //0< gen matching succeeded &&
      //first bjet is out of selection -> 1 bit
      //first wjet is out of selection -> 2 bit
      //second wjet is out of selection -> 3 bit
      //second bjet is out of selection -> 4 bit        
      FillHist(param.Name+"/Objects_In_Sel_Jet", chk_included, weight, 18, -2, 16);
    }//if(!IsData)

  //kinematic fitter
  vector<float> vec_resolution_pt;
  for(auto& jet : vec_sel_jet)
    {
      //JEC
      float resolution_pt = jet_resolution.getResolution({{JME::Binning::JetPt, jet.Pt()}, {JME::Binning::JetEta, jet.Eta()}, {JME::Binning::Rho, Rho}});
      float resolution_pt_sf = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jet.Pt()},{JME::Binning::JetEta, jet.Eta()}}, Variation::NOMINAL);

      vec_resolution_pt.push_back(resolution_pt*resolution_pt_sf);
    }

  fitter_driver->Set_Objects(vec_sel_jet, vec_resolution_pt, vec_btag, muon, met, chk_matched_jets_only, index_matched_jet);
  fitter_driver->Scan();

  if(!IsData && run_permutation_tree)
    {
      vector<Gen> vec_gen = GetGens();
      Gen gen_neutrino = Neutrino(vec_gen);
      
      gen_neutrino_px = gen_neutrino.Px();
      gen_neutrino_py = gen_neutrino.Py();
      gen_neutrino_pz = gen_neutrino.Pz();

      met_px = met.Px();
      met_py = met.Py();

      Results_Container results_container = fitter_driver->Get_Results();
    
      //search correct index first
      int index_correct = -999;
      float diff_pz = 99999.;
      for(unsigned int i=0; i<results_container.vec_results.size(); ++i)
       	{
       	  Results results = results_container.vec_results[i];

	  int index_had_t_b = results.index_had_t_b;
       	  int index_w_u = results.index_w_u;
       	  int index_w_d = results.index_w_d;
       	  int index_lep_t_b = results.index_lep_t_b;
	  
      	  if(index_matched_jet[0]==index_had_t_b && index_matched_jet[1]==index_w_u && 
      	     index_matched_jet[2]==index_w_d && index_matched_jet[3]==index_lep_t_b)
	    {
	      float diff_current = TMath::Abs(results.neutrino_pz_sol - gen_neutrino_pz);
	      
	      if(diff_current<diff_pz) 
		{
		  index_correct = i;
		  diff_pz = diff_current;
		}
	    }
	}//search done
      
      for(unsigned int i=0; i<results_container.vec_results.size(); ++i)
	{
	  Results results = results_container.vec_results[i];
	  
	  met_rebalance_px = results.met_rebalance_px;
	  met_rebalance_py = results.met_rebalance_py;
	  neutrino_pz_sol = results.neutrino_pz_sol;

	  mt_gen = Calculate_Mt(muon, gen_neutrino_px, gen_neutrino_py); 
	  mt_met = Calculate_Mt(muon, met_px, met_py);
	  mt_met_rebalance = Calculate_Mt(muon, met_rebalance_px, met_rebalance_py);

          int index_had_t_b = results.index_had_t_b;
          int index_w_u = results.index_w_u;
          int index_w_d = results.index_w_d;
          int index_lep_t_b = results.index_lep_t_b;

       	  had_t_b_pt = vec_sel_jet[index_had_t_b].Pt();
	  w_u_pt = vec_sel_jet[index_w_u].Pt();
       	  w_d_pt = vec_sel_jet[index_w_d].Pt();
       	  lep_t_b_pt = vec_sel_jet[index_lep_t_b].Pt();

       	  had_t_b_bscore = vec_sel_jet[index_had_t_b].GetTaggerResult(JetTagging::DeepJet);
       	  lep_t_b_bscore = vec_sel_jet[index_lep_t_b].GetTaggerResult(JetTagging::DeepJet);

       	  tt_delta_phi = results.tt_delta_phi;
	  
       	  had_t_mass = results.initial_had_t_mass;
       	  had_w_mass = results.initial_had_w_mass;
      	  lep_t_mass = results.initial_lep_t_mass;
	  
      	  TLorentzVector lep_t_partial = vec_sel_jet[index_lep_t_b] + muon;
      	  lep_t_partial_mass = lep_t_partial.M();

      	  chi2_jet_had_t_b = results.chi2_jet_had_t_b;
      	  chi2_jet_w_u = results.chi2_jet_w_u;
      	  chi2_jet_w_d = results.chi2_jet_w_d;
      	  chi2_jet_lep_t_b = results.chi2_jet_lep_t_b;
      	  chi2_jet_extra = results.chi2_jet_extra;

      	  chi2_constraint_had_t = results.chi2_constraint_had_t;
      	  chi2_constraint_had_w = results.chi2_constraint_had_w;
      	  chi2_constraint_lep_t = results.chi2_constraint_lep_t;
          chi2_constraint_lep_w = results.chi2_constraint_lep_w;

      	  chi2 = results.chi2;
	  
       	  if(i==index_correct) permutation_tree_correct->Fill();
       	  else permutation_tree_wrong->Fill();
	  
	}//for(unsigned int i=0; i<result_container.vec_results.size(); ++i)
      
      return;
    }//if(run_permutation)
    
  bool chk_fitter_status = fitter_driver->Check_Status();
  bool fitter_correct = false;
  if(chk_fitter_status)
    {
      FillHist(param.Name+"/Cut_Flow", Cut_Flow::KF_Pass, weight_cut_flow, n_cut_flow, 0, n_cut_flow);
      
      FillHist(param.Name+"/Fitter_Succeeded", 0., weight, 1, 0., 1.);
      FillHist(param.Name+"/Fitter_Succeeded_NW", 0., 1, 1, 0., 1.);

      Results_Container results_container = fitter_driver->Get_Results();
      
      chi2 = results_container.best_chi2;

      int index_had_t_b = results_container.best_index_had_t_b;
      int index_w_u = results_container.best_index_w_u;
      int index_w_d = results_container.best_index_w_d;
      int index_lep_t_b = results_container.best_index_lep_t_b;
      
      Jet jet_had_t_b = vec_sel_jet.at(index_had_t_b);
      Jet jet_w_u = vec_sel_jet.at(index_w_u);
      Jet jet_w_d = vec_sel_jet.at(index_w_d);
      Jet jet_lep_t_b = vec_sel_jet.at(index_lep_t_b);
      
      float initial_had_t_m = results_container.best_initial_had_t_mass;
      float initial_had_w_m = results_container.best_initial_had_w_mass;
      float initial_lep_t_m = results_container.best_initial_lep_t_mass;
      float initial_lep_w_m = results_container.best_initial_lep_w_mass;

      float fitted_had_t_m = results_container.best_fitted_had_t_mass;
      float fitted_had_w_m = results_container.best_fitted_had_w_mass;
      float fitted_lep_t_m = results_container.best_fitted_lep_t_mass;
      float fitted_lep_w_m = results_container.best_fitted_lep_w_mass;
      
      FillHist(param.Name+"/Chi2", results_container.best_chi2, weight, 30, 0, 30);

      FillHist(param.Name+"/Had_T_M_Initial", chi2, initial_had_t_m, weight, 500, 0, 250, 50, 0, 350);
      FillHist(param.Name+"/Had_W_M_Initial", chi2, initial_had_w_m, weight, 500, 0, 250, 50, 0, 200);
      FillHist(param.Name+"/Lep_T_M_Initial", chi2, initial_lep_t_m, weight, 500, 0, 250, 50, 0, 350);
      FillHist(param.Name+"/Lep_W_M_Initial", chi2, initial_lep_w_m, weight, 500, 0, 250, 50, 0, 200);

      FillHist(param.Name+"/Had_T_M_Fitted", chi2, fitted_had_t_m, weight, 500, 0, 250, 50, 0, 350);
      FillHist(param.Name+"/Had_W_M_Fitted", chi2, fitted_had_w_m, weight, 500, 0, 250, 50, 0, 200);
      FillHist(param.Name+"/Lep_T_M_Fitted", chi2, fitted_lep_t_m, weight, 500, 0, 250, 50, 0, 350);
      FillHist(param.Name+"/Lep_W_M_Fitted", chi2, fitted_lep_w_m, weight, 500, 0, 250, 50, 0, 200);

      if(!IsData)
	{
	  //check gen-jet matching and KF give same result or not
	  if(index_matched_jet_match[1]!=-999 && index_matched_jet_match[2]!=-999)
	    {
	      Jet jet_gen_matched[2] = {vec_sel_jet_match.at(index_matched_jet_match[1]), vec_sel_jet_match.at(index_matched_jet_match[2])};
	      Jet jet_kf_matched[2] = {vec_sel_jet.at(index_w_u), vec_sel_jet.at(index_w_d)};
	      fitter_correct = Compare_Jet_Pair(jet_gen_matched, jet_kf_matched);
	    
	      //event counter for KF performance evaluation
	      if(run_kf_eval)
		{
		  if(fitter_correct) FillHist(param.Name+"/KF_Eval_Correct", n_sel_jet, n_btag, 1, 10, 0, 10, 8, 0, 8);
		  else FillHist(param.Name+"/KF_Eval_Wrong", n_sel_jet, n_btag, 1, 10, 0, 10, 8, 0, 8);
		}
	    }     
	  else fitter_correct = false;
	    
	  if(run_debug)
	    {
	      cout << "Fitter converged: " << index_had_t_b << ",  " << index_w_u << ", " << index_w_d << ", " << index_lep_t_b << endl;
	      for(unsigned int i=0; i<vec_sel_jet.size(); i++)
		{
		  Jet jet = vec_sel_jet.at(i);
		  
		  int jet_flavour = jet.GenHFHadronMatcherFlavour();
		  int jet_origin = jet.GenHFHadronMatcherOrigin();
		  
		  float pt = jet.Pt();

		  cout << i << ", " << jet_flavour << ", " << jet_origin << ", " << pt << endl;   
		}

	      cout << "Fitter Correct = " << fitter_correct << endl;
	      cout << "Included = " << chk_included << endl;
 	    }
	  
	  //Histo for HF contamination
	  int hf_contamination = 0;
	  for(int i=0; i<2; i++)
	    {
	      Jet jet;
	      if(i==0) jet = vec_sel_jet.at(index_w_u);
	      else jet = vec_sel_jet.at(index_w_d);

	      int jet_flavour = jet.GenHFHadronMatcherFlavour();
	      int jet_origin = jet.GenHFHadronMatcherOrigin();
	      
	      int cases = 0;
	      if(Abs(jet_origin)==6) cases = 1;
	      else if(Abs(jet_origin)==21) cases = 2;
	      
	      hf_contamination += (2*i+1)*cases;
	    }
	  FillHist(param.Name+"/HF_Contamination", hf_contamination, weight, 9, 0, 9);

	  int fitter_score;
	  if(fitter_correct) fitter_score = 0;
	  else
	    {
	      if(chk_included==0) fitter_score = 1;
	      else fitter_score = 2;
	    }
	  FillHist(param.Name+"/Fitter_Correct", fitter_score, chi2, weight, 6, 0, 6, 100, 0, 20);
	}
      
      /* Template */
      bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_u = Get_CvsB(jet_w_u);
      cvsl_w_u = Get_CvsL(jet_w_u);
      
      bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_d = Get_CvsB(jet_w_d);
      cvsl_w_d = Get_CvsL(jet_w_d);

      result_tree->Fill();

      FillHist(param.Name+"/B_vs_C_W_U", bvsc_w_u, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_B_W_U", cvsb_w_u, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_L_W_U", cvsl_w_u, weight, 10, 0, 1);
    
      FillHist(param.Name+"/B_vs_C_W_D", bvsc_w_d, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_B_W_D", cvsb_w_d, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_L_W_D", cvsl_w_d, weight, 10, 0, 1);

      //n b-tagged
      float n_btag_w = 0;
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < bvsc_w_u) n_btag_w++;
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < bvsc_w_d) n_btag_w++;
            
      FillHist(param.Name+"/N_BTag_W_Candidate", n_btag_w, weight, 3, 0, 3);

      //n c-tagged
      
      //naive template
      FillHist(param.Name+"/Naive_Template", bvsc_w_u, cvsb_w_u, cvsl_w_u, weight, 20, 0, 1, 20, 0, 1, 20, 0, 1);
      FillHist(param.Name+"/Naive_Template", bvsc_w_d, cvsb_w_d, cvsl_w_d, weight, 20, 0, 1, 20, 0, 1, 20, 0, 1);
    }//if(fitter_driver->Check_Status()==true)
  else
    {
      FillHist(param.Name+"/Fitter_Diverged", 0., weight, 1, 0., 1.);
      FillHist(param.Name+"/Fitter_Diverged_NW", 0., 1, 1, 0., 1.);

      if(run_debug)
	{
	  cout << "Fitter failed" << endl;
	  cout << "Included = " << chk_included << endl;
	}
      
      int fitter_score;
      if(chk_included==0) fitter_score = 3;
      else fitter_score = 4;
      chi2 = 999;
      FillHist(param.Name+"/Fitter_Correct", fitter_score, chi2, weight, 6, 0, 6, 100, 0, 20);
    }
  
  return;
}//Vcb_Mu::executeEventFromParameter(AnalyzerParameter param)

//////////

float Vcb_Mu::Calculate_Mt(const Particle& lepton, const float& neu_px, const float& neu_py)
{
  float px = lepton.Px() + neu_px;
  float py = lepton.Py() + neu_py;
  
  float mt = TMath::Sqrt(W_MASS*W_MASS + px*px + py*py);

  return mt;
}//float Vcb_Mu::Calculate_Mt(const Particle& lepton, const Particle& met)

//////////

int Vcb_Mu::Chk_Included(const int index_matched_jet[4])
{
  int result = 0;
  unsigned int tmp = 1;
  for(int i=0; i<4; ++i)
    {     
      if(index_matched_jet[i]==-999) result += tmp;

      tmp = tmp << 1;
    }//for(int i=0; i<4; i++)

  return result;
}//int Vcb_Mu::Chk_Included(const int index_matched_jet[4])

//////////

bool Vcb_Mu::Compare_Jet_Pair(const Jet jet0[2], const Jet jet1[2])
{
  Double_t pt0[2] = {jet0[0].Pt(), jet0[1].Pt()};
  // float eta0[2] = {jet0[0].Eta(), jet0[1].Eta()};
  // float phi0[2] = {jet0[0].Phi(), jet0[1].Phi()};

  Double_t pt1[2] = {jet1[0].Pt(), jet1[1].Pt()};
  // float eta1[2] = {jet1[0].Eta(), jet1[1].Eta()};
  // float phi1[2] = {jet1[0].Phi(), jet1[1].Phi()};

  if(Abs(pt0[0]-pt1[0])<1e-8 && Abs(pt0[1]-pt1[1])<1e-8) return true;
  else if(Abs(pt0[0]-pt1[1])<1e-8 && Abs(pt0[1]-pt1[0])<1e-8) return true; 
  
  return false;
}//bool Vcb_Mu::Compare_Jet_Pair(const Jet& jet0, const Jet& jet1)

//////////

void Vcb_Mu::Gen_Match_TT(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, const vector<float>& vec_jer, int index_gen[4], int index_matched_jet[4], bool surely_matched[4], float dr_return[4])
{
  for(int i=0; i<4; i++)
    {
      index_gen[i] = -999;
      index_matched_jet[i] = -999;
      surely_matched[i] = false;
      dr_return[i] = -999;
    }
  
  int selected_w = Gen_Match_W(vec_jet, vec_gen, vec_hf_flavour, vec_hf_origin, vec_jer, &index_gen[1], &index_matched_jet[1], &surely_matched[1], &dr_return[1]);
  
  int index_last_t = -999;
  int index_last_at = -999;
  int index_first_b = -999;
  int index_first_ab = -999;
  for(unsigned int i=0; i<vec_gen.size(); i++)
    {
      Gen gen = vec_gen.at(i);

      int pid = gen.PID();
      int m_index = gen.MotherIndex();

      //find last index of t and tbar
      if(pid==6) index_last_t = i;
      if(pid==-6) index_last_at = i;
	
      //find b from t decay
      if(m_index==index_last_t && pid==5) index_first_b = i;
      if(m_index==index_last_at && pid==-5) index_first_ab = i;
    }
 
  //W+ decays hadronically
  if(selected_w==1)
    {
      index_gen[0] = index_first_b;
      index_gen[3] = index_first_ab;
    }
  //W- decays hadronically
  else if(selected_w==-1)
    {
      index_gen[0] = index_first_ab;
      index_gen[3] = index_first_b;
    }
  //no W decays hadronically
  else 
    {
      if(run_debug) cout << "Can't find hadronically decayed W" << endl;

      return;
    }

  for(unsigned int i=0; i<vec_jet.size(); i++)
    {
      int jet_flavour = vec_hf_flavour.at(i);
      int jet_origin = vec_hf_origin.at(i);
      
      float jet_eta = vec_jet.at(i).Eta();
      float jet_phi = vec_jet.at(i).Phi();
      
      //b from t
      if(jet_flavour==5 && jet_origin==6)
	{
	  Gen gen_b = vec_gen.at(index_first_b);

	  float del_eta = jet_eta - gen_b.Eta();
	  float del_phi = jet_phi - gen_b.Phi();

	  float dr = Sqrt(Power(del_eta, 2.0) + Power(del_phi, 2.0));

	  if(selected_w==1)
	    {
	      index_matched_jet[0] = i;
	      surely_matched[0]  = true;
	      dr_return[0] = dr;
	    }
	  else
	    {
	      index_matched_jet[3] = i;
              surely_matched[3]  = true;
              dr_return[3] = dr;
	    }
	}
      
      //ab from at
      if(jet_flavour==5 && jet_origin==-6)
	{
	  Gen gen_ab = vec_gen.at(index_first_ab);

	  float del_eta = jet_eta - gen_ab.Eta();
          float del_phi = jet_phi - gen_ab.Phi();

	  float dr = Sqrt(Power(del_eta, 2.0) + Power(del_phi, 2.0));

	  if(selected_w==1)
            {
	      index_matched_jet[3] = i;
              surely_matched[3]  = true;
              dr_return[3] = dr;
            }
          else
            {
	      index_matched_jet[0] = i;
              surely_matched[0]  = true;
              dr_return[0] = dr;
            }
	}
    }//for(unsigned int i=0; i<vec_jet.size(); i++)


  if(run_debug)
    {
      cout << endl;
      cout << endl;
      if(selected_w == 1) cout << "W+ decayed hadronically" << endl;
      else cout << "W- decayed hadronically" << endl;
      
      for(unsigned int i=0; i<4; i++)
	{
	  if(index_gen[i]!=-999)
	    {
	      Gen gen = vec_gen.at(index_gen[i]);
	      int pid = gen.PID();
	      
	      Gen gen_mother = vec_gen.at(gen.MotherIndex());
	      int m_pid = gen_mother.PID();

	      float gen_eta = gen.Eta();
	      float  gen_phi = gen.Phi();
	      float gen_pt = gen.Pt();
	      
	      cout << "(" << index_gen[i] << ", " << pid << ", " << m_pid << ", " << gen_eta << ", " << gen_phi << ", " << gen_pt << "), ";
	      
	    }
	  else 
	    {
	      cout << "(), ";
	    }
	}
      cout << endl;
           
      for(unsigned int i=0; i<vec_jet.size(); i++)
	{
	  int jet_flavour = vec_hf_flavour.at(i);
          int jet_origin = vec_hf_origin.at(i);
	  
	  Jet jet = vec_jet.at(i);
	  
	  float jet_eta = jet.Eta();
	  float jet_phi = jet.Phi();
	  float jet_pt =  jet.Pt();
	  
	  float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
	  bool b_tagged = false;
	  if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score) b_tagged = true;
	  
	  cout << "(" << jet_flavour << ", " << jet_origin << ", " << jet_eta << ", " << jet_phi << ", " << jet_pt << ", " << b_tagged << ")" << endl;
	}
      cout << endl;
      
      for(int i=0; i<4; i++)
	{
	  cout << "Matched jet index = " << index_matched_jet[i] << ", surely matched = " << surely_matched[i] << ", d_r = " << dr_return[i] << endl;
	}
    }
  
  return;
}//void Vcb_Mu::Gen_Match_TT(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, int index_matched_jet[4], bool surely_matched[4], float dr_return[4])

//////////

int Vcb_Mu::Gen_Match_W(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, const vector<float>& vec_jer, int index_gen[2], int index_matched_jet[2], bool surely_matched[2], float dr_return[2])
{
  int index_last_w = -999;
  int index_last_aw = -999; 
  int index_d0_w = -999;
  int index_d1_w = -999;
  int index_d0_aw = -999;
  int index_d1_aw = -999;
  
  int selected_w = -999;

  //scan gen to find W
  for(unsigned int i=0; i<vec_gen.size(); i++)
    {
      Gen gen = vec_gen.at(i);

      int pid = gen.PID();
      int m_index = gen.MotherIndex();

      //find last index of W+ and W-
      if(pid==24) index_last_w = i;
      if(pid==-24) index_last_aw = i;
      
      //find decay products of W+
      if(m_index == index_last_w && pid!=24)
	{
	  if(index_d0_w==-999)  index_d0_w = i;
	  else index_d1_w = i;
	}
      
      //find decay products of W-
      if(m_index == index_last_aw && pid!=-24)
	{
	  if(index_d0_aw==-999)  index_d0_aw = i;
          else index_d1_aw = i;
	}
    }
  
  //if input sample is not TT, so both of W couldn't be found 
  if(index_last_w==-999 || index_last_aw==-999)
    {
      if(run_debug) cout << "Can't find both of W" << endl;
      return -999;
    }
  
  int index_had_w[2];
  //W+ decay hadronically
  if(Abs(vec_gen.at(index_d0_w).PID())<10)
    {
      index_gen[0] = index_d0_w;
      index_gen[1] = index_d1_w;

      index_had_w[0] = index_d0_w;
      index_had_w[1] = index_d1_w;
      
      selected_w = 1;//mean W+ decayed hadronically
    }
  //W- decay hadronically
  else if(Abs(vec_gen.at(index_d0_aw).PID()<10))
    {
      index_gen[0] = index_d0_aw;
      index_gen[1] = index_d1_aw;
      
      index_had_w[0] = index_d0_aw;
      index_had_w[1] = index_d1_aw;
   
      selected_w = -1;//mean W- decayed hadronically
    }
  //no W decays hadronically
  else 
    {
      if(run_debug) cout << "No W decays hadronically" << endl;

      return -999;
    }
  
  //calculate del_r for W
  vector<float> vec_dr[2];
  vector<bool> vec_chk_pt_diff[2];
  for(int i=0; i<2; i++)
    {
      Gen gen = vec_gen.at(index_had_w[i]);
      
      float gen_eta = gen.Eta();
      float gen_phi = gen.Phi();
      float gen_pt = gen.Pt();

      for(unsigned int j=0; j<vec_jet.size(); j++)
	{
	  Jet jet = vec_jet.at(j);
	  
	  float jet_eta = jet.Eta();
	  float jet_phi = jet.Phi();
	  float jet_pt = jet.Pt();

	  float del_eta = jet_eta - gen_eta;
	  float del_phi = jet_phi - gen_phi;

	  float del_r = Sqrt(Power(del_eta, 2.)+Power(del_phi, 2.));
	  vec_dr[i].push_back(del_r);

	  float del_pt = Abs(gen_pt - jet_pt);
	  float max_pt_diff = vec_jer.at(j)*jet_pt;
	  
	  //maximum allowed pt diff for Gen matching
	  if(del_pt < 1.6*max_pt_diff) vec_chk_pt_diff[i].push_back(true);
	  else vec_chk_pt_diff[i].push_back(false);
 	}//for(int j=0; j<vec_jet.size(); j++)
    }//for(int i=0; i<2; i++)
  
  //matching with GenHFHadronMatcher information
  for(unsigned int i=0; i<vec_jet.size(); i++)
    {
      int hf_flavour = vec_hf_flavour.at(i);
      int hf_origin = vec_hf_origin.at(i); 
      
      //w is tagged by GenHFHadronMatcher
      if(Abs(hf_origin)==24)
	{
	  for(int j=0; j<2; j++)
	    {
	      if(Abs(hf_flavour)==Abs(vec_gen.at(index_had_w[j]).PID()))
		{
		  index_matched_jet[j] = i;
		  surely_matched[j] = true;
		} 
	    }
	}
    }
  
  //for light jet, try DR matching + maximum pt diff
  for(int i=0; i<2; i++)
    {
      //skip if W parton is matched surely with GenHFMatcher
      if(surely_matched[i]==true)
	{
	  dr_return[i] = vec_dr[i].at(index_matched_jet[i]);
	  continue;
	}
      
      float closest_dr = 999; 
      
      for(unsigned int j=0; j<vec_jet.size(); j++)
	{
	  int hf_origin = vec_hf_origin.at(j);
	  if(hf_origin!=-999) continue;

	  float dr = vec_dr[i].at(j);
	  bool chk_pt_diff = vec_chk_pt_diff[i].at(j);
	  if(dr<closest_dr && chk_pt_diff==true)
	    {
	      closest_dr = dr;
	      index_matched_jet[i] = j;
	    }
	}
  
      if(index_matched_jet[i]!=-999) dr_return[i] = vec_dr[i].at(index_matched_jet[i]);
      else dr_return[i] = -999;
    }
  
  //if same two jets are matched
  if(index_matched_jet[0]==index_matched_jet[1] && index_matched_jet[1]!=-999)
    {
      int index;
      if(surely_matched[0]==true && surely_matched[1]==false) index = 1;
      else if(surely_matched[0]==false && surely_matched[1]==true) index = 0;
      else
	{
	  if(dr_return[0]<dr_return[1]) index = 1;
	  else index = 0;
	} 

      vector<float> vec_dr_copy;
      for(unsigned int i=0; i<vec_dr[index].size(); i++)
	{
	  int hf_origin = vec_hf_origin.at(i);
	  bool chk_pt_diff = vec_chk_pt_diff[index].at(i);
	  if(hf_origin==-999 && chk_pt_diff==true) vec_dr_copy.push_back(vec_dr[index].at(i));
	}
      
      nth_element(vec_dr_copy.begin(), vec_dr_copy.begin()+1, vec_dr_copy.end());
      
      float second_smallest_dr = vec_dr_copy[1];
      
      //means there is no second smallest dr
      if(second_smallest_dr==0) 
	{
	 index_matched_jet[index] = -999;
	 dr_return[index] = -999;
	}
      else
	{
	  for(unsigned int i=0; i<vec_dr[index].size(); i++)
	    {
	      if(vec_dr[index].at(i)==second_smallest_dr)
		{
		  index_matched_jet[index] = i;
		  dr_return[index] = vec_dr[index].at(i);
		  break;
		}
	    }
	}
    }
  
  cout << endl;
  cout << endl;
  cout << "Daughters from W"; 
  if(selected_w>0) cout << "+" << endl;
  else cout << "-" << endl;
  cout << "(" << index_had_w[0] << ", " << vec_gen.at(index_had_w[0]).PID() <<  "), (" << index_had_w[1] << ", " <<  vec_gen.at(index_had_w[1]).PID() << ")" << endl;
  
  for(unsigned int i=0; i<vec_jet.size(); i++)
    {
      int jet_flavour = vec_hf_flavour.at(i);
      int jet_origin = vec_hf_origin.at(i);
      
      cout << "(" << jet_flavour << "," <<  jet_origin << "), ";
    }
  cout << endl;
  for(int i=0; i<2; i++)
    {
      for(unsigned int j=0; j<vec_dr[i].size(); j++)
  	{
  	  cout << i << " " << j << " " << vec_dr[i].at(j) << " " << vec_chk_pt_diff[i].at(j) << endl;
  	}
      cout << "Matched jet index = " << index_matched_jet[i] << ", surely matched = " << surely_matched[i] << ", d_r = " << dr_return[i] << endl;
    }
  
  return selected_w;
}//int Vcb_Mu::Gen_Match_W(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour , const vector<int>& vec_hf_origin, int index_matched_jet[2])

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

void Vcb_Mu::Index_Converter(const vector<Jet>& vec_sel_jet, const vector<Jet>& vec_sel_jet_match, const int index_matched_jet_match[4], int index_matched_jet[4])
{
  for(int i=0; i<4; i++)
    {
      index_matched_jet[i] = -999;

      if(index_matched_jet_match[i]==-999) continue;
      
      Jet jet_match = vec_sel_jet_match[index_matched_jet_match[i]];

      float jet_match_pt = jet_match.Pt();
      float jet_match_eta = jet_match.Eta();
      float jet_match_phi = jet_match.Phi();

      for(unsigned int j=0; j<vec_sel_jet.size(); j++)
        {
          Jet jet = vec_sel_jet.at(j);

          float jet_pt = jet.Pt();
          float jet_eta = jet.Eta();
          float jet_phi = jet.Phi();

          if(abs(jet_match_pt-jet_pt)<1e-8 && abs(jet_match_eta-jet_eta)<1e-8 && abs(jet_match_phi-jet_phi)<1e-8)
            {
              index_matched_jet[i] = j;

              break;
            }
        }//for(unsigned int j=0; j<vec_sel_jet.size(); j++)
    }//for(int i=0; i<4; i++)

  return;
}//void Vcb_Mu::Index_Converter(const vector<Jet>& vec_sel_jet, const vector<Jet>& vec_sel_jet_match, const int index_matched_jet_match[4], int index_matched_jet[4])

//////////

Gen Vcb_Mu::Neutrino(const vector<Gen>& vec_gen)
{
  //find neutrino from W
  int index_last_w = -999;
  int index_last_aw = -999;
  int index_d_w[4] = {-999, -999, -999, -999};
  
  for(unsigned int i=0; i<vec_gen.size(); ++i)
    {
      Gen gen = vec_gen[i];

      int pid = gen.PID();
      int m_index = gen.MotherIndex();

      //find last index of W+ and W-
      if(pid==24) index_last_w = i;
      if(pid==-24) index_last_aw = i;

      //find decay products of W+
      if(m_index == index_last_w && pid!=24)
        {
          if(index_d_w[0]==-999) index_d_w[0] = i;
	  else index_d_w[1] = i;
        }

      //find decay products of W-
      if(m_index == index_last_aw && pid!=-24)
        {
	  if(index_d_w[2]==-999) index_d_w[2] = i;
	  else index_d_w[3] = i;
        }
    }//for(unsigned int i=0; i<vec_gen.size(); ++i)
  
  //Find neutrino within daughters of leptonic decayed W
  int index_neutrino = 0;
  for(int i=0; i<4; ++i)
    {
      if(Abs(vec_gen[index_d_w[i]].PID())==12 || Abs(vec_gen[index_d_w[i]].PID())==14)
	{
	  index_neutrino = index_d_w[i];
	  break;
	}
    }
 
  Gen gen = vec_gen[index_neutrino];
  
  return gen;
}//Gen Neutrino(const vector<Gen>& vec_gen)

//////////
