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

  if(run_hf_contamination_tree)
    {
      outfile->cd();
      hf_contamination_tree_correct->Write();
      hf_contamination_tree_wrong->Write();
    }

  if(run_reco_eval)
    {
      outfile->cd();
      reco_eval_tree_correct->Write();
      reco_eval_tree_wrong->Write();
    }

  if(run_template)
    {
      outfile->cd();
      
      for(int i=0; i<5; i++)
	{
	  template_tree[i]->Write();
	}
    }
  
  if(run_template_truth)
    {
      outfile->cd();
      
      for(int i=0; i<5; i++)
	{
	  template_truth_tree[i]->Write();
	}
    }

}//Vcb_Mu::~Vcb_Mu()

//////////

void Vcb_Mu::initializeAnalyzer()
{
  run_debug = HasFlag("RunDebug");
  cout << "[Vcb_Mu::initializeAnalyzer] RunDebug = " << run_debug << endl;

  run_permutation_tree = HasFlag("RunPermutationTree");
  cout << "[Vcb_Mu::initializeAnalyzer] RunPermutationTree = " << run_permutation_tree << endl;

  run_hf_contamination_tree = HasFlag("RunHFContaminationTree");
  cout << "[Vcb__mu::initializeAnalyzer] RunHFContaminationTree = " << run_hf_contamination_tree << endl;

  run_reco_eval = HasFlag("RunRecoEval");
  cout << "[Vcb_Mu::initializeAnalyzer] RunRecoEval = " << run_reco_eval << endl;

  run_chi = HasFlag("RunChi");
  cout << "[Vcb_Mu::initializeAnalyzer] RunChi = " << run_chi << endl;

  run_template = HasFlag("RunTemplate");
  cout << "[Vcb_Mu::initializeAnalyzer] RunTemplate = " << run_template << endl;
  
  run_template_truth = HasFlag("RunTemplateTruth");
  cout << "[Vcb_Mu::initializeAnalyzer] RuTemplateTruth = " << run_template_truth << endl;

  run_syst = HasFlag("RunSyst");
  cout << "[Vcb_Mu::initializeAnalyzer] RunSyst = " << run_syst << endl;

  //to check additional figure of merit using hadronic w mass constrain to suppress background  
  rm_wm_constraint = HasFlag("RM_WM");
  cout << "[Vcb_Mu::initializeAnalyzer] RM_WM = " << rm_wm_constraint << endl;

  rm_bjet_energy_reg_nn = HasFlag("RM_REG_NN");
  cout << "[Vcb_Mu::initializeAnalyzer] RM_REG_NN = " << rm_bjet_energy_reg_nn << endl;
  
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

  fitter_driver = new TKinFitterDriver(DataYear, run_permutation_tree, run_chi, rm_wm_constraint, rm_bjet_energy_reg_nn);
    
  if(run_permutation_tree)
    {
      chk_matched_jets_only = false;
      
      permutation_tree_correct = new TTree("Permutation_Correct", "Permutation_Correct");
      permutation_tree_correct->Branch("weight", &weight);
      permutation_tree_correct->Branch("n_jets", &n_sel_jet);
      permutation_tree_correct->Branch("n_bjets", &n_b_jet);
      permutation_tree_correct->Branch("n_cjets", &n_c_jet);
      permutation_tree_correct->Branch("n_matched_jets", &n_matched_jets);
      permutation_tree_correct->Branch("gen_neutrino_px", &gen_neutrino_px);
      permutation_tree_correct->Branch("gen_neutrino_py", &gen_neutrino_py);
      permutation_tree_correct->Branch("gen_neutrino_pz", &gen_neutrino_pz);
      permutation_tree_correct->Branch("met_px", &met_px);
      permutation_tree_correct->Branch("met_py", &met_py);
      permutation_tree_correct->Branch("met_rebalance_px", &met_rebalance_px);
      permutation_tree_correct->Branch("met_rebalance_py", &met_rebalance_py);
      permutation_tree_correct->Branch("neutrino_pz_sol", &neutrino_pz_sol);
      permutation_tree_correct->Branch("neutrino_pz_sol_unrebal", &neutrino_pz_sol_unrebal);
      permutation_tree_correct->Branch("chk_real_neu_pz", &chk_real_neu_pz);
      permutation_tree_correct->Branch("mt_gen", &mt_gen);
      permutation_tree_correct->Branch("mt_met", &mt_met);
      permutation_tree_correct->Branch("mt_met_rebalance", &mt_met_rebalance);
      permutation_tree_correct->Branch("had_t_b_pt", &had_t_b_pt);
      permutation_tree_correct->Branch("w_u_pt", &w_u_pt);
      permutation_tree_correct->Branch("w_d_pt", &w_d_pt);
      permutation_tree_correct->Branch("lep_t_b_pt", &lep_t_b_pt);
      permutation_tree_correct->Branch("bvsc_had_t_b", &bvsc_had_t_b);
      permutation_tree_correct->Branch("cvsb_had_t_b", &cvsb_had_t_b);
      permutation_tree_correct->Branch("cvsl_had_t_b", &cvsl_had_t_b);
      permutation_tree_correct->Branch("bvsc_w_u", &bvsc_w_u);
      permutation_tree_correct->Branch("cvsb_w_u", &cvsb_w_u);
      permutation_tree_correct->Branch("cvsl_w_u", &cvsl_w_u);
      permutation_tree_correct->Branch("bvsc_w_d", &bvsc_w_d);
      permutation_tree_correct->Branch("cvsb_w_d", &cvsb_w_d);
      permutation_tree_correct->Branch("cvsl_w_d", &cvsl_w_d);
      permutation_tree_correct->Branch("bvsc_lep_t_b", &bvsc_lep_t_b);
      permutation_tree_correct->Branch("cvsb_lep_t_b", &cvsb_lep_t_b);
      permutation_tree_correct->Branch("cvsl_lep_t_b", &cvsl_lep_t_b);
      permutation_tree_correct->Branch("met_rebalance_px", &met_rebalance_px);
      permutation_tree_correct->Branch("met_rebalance_py", &met_rebalance_py);
      permutation_tree_correct->Branch("del_phi_w_u_w_d", &del_phi_w_u_w_d);
      permutation_tree_correct->Branch("del_phi_had_w_had_t_b", &del_phi_had_w_had_t_b);
      permutation_tree_correct->Branch("del_phi_lep_neu", &del_phi_lep_neu);
      permutation_tree_correct->Branch("del_phi_lep_w_lep_t_b", &del_phi_lep_w_lep_t_b);
      permutation_tree_correct->Branch("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
      permutation_tree_correct->Branch("del_eta_w_u_w_d", &del_eta_w_u_w_d);
      permutation_tree_correct->Branch("del_eta_had_w_had_t_b", &del_eta_had_w_had_t_b);
      permutation_tree_correct->Branch("del_eta_lep_neu", &del_eta_lep_neu);
      permutation_tree_correct->Branch("del_eta_lep_w_lep_t_b", &del_eta_lep_w_lep_t_b);
      permutation_tree_correct->Branch("del_eta_had_t_lep_t", &del_eta_had_t_lep_t);
      permutation_tree_correct->Branch("del_r_w_u_w_d", &del_r_w_u_w_d);
      permutation_tree_correct->Branch("del_r_had_w_had_t_b", &del_r_had_w_had_t_b);
      permutation_tree_correct->Branch("del_r_lep_neu", &del_r_lep_neu);
      permutation_tree_correct->Branch("del_r_lep_w_lep_t_b", &del_r_lep_w_lep_t_b);
      permutation_tree_correct->Branch("del_r_had_t_lep_t", &del_r_had_t_lep_t);
      permutation_tree_correct->Branch("theta_w_u_w_d", &theta_w_u_w_d);
      permutation_tree_correct->Branch("theta_had_w_had_t_b", &theta_had_w_had_t_b);
      permutation_tree_correct->Branch("theta_lep_neu", &theta_lep_neu);
      permutation_tree_correct->Branch("theta_lep_w_lep_t_b", &theta_lep_w_lep_t_b);
      permutation_tree_correct->Branch("theta_had_t_lep_t", &theta_had_t_lep_t);
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
      permutation_tree_wrong->Branch("n_jets", &n_sel_jet);
      permutation_tree_wrong->Branch("n_bjets", &n_b_jet);
      permutation_tree_wrong->Branch("n_cjets", &n_c_jet);
      permutation_tree_wrong->Branch("n_matched_jets", &n_matched_jets);
      permutation_tree_wrong->Branch("gen_neutrino_px", &gen_neutrino_px);
      permutation_tree_wrong->Branch("gen_neutrino_py", &gen_neutrino_py);
      permutation_tree_wrong->Branch("gen_neutrino_pz", &gen_neutrino_pz);
      permutation_tree_wrong->Branch("met_px", &met_px);
      permutation_tree_wrong->Branch("met_py", &met_py);
      permutation_tree_wrong->Branch("met_rebalance_px", &met_rebalance_px);
      permutation_tree_wrong->Branch("met_rebalance_py", &met_rebalance_py);
      permutation_tree_wrong->Branch("neutrino_pz_sol", &neutrino_pz_sol);
      permutation_tree_wrong->Branch("neutrino_pz_sol_unrebal", &neutrino_pz_sol_unrebal);
      permutation_tree_wrong->Branch("chk_real_neu_pz", &chk_real_neu_pz);
      permutation_tree_wrong->Branch("mt_gen", &mt_gen);
      permutation_tree_wrong->Branch("mt_met", &mt_met);
      permutation_tree_wrong->Branch("mt_met_rebalance", &mt_met_rebalance);
      permutation_tree_wrong->Branch("had_t_b_pt", &had_t_b_pt);
      permutation_tree_wrong->Branch("w_u_pt", &w_u_pt);
      permutation_tree_wrong->Branch("w_d_pt", &w_d_pt);
      permutation_tree_wrong->Branch("lep_t_b_pt", &lep_t_b_pt);
      permutation_tree_wrong->Branch("bvsc_had_t_b", &bvsc_had_t_b);
      permutation_tree_wrong->Branch("cvsb_had_t_b", &cvsb_had_t_b);
      permutation_tree_wrong->Branch("cvsl_had_t_b", &cvsl_had_t_b);
      permutation_tree_wrong->Branch("bvsc_w_u", &bvsc_w_u);
      permutation_tree_wrong->Branch("cvsb_w_u", &cvsb_w_u);
      permutation_tree_wrong->Branch("cvsl_w_u", &cvsl_w_u);
      permutation_tree_wrong->Branch("bvsc_w_d", &bvsc_w_d);
      permutation_tree_wrong->Branch("cvsb_w_d", &cvsb_w_d);
      permutation_tree_wrong->Branch("cvsl_w_d", &cvsl_w_d);
      permutation_tree_wrong->Branch("bvsc_lep_t_b", &bvsc_lep_t_b);
      permutation_tree_wrong->Branch("cvsb_lep_t_b", &cvsb_lep_t_b);
      permutation_tree_wrong->Branch("cvsl_lep_t_b", &cvsl_lep_t_b);
      permutation_tree_wrong->Branch("met_rebalance_px", &met_rebalance_px);
      permutation_tree_wrong->Branch("met_rebalance_py", &met_rebalance_py);
      permutation_tree_wrong->Branch("del_phi_w_u_w_d", &del_phi_w_u_w_d);
      permutation_tree_wrong->Branch("del_phi_had_w_had_t_b", &del_phi_had_w_had_t_b);
      permutation_tree_wrong->Branch("del_phi_lep_neu", &del_phi_lep_neu);
      permutation_tree_wrong->Branch("del_phi_lep_w_lep_t_b", &del_phi_lep_w_lep_t_b);
      permutation_tree_wrong->Branch("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
      permutation_tree_wrong->Branch("del_eta_w_u_w_d", &del_eta_w_u_w_d);
      permutation_tree_wrong->Branch("del_eta_had_w_had_t_b", &del_eta_had_w_had_t_b);
      permutation_tree_wrong->Branch("del_eta_lep_neu", &del_eta_lep_neu);
      permutation_tree_wrong->Branch("del_eta_lep_w_lep_t_b", &del_eta_lep_w_lep_t_b);
      permutation_tree_wrong->Branch("del_eta_had_t_lep_t", &del_eta_had_t_lep_t);
      permutation_tree_wrong->Branch("del_r_w_u_w_d", &del_r_w_u_w_d);
      permutation_tree_wrong->Branch("del_r_had_w_had_t_b", &del_r_had_w_had_t_b);
      permutation_tree_wrong->Branch("del_r_lep_neu", &del_r_lep_neu);
      permutation_tree_wrong->Branch("del_r_lep_w_lep_t_b", &del_r_lep_w_lep_t_b);
      permutation_tree_wrong->Branch("del_r_had_t_lep_t", &del_r_had_t_lep_t);
      permutation_tree_wrong->Branch("theta_w_u_w_d", &theta_w_u_w_d);
      permutation_tree_wrong->Branch("theta_had_w_had_t_b", &theta_had_w_had_t_b);
      permutation_tree_wrong->Branch("theta_lep_neu", &theta_lep_neu);
      permutation_tree_wrong->Branch("theta_lep_w_lep_t_b", &theta_lep_w_lep_t_b);
      permutation_tree_wrong->Branch("theta_had_t_lep_t", &theta_had_t_lep_t);
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
    }//if(run_permutation_tree)
  
  if(run_hf_contamination_tree)
    {
      chk_matched_jets_only = false;

      hf_contamination_tree_correct = new TTree("HF_Contamination_Tree_Correct", "HF_Contamination_Tree_Correct");
      hf_contamination_tree_correct->Branch("n_jets", &n_sel_jet);
      hf_contamination_tree_correct->Branch("n_bjets", &n_b_jet);
      hf_contamination_tree_correct->Branch("n_cjets", &n_c_jet);
      hf_contamination_tree_correct->Branch("best_mva_score", &best_mva_score);
      hf_contamination_tree_correct->Branch("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
      hf_contamination_tree_correct->Branch("theta_b_b", &theta_b_b);
      hf_contamination_tree_correct->Branch("theta_c_c", &theta_c_c);
      hf_contamination_tree_correct->Branch("theta_p_had_w", &theta_p_had_w);
      hf_contamination_tree_correct->Branch("theta_w_u_b", &theta_w_u_b);
      hf_contamination_tree_correct->Branch("theta_w_d_b", &theta_w_d_b);
      hf_contamination_tree_correct->Branch("w_u_b_bscore", &w_u_b_bscore);
      hf_contamination_tree_correct->Branch("w_d_b_bscore", &w_d_b_bscore);
      hf_contamination_tree_correct->Branch("m_w_u_b", &m_w_u_b);
      hf_contamination_tree_correct->Branch("m_w_d_b", &m_w_d_b);

      hf_contamination_tree_wrong = new TTree("HF_Contamination_Tree_Wrong", "HF_Contamination_Tree_Wrong");
      hf_contamination_tree_wrong->Branch("n_jets", &n_sel_jet);
      hf_contamination_tree_wrong->Branch("n_bjets", &n_b_jet);
      hf_contamination_tree_wrong->Branch("n_cjets", &n_c_jet);
      hf_contamination_tree_wrong->Branch("best_mva_score", &best_mva_score);
      hf_contamination_tree_wrong->Branch("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
      hf_contamination_tree_wrong->Branch("theta_p_had_w", &theta_p_had_w);
      hf_contamination_tree_wrong->Branch("theta_b_b", &theta_b_b);
      hf_contamination_tree_wrong->Branch("theta_c_c", &theta_c_c);
      hf_contamination_tree_wrong->Branch("theta_w_u_b", &theta_w_u_b);
      hf_contamination_tree_wrong->Branch("theta_w_d_b", &theta_w_d_b);
      hf_contamination_tree_wrong->Branch("theta_b_b", &theta_b_b);
      hf_contamination_tree_wrong->Branch("w_u_b_bscore", &w_u_b_bscore);
      hf_contamination_tree_wrong->Branch("w_d_b_bscore", &w_d_b_bscore);
      hf_contamination_tree_wrong->Branch("m_w_u_b", &m_w_u_b);
      hf_contamination_tree_wrong->Branch("m_w_d_b", &m_w_d_b);
    }

  if(run_reco_eval)
    {
      //chk_matched_jets_only = true;
      chk_matched_jets_only = false;
      
      reco_eval_tree_correct = new TTree("Reco_Eval_Tree_Correct", "Reco_Eval_Tree_Correct");
      reco_eval_tree_correct->Branch("weight", &weight);
      reco_eval_tree_correct->Branch("n_jets", &n_sel_jet);
      reco_eval_tree_correct->Branch("n_bjets", &n_b_jet);
      reco_eval_tree_correct->Branch("n_cjets", &n_c_jet);
      reco_eval_tree_correct->Branch("decay_mode", &decay_mode);
      reco_eval_tree_correct->Branch("best_mva_score", &best_mva_score);
      reco_eval_tree_correct->Branch("best_chi2", &best_chi2); 
      reco_eval_tree_correct->Branch("mt", &mt);
      reco_eval_tree_correct->Branch("chk_included", &chk_included);
      reco_eval_tree_correct->Branch("chk_hf_contamination", &chk_hf_contamination);
      reco_eval_tree_correct->Branch("mva_hf_score", &mva_hf_score);
      reco_eval_tree_correct->Branch("bvsc_had_t_b", &bvsc_had_t_b);
      reco_eval_tree_correct->Branch("cvsb_had_t_b", &cvsb_had_t_b);
      reco_eval_tree_correct->Branch("cvsl_had_t_b", &cvsl_had_t_b);
      reco_eval_tree_correct->Branch("had_t_b_pt", &had_t_b_pt);
      reco_eval_tree_correct->Branch("bvsc_w_u", &bvsc_w_u);
      reco_eval_tree_correct->Branch("cvsb_w_u", &cvsb_w_u);
      reco_eval_tree_correct->Branch("cvsl_w_u", &cvsl_w_u);
      reco_eval_tree_correct->Branch("w_u_pt", &w_u_pt);
      reco_eval_tree_correct->Branch("m_w_u", &m_w_u);
      reco_eval_tree_correct->Branch("bvsc_w_d", &bvsc_w_d);
      reco_eval_tree_correct->Branch("cvsb_w_d", &cvsb_w_d);
      reco_eval_tree_correct->Branch("cvsl_w_d", &cvsl_w_d);
      reco_eval_tree_correct->Branch("w_d", &w_d_pt);
      reco_eval_tree_correct->Branch("m_w_d", &m_w_d);
      reco_eval_tree_correct->Branch("bvsc_lep_t_b", &bvsc_lep_t_b);
      reco_eval_tree_correct->Branch("cvsb_lep_t_b", &cvsb_lep_t_b);
      reco_eval_tree_correct->Branch("cvsl_lep_t_b", &cvsl_lep_t_b);
      reco_eval_tree_correct->Branch("lep_t_b_pt", &lep_t_b_pt);
      reco_eval_tree_correct->Branch("swapped", &swapped);
      
      reco_eval_tree_wrong = new TTree("Reco_Eval_Tree_Wrong", "Reco_Eval_Tree_Wrong");
      reco_eval_tree_wrong->Branch("weight", &weight);
      reco_eval_tree_wrong->Branch("n_jets", &n_sel_jet);
      reco_eval_tree_wrong->Branch("n_bjets", &n_b_jet);
      reco_eval_tree_wrong->Branch("n_cjets", &n_c_jet);
      reco_eval_tree_wrong->Branch("decay_mode", &decay_mode);
      reco_eval_tree_wrong->Branch("best_mva_score", &best_mva_score);
      reco_eval_tree_wrong->Branch("best_chi2", &best_chi2);
      reco_eval_tree_wrong->Branch("mt", &mt);
      reco_eval_tree_wrong->Branch("chk_included", &chk_included);
      reco_eval_tree_wrong->Branch("chk_hf_contamination", &chk_hf_contamination);
      reco_eval_tree_wrong->Branch("mva_hf_score", &mva_hf_score);
      reco_eval_tree_wrong->Branch("bvsc_had_t_b", &bvsc_had_t_b);
      reco_eval_tree_wrong->Branch("cvsb_had_t_b", &cvsb_had_t_b);
      reco_eval_tree_wrong->Branch("cvsl_had_t_b", &cvsl_had_t_b);
      reco_eval_tree_wrong->Branch("had_t_b_pt", &had_t_b_pt);
      reco_eval_tree_wrong->Branch("bvsc_w_u", &bvsc_w_u);
      reco_eval_tree_wrong->Branch("cvsb_w_u", &cvsb_w_u);
      reco_eval_tree_wrong->Branch("cvsl_w_u", &cvsl_w_u);
      reco_eval_tree_wrong->Branch("w_u_pt", &w_u_pt);
      reco_eval_tree_wrong->Branch("m_w_u", &m_w_u);
      reco_eval_tree_wrong->Branch("bvsc_w_d", &bvsc_w_d);
      reco_eval_tree_wrong->Branch("cvsb_w_d", &cvsb_w_d);
      reco_eval_tree_wrong->Branch("cvsl_w_d", &cvsl_w_d);
      reco_eval_tree_wrong->Branch("w_d_pt", &w_d_pt);
      reco_eval_tree_wrong->Branch("m_w_d", &m_w_d);
      reco_eval_tree_wrong->Branch("bvsc_lep_t_b", &bvsc_lep_t_b);
      reco_eval_tree_wrong->Branch("cvsb_lep_t_b", &cvsb_lep_t_b);
      reco_eval_tree_wrong->Branch("cvsl_lep_t_b", &cvsl_lep_t_b);
      reco_eval_tree_wrong->Branch("lep_t_b_pt", &lep_t_b_pt);
      reco_eval_tree_wrong->Branch("swapped", &swapped);
    }//run_reco_eval

  if(run_template_truth)
    {
      template_truth_tree[0] = new TTree("Template_Truth_21", "Template_Truth_21");
      template_truth_tree[1] = new TTree("Template_Truth_23", "Template_Truth_23");
      template_truth_tree[2] = new TTree("Template_Truth_41", "Template_Truth_41");
      template_truth_tree[3] = new TTree("Template_Truth_43", "Template_Truth_43");
      template_truth_tree[4] = new TTree("Template_Truth_45", "Template_Truth_45");
      
      for(int i=0; i<5; i++)
	{
	  template_truth_tree[i]->Branch("best_mva_score", &best_mva_score);//dummy
	  template_truth_tree[i]->Branch("bvsc_w_u", &bvsc_w_u);
	  template_truth_tree[i]->Branch("cvsb_w_u", &cvsb_w_u);
	  template_truth_tree[i]->Branch("cvsl_w_u", &cvsl_w_u);
	  template_truth_tree[i]->Branch("m_w_u", &m_w_u);
	  template_truth_tree[i]->Branch("bvsc_w_d", &bvsc_w_d);
          template_truth_tree[i]->Branch("cvsb_w_d", &cvsb_w_d);
          template_truth_tree[i]->Branch("cvsl_w_d", &cvsl_w_d);
	  template_truth_tree[i]->Branch("m_w_d", &m_w_d);
	}
    }//run_template_truth

  result_tree = new TTree("Result_Tree", "Result_Tree");
  result_tree->Branch("weight", &weight);
  result_tree->Branch("decay_mode", &decay_mode);
  result_tree->Branch("n_jets", &n_sel_jet);
  result_tree->Branch("n_bjets", &n_b_jet);
  result_tree->Branch("n_cjets", &n_c_jet);
  result_tree->Branch("best_chi2", &best_chi2);
  result_tree->Branch("best_mva_score", &best_mva_score);
  result_tree->Branch("mt", &mt);
  result_tree->Branch("bvsc_w_u", &bvsc_w_u);
  result_tree->Branch("cvsb_w_u", &cvsb_w_u);
  result_tree->Branch("cvsl_w_u", &cvsl_w_u);
  result_tree->Branch("bvsc_w_d", &bvsc_w_d);
  result_tree->Branch("cvsb_w_d", &cvsb_w_d);
  result_tree->Branch("cvsl_w_d", &cvsl_w_d);
  result_tree->Branch("jet_mass_w_u", &jet_mass_w_u);
  result_tree->Branch("jet_mass_w_d", &jet_mass_w_d);
  result_tree->Branch("mva_hf_score", &mva_hf_score);

  reader_hf_contamination_lessthantwo = new TMVA::Reader("!Color:!Silent");
  reader_hf_contamination_lessthantwo->AddSpectator("n_jets", &n_sel_jet);
  reader_hf_contamination_lessthantwo->AddVariable("n_bjets", &n_b_jet_f);//strange...
  reader_hf_contamination_lessthantwo->AddVariable("n_cjets", &n_c_jet_f);
  reader_hf_contamination_lessthantwo->AddVariable("best_mva_score", &best_mva_score);
  reader_hf_contamination_lessthantwo->AddVariable("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
  reader_hf_contamination_lessthantwo->AddVariable("theta_b_b", &theta_b_b);
  reader_hf_contamination_lessthantwo->AddVariable("theta_p_had_w", &theta_p_had_w);
  reader_hf_contamination_lessthantwo->AddVariable("theta_w_u_b", &theta_w_u_b);
  reader_hf_contamination_lessthantwo->AddVariable("theta_w_d_b", &theta_w_d_b);
  reader_hf_contamination_lessthantwo->AddVariable("w_u_b_bscore", &w_u_b_bscore);
  reader_hf_contamination_lessthantwo->AddVariable("w_d_b_bscore", &w_d_b_bscore);
  reader_hf_contamination_lessthantwo->AddVariable("m_w_u_b", &m_w_u_b);
  reader_hf_contamination_lessthantwo->AddVariable("m_w_d_b", &m_w_d_b);
  
  TString weight_file_base = getenv("SKFlat_WD");
  weight_file_base += "/data/Run2UltraLegacy_v2/";
  weight_file_base += DataEra;
  weight_file_base += "/HF_Contamination/";
  
  TString weight_file = weight_file_base + "LessThanTwo/weights/TMVAClassification_BDTG.weights.xml";

  cout << weight_file << endl;
  reader_hf_contamination_lessthantwo->BookMVA("LessThanTwo", weight_file);

  reader_hf_contamination_morethantwo = new TMVA::Reader("!Color:!Silent");
  reader_hf_contamination_morethantwo->AddSpectator("n_jets", &n_sel_jet);
  reader_hf_contamination_morethantwo->AddVariable("n_bjets", &n_b_jet_f);//strange...
  reader_hf_contamination_morethantwo->AddVariable("n_cjets", &n_c_jet_f);
  reader_hf_contamination_morethantwo->AddVariable("best_mva_score", &best_mva_score);
  reader_hf_contamination_morethantwo->AddVariable("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
  reader_hf_contamination_morethantwo->AddVariable("theta_b_b", &theta_b_b);
  reader_hf_contamination_morethantwo->AddVariable("theta_c_c", &theta_c_c);
  reader_hf_contamination_morethantwo->AddVariable("theta_p_had_w", &theta_p_had_w);
  reader_hf_contamination_morethantwo->AddVariable("theta_w_u_b", &theta_w_u_b);
  reader_hf_contamination_morethantwo->AddVariable("theta_w_d_b", &theta_w_d_b);
  reader_hf_contamination_morethantwo->AddVariable("w_u_b_bscore", &w_u_b_bscore);
  reader_hf_contamination_morethantwo->AddVariable("w_d_b_bscore", &w_d_b_bscore);
  reader_hf_contamination_morethantwo->AddVariable("m_w_u_b", &m_w_u_b);
  reader_hf_contamination_morethantwo->AddVariable("m_w_d_b", &m_w_d_b);
  
  weight_file = weight_file_base + "MoreThanTwo/weights/TMVAClassification_BDTG.weights.xml";

  cout << weight_file << endl;
  reader_hf_contamination_morethantwo->BookMVA("MoreThanTwo", weight_file);

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
  weight = 1;

  Event ev = GetEvent();

  decay_mode = 999;
  
  if(!IsData)
    {
      //separate W decay mode
      vector<Gen> vec_gen = GetGens();
      decay_mode = Get_W_Decay_Mode(vec_gen);
      
      //lumi
      //float lumi_weight = weight_norm_1invpb*ev.GetTriggerLumi("Full");
      //weight *= lumi_weight;

      //MCweight +1 or -1
      float mc_weight = MCweight()*ev.GetTriggerLumi("Full");
      weight *= mc_weight;

      //pileup reweight
      float pileup_weight = mcCorr->GetPileUpWeight(nPileUp, 0);
      weight *= pileup_weight;

      //L1 prefire
      float prefire_weight = GetPrefireWeight(0);
      weight *= prefire_weight;
    
      // if(!TMath::Finite(lumi_weight)) cout << lumi_weight << endl;
      // if(!TMath::Finite(mc_weight)) cout << mc_weight << endl;
      // if(!TMath::Finite(pileup_weight)) cout << pileup_weight << endl;
      // if(!TMath::Finite(prefire_weight)) cout << prefire_weight << endl;
    }

  //no cut
  FillHist(param.Name+"/Cut_Flow", Cut_Flow::No_Cut, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name+Form("/Cut_Flow_%d", decay_mode), Cut_Flow::No_Cut, weight, n_cut_flow, 0, n_cut_flow);

  //met filter
  if(!PassMETFilter()) return;
  FillHist(param.Name+"/Cut_Flow", Cut_Flow::Met_Filter, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name+Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Met_Filter, weight, n_cut_flow, 0, n_cut_flow);

  //set objects
  vector<Muon> vec_this_muon = vec_muon;
  vector<Electron> vec_this_electron = vec_electron;
  vector<Jet> vec_this_jet = vec_jet;
  
  Particle met = ev.GetMETVector();

  //syst basics for objects
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
  
  //setup objects
  vector<Muon> vec_sel_muon = SelectMuons(vec_this_muon, param.Muon_Tight_ID, trig_safe_pt_cut, MUON_ETA);

  //for lepton veto
  vector<Muon> vec_muon_veto = SelectMuons(vec_this_muon, param.Muon_Tight_ID, MUON_PT_VETO, MUON_ETA);
  vector<Electron> vec_electron_veto = SelectElectrons(vec_this_electron, param.Electron_Loose_ID, ELECTRON_PT_VETO, ELECTRON_ETA);

  float jet_eta_cut = 999;
  if(DataYear==2016) jet_eta_cut = JET_ETA_2016;
  else if(DataYear==2017||DataYear==2017) jet_eta_cut = JET_ETA;
  
  //Jet selection
  vector<Jet> vec_sel_jet = SelectJets(vec_this_jet, param.Jet_ID, JET_PT, jet_eta_cut);
  vec_sel_jet = SelectJets(vec_sel_jet, "LoosePileupJetVeto", JET_PT, jet_eta_cut);
  vec_sel_jet = JetsVetoLeptonInside(vec_sel_jet, vec_electron_veto, vec_muon_veto, DR_LEPTON_VETO);
  n_sel_jet = vec_sel_jet.size();

  vector<Jet> vec_sel_jet_match = SelectJets(vec_this_jet, param.Jet_ID, JET_PT_MATCH, JET_ETA_MATCH);
  vec_sel_jet_match = SelectJets(vec_sel_jet_match, "LoosePileupJetVeto", JET_PT_MATCH, JET_ETA_MATCH);

  //sort jet as pt ordering
  sort(vec_sel_jet.begin(), vec_sel_jet.end(), PtComparing);
  sort(vec_sel_jet_match.begin(), vec_sel_jet_match.end(), PtComparing);
  
  //single lepton trigger
  if(!ev.PassTrigger(vec_mu_trig)) return;

  if(!IsData)
    {
      //SF for muon trigger effi
      float sf_mu_trig_effi;
      
      if(vec_sel_muon.size()==0) sf_mu_trig_effi = 1;
      else sf_mu_trig_effi = mcCorr->MuonTrigger_SF("POGTight", mu_trig, vec_sel_muon, 0);
      
      weight *= sf_mu_trig_effi;

      // cout << "n_muon = " << vec_sel_muon.size() << endl;
      // for(int i=0; i<vec_sel_muon.size(); i++)
      //   {
      //     Muon muon = vec_sel_muon.at(i);
      //     cout << muon.Eta() << "\t" << muon.MiniAODPt() << endl;
      //   }
      // cout << sf_mu_trig_effi << endl;
    }

  FillHist(param.Name+"/Cut_Flow", Cut_Flow::Trigger, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name+Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Trigger, weight, n_cut_flow, 0, n_cut_flow);
  
  //cut on lepton
  if(vec_sel_muon.size()!=1) return;
  if(vec_muon_veto.size()!=1) return;
  if(vec_electron_veto.size()!=0) return;
  
  Muon muon = vec_sel_muon.at(0);

  if(!IsData)
    {
      //SF for muon id
      float sf_mu_id_effi = mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muon.Eta(), muon.MiniAODPt());
      weight *= sf_mu_id_effi;

      //if(run_debug) cout << "M ID SF = " << muon.Eta() << "\t" << muon.MiniAODPt() << "\t" << sf_mu_id_effi << endl;

      //SF for muon iso
      float sf_mu_iso_effi = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muon.Eta(), muon.MiniAODPt());
      weight *= sf_mu_iso_effi;
      
      //if(run_debug) cout << "M ISO SF = " << muon.Eta() << "\t" << muon.MiniAODPt() << "\t" << sf_mu_iso_effi << endl;
    }

  FillHist(param.Name+"/Cut_Flow", Cut_Flow::Single_Lepton, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name+Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Single_Lepton, weight, n_cut_flow, 0, n_cut_flow);
  
  //cut on jet
  if(n_sel_jet<4) return;

  if(!IsData)
    {
      //SF for PUJet Veto
      float sf_pujet_veto = mcCorr->PileupJetVetoReweight(vec_sel_jet, "Loose", "");
      weight *= sf_pujet_veto;
    }


  FillHist(param.Name+"/Cut_Flow", Cut_Flow::At_Least_Four_Jets, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name+Form("/Cut_Flow_%d", decay_mode), Cut_Flow::At_Least_Four_Jets, weight, n_cut_flow, 0, n_cut_flow);
  
  //n of btag
  n_b_jet = 0;
  vector<bool> vec_btag;
  for(auto& jet : vec_sel_jet)
    {
      float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score)
        {
          n_b_jet++;
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
  
  //n of ctag
  n_c_jet = 0;
  vector<bool> vec_ctag;
  for(auto& jet : vec_sel_jet)
    {
      float cvsb = jet.GetTaggerResult(JetTagging::DeepJet_CvsB);
      float cvsl = jet.GetTaggerResult(JetTagging::DeepJet_CvsL);
      
      float cvsl_wp = -1;
      float cvsb_wp = -1;

      if(DataYear==2017)
        {
          cvsl_wp = cvsl_2017_m;
          cvsb_wp = cvsb_2017_m;
        }
      else if(DataYear==2018)
        {
          cvsl_wp = cvsl_2018_m;
          cvsb_wp = cvsb_2018_m;
        }

      if(cvsl_wp<cvsl && cvsb_wp<cvsb)
        {
          n_c_jet++;
          vec_ctag.push_back(true);
        }
      else vec_ctag.push_back(false);
    }
  //cout << "test n_of_ctag" << n_c_jet << endl;

  if(!IsData)
    {
      //SF for b-tagging
      float sf_b_tag = 1;
      
      if(run_debug) sf_b_tag *= mcCorr->GetBTaggingReweight_1a(vec_sel_jet, vec_jet_tagging_para.at(0));
      else sf_b_tag *= mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0));
      
      weight *= sf_b_tag;
    }

  if(n_b_jet<2) return;
  
  FillHist(param.Name+"/Cut_Flow", Cut_Flow::At_Least_Two_B_Tagged, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name+Form("/Cut_Flow_%d", decay_mode), Cut_Flow::At_Least_Two_B_Tagged, weight, n_cut_flow, 0, n_cut_flow);
  
  //cut on MET
  if(met.Pt()<MET_PT) return;

  FillHist(param.Name+"/Cut_Flow", Cut_Flow::MET, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name+Form("/Cut_Flow_%d", decay_mode), Cut_Flow::MET, weight, n_cut_flow, 0, n_cut_flow);

  //Data and MC comparison
  if(n_b_jet==2)
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
      FillHist(param.Name+"/TwoB/N_BJet", n_b_jet, weight, 10, 0, 10);
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
      FillHist(param.Name+"/ThreeB/N_BJet", n_b_jet, weight, 10, 0, 10);
    }
  
  //Gen matching
  vector<int> vec_hf_flavour;
  vector<int> vec_hf_origin;
  int index_matched_jet_match[4];//index of sel_jet_match matched to gen
  int index_matched_jet[4];//index of sel_jet matched to gen
  int index_gen[4];
  bool surely_matched[4];
  float matched_jet_dr[4];
  int switch_included = -1;
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

      //PrintGen(vec_gen);
      
      Gen_Match_TT(vec_sel_jet_match, vec_gen, vec_hf_flavour, vec_hf_origin, vec_jer_match, index_gen, index_matched_jet_match, surely_matched, matched_jet_dr);
      
      bool chk_gen_acceptance = true;
      for(int i=0; i<4; ++i)
	{
	  int index = index_gen[i];
	  
	  float pt = vec_gen[index].Pt();
	  float eta = vec_gen[index].Eta();
	
	  if(pt<JET_PT_MATCH || JET_ETA_MATCH<TMath::Abs(eta))
	    {
	      chk_gen_acceptance = false;
	      break;
	    }
	}
      
      if(chk_gen_acceptance) FillHist(param.Name+"/Gen_Acceptance", 1, weight, 2, 0, 2);
      else FillHist(param.Name+"/Gen_Acceptance", 0, weight, 2, 0, 2);

      cout << "test index_matched_jet_match[i]" << endl;
      for(int i=0; i<4; i++){ cout << index_matched_jet_match[i] << " " << matched_jet_dr[i] << endl; }

      //if all four jets are successfully match
      if(-1<index_matched_jet_match[0] && -1<index_matched_jet_match[1] && -1<index_matched_jet_match[2] && -1<index_matched_jet_match[3])
	{
	  FillHist(param.Name+"/PF_Gen_Matched", 1, weight, 2, 0, 2);
 
	  //BJetRegression Test
	  TLorentzVector w_gen_matched = vec_sel_jet_match.at(index_matched_jet_match[1]) + vec_sel_jet_match.at(index_matched_jet_match[2]);
	  TLorentzVector t_gen_matched = vec_sel_jet_match.at(index_matched_jet_match[0]) + w_gen_matched;

	  FillHist(param.Name+"/W_Gen_Matched_Mass", w_gen_matched.M(), weight, 100, 0, 400);
	  FillHist(param.Name+"/T_Gen_Matched_Mass", t_gen_matched.M(), weight, 100, 0, 600);

	  //correction
	  Jet j0 = vec_sel_jet_match.at(index_matched_jet_match[0]);
	  Jet j1 = vec_sel_jet_match.at(index_matched_jet_match[1]);
          Jet j2 = vec_sel_jet_match.at(index_matched_jet_match[2]);

	  float j0_corr = j0.BJetNNCorrection();
	  j0.SetPtEtaPhiM(j0_corr*j0.Pt(), j0.Eta(), j0.Phi(), j0.M());
	  
	  float j1_corr = 1;
	  if(vec_hf_flavour[index_matched_jet_match[1]]==4) j1_corr = j1.CJetNNCorrection();
	  else if(vec_hf_flavour[index_matched_jet_match[1]]==5) j1_corr = j1.BJetNNCorrection();
	  
	  j1.SetPtEtaPhiM(j1_corr*j1.Pt(), j1.Eta(), j1.Phi(), j1.M());

	  float j2_corr = 1;
	  if(vec_hf_flavour[index_matched_jet_match[2]]==4) j2_corr = j2.CJetNNCorrection();
          else if(vec_hf_flavour[index_matched_jet_match[2]]==5) j2_corr = j2.BJetNNCorrection();
	  
	  j2.SetPtEtaPhiM(j2_corr*j2.Pt(), j2.Eta(), j2.Phi(), j2.M());

	  TLorentzVector w_gen_matched_corr = j1 + j2;
          TLorentzVector t_gen_matched_corr = j0 + w_gen_matched_corr;

	  FillHist(param.Name+"/W_Gen_Matched_Mass_CorrNN", w_gen_matched_corr.M(), weight, 100, 0, 400);
          FillHist(param.Name+"/T_Gen_Matched_Mass_CorrNN", t_gen_matched_corr.M(), weight, 100, 0, 600);

	  //Delta R
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
	  switch_included = Chk_Included(index_matched_jet);
	}
      //at least of matching failed
      else 
	{
	  cout << "Test Fill 0" << endl;
	  FillHist(param.Name+"/PF_Gen_Matched", 0, weight, 2, 0, 2);
	}
      
      //-1 gen matching failed
      // 0 gen matching succeeded && four jets are in selection
      //0< gen matching succeeded &&
      //first bjet is out of selection -> 1 bit
      //first wjet is out of selection -> 2 bit
      //second wjet is out of selection -> 3 bit
      //second bjet is out of selection -> 4 bit        
      FillHist(param.Name+"/Objects_In_Sel_Jet", switch_included, weight, 18, -2, 16);
    }//if(!IsData)
  
  //To match index of jet_match (loose jet wo/ lepton veto) and jet. If 
  Index_Converter(vec_sel_jet, vec_sel_jet_match, index_matched_jet_match, index_matched_jet);

  cout << "test index_matched_jet[i]" << endl;
  for(int i=0; i<4; i++){ cout << index_matched_jet[i] << endl; }
  
  //To estimate the fraction of the included events
  if(Chk_Included(index_matched_jet)==0) FillHist(param.Name+"/Included", 1, weight, 2, 0, 2);
  else FillHist(param.Name+"/Included", 0, weight, 2, 0, 2);

  // cout << "n_size sel_jet = " << vec_sel_jet.size() << ", Index of matched_sel_jet[4]" << endl;
  // for(int i=0; i<4; i++)
  //        {
  //          cout << index_matched_jet[i] << " ";
  //        }
  // cout << endl;

  //
  if(!IsData && run_template_truth)
    {
      if(index_matched_jet[1]<0 || index_matched_jet[2]<0) return;

      Jet jet_w_u = vec_sel_jet[index_matched_jet[1]];
      Jet jet_w_d = vec_sel_jet[index_matched_jet[2]];

      best_mva_score = 1;

      bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);
      m_w_u = jet_w_u.GetM();

      bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);
      m_w_d = jet_w_d.GetM();

      int index = -1;
      if(decay_mode==21) index = 0;
      else if(decay_mode==23) index = 1;
      else if(decay_mode==41) index = 2;
      else if(decay_mode==43) index = 3;
      else if(decay_mode==45) index = 4;
      if(index!=-1) template_truth_tree[index]->Fill();

      return;
    }//if(!IsData && run_template_truth)

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
  
  //HF contamination score
  Results_Container results_container = fitter_driver->Get_Results();

  n_b_jet_f = n_b_jet;
  n_c_jet_f = n_c_jet;
  best_mva_score = results_container.best_mva_score;
  del_phi_had_t_lep_t = results_container.best_del_phi_had_t_lep_t;

  theta_b_b = 999;
  for(int i=0; i<n_sel_jet; i++)
    {
      if(!vec_btag[i])  continue;

      Jet b0 = vec_sel_jet[i];

      for(int j=i+1; j<n_sel_jet; j++)
        {
          if(!vec_btag[j]) continue;

          Jet b1 = vec_sel_jet[j];

          float theta_b_b_new = b0.Angle(b1.Vect());

          if(theta_b_b_new<theta_b_b) theta_b_b = theta_b_b_new;
        }
    }
  
  theta_c_c = 999;
  for(int i=0; i<n_sel_jet; i++)
    {
      if(!vec_ctag[i]) continue;
      
      Jet c0 = vec_sel_jet[i];
      
      for(int j=i+1; j<n_sel_jet; j++)
	{
	  if(!vec_ctag[j]) continue;

	  Jet c1 = vec_sel_jet[j];

	  float theta_c_c_new = c0.Angle(c1.Vect());

	  if(theta_c_c_new<theta_c_c) theta_c_c = theta_c_c_new;
	}
    }

  int index[2];
  index[0] = results_container.best_index_w_u;
  index[1] = results_container.best_index_w_d;

  for(int i=0; i<2; i++)
    {
      Jet jet_w_candi = vec_sel_jet[index[i]];

      int index_closest_b = -1;
      float theta_b = 999;
      //float mass_b = 999;
      for(int j=0; j<n_sel_jet; j++)
	{
	  if(!vec_btag[j]) continue;
	  
	  Jet bjet = vec_sel_jet[j];

	  //if two jets are same
	  if(Compare_Jet(jet_w_candi, bjet)) continue;

	  float theta_b_new = jet_w_candi.Angle(bjet.Vect());
	  
	  if(theta_b_new<theta_b)
	    {
	      index_closest_b = j;
	      theta_b = theta_b_new;
	    }
	}//loop over bjets
      
      if(i==0)
	{
	  theta_w_u_b = theta_b;

	  w_u_b_bscore = vec_sel_jet[index_closest_b].GetTaggerResult(JetTagging::DeepJet);
	  m_w_u_b = (jet_w_candi+vec_sel_jet[index_closest_b]).M();
	}
      else if(i==1)
	{
	  theta_w_d_b = theta_b;

	  w_d_b_bscore = vec_sel_jet[index_closest_b].GetTaggerResult(JetTagging::DeepJet);
	  m_w_d_b = (jet_w_candi+vec_sel_jet[index_closest_b]).M();
	}
    }//loop over w candidate
  
  if(n_c_jet<=1) mva_hf_score = reader_hf_contamination_lessthantwo->EvaluateMVA("LessThanTwo");
  else if(2<=n_c_jet) mva_hf_score = reader_hf_contamination_morethantwo->EvaluateMVA("MoreThanTwo");
  
  //Permutation Tree for signal study 
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
      unsigned int index_correct = 999999;//should be very large number
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
      
      //calculate neutrino pz solution with unrebalanced met and find best one
      float diff_pz_unrebal[2];
      float neu_pz_sol_unrebal[2];
      
      Sol_Neutrino_Pz(muon, met, neu_pz_sol_unrebal);
      
      diff_pz_unrebal[0] = TMath::Abs(gen_neutrino_pz - neu_pz_sol_unrebal[0]);
      diff_pz_unrebal[1] = TMath::Abs(gen_neutrino_pz - neu_pz_sol_unrebal[1]);
      
      if(diff_pz_unrebal[0]<diff_pz_unrebal[1]) neutrino_pz_sol_unrebal = neu_pz_sol_unrebal[0];
      else neutrino_pz_sol_unrebal = neu_pz_sol_unrebal[1];

      for(unsigned int i=0; i<results_container.vec_results.size(); ++i)
	{
	  Results results = results_container.vec_results[i];
	  
	  met_rebalance_px = results.met_rebalance_px;
	  met_rebalance_py = results.met_rebalance_py;
	  neutrino_pz_sol = results.neutrino_pz_sol;
	  chk_real_neu_pz = results.chk_real_neu_pz;

	  mt_gen = Calculate_Mt(muon, gen_neutrino_px, gen_neutrino_py); 
	  mt_met = Calculate_Mt(muon, met_px, met_py);
	  mt_met_rebalance = Calculate_Mt(muon, met_rebalance_px, met_rebalance_py);

          int index_had_t_b = results.index_had_t_b;
          int index_w_u = results.index_w_u;
          int index_w_d = results.index_w_d;
          int index_lep_t_b = results.index_lep_t_b;

	  n_matched_jets = 0;
	  if(index_matched_jet[0]==index_had_t_b) n_matched_jets++;
	  if(index_matched_jet[1]==index_w_u) n_matched_jets++;
          if(index_matched_jet[2]==index_w_d) n_matched_jets++;
	  if(index_matched_jet[3]==index_lep_t_b) n_matched_jets++; 
	  
	  Jet jet_had_t_b = vec_sel_jet[index_had_t_b];
	  Jet jet_w_u = vec_sel_jet[index_w_u];
	  Jet jet_w_d = vec_sel_jet[index_w_d];
	  Jet jet_lep_t_b = vec_sel_jet[index_lep_t_b];

	  //pt
	  had_t_b_pt = jet_had_t_b.Pt();
	  w_u_pt = jet_w_u.Pt();
       	  w_d_pt = jet_w_d.Pt();
       	  lep_t_b_pt = jet_lep_t_b.Pt();

	  //tagging scores
	  bvsc_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet);
	  cvsb_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
	  cvsl_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);
       	  
	  bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
	  cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
	  cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);

	  bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
	  cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
          cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);

	  bvsc_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet);
	  cvsb_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
	  cvsl_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);
	  
	  //angles
	  del_phi_w_u_w_d = results.del_phi_w_u_w_d;
	  del_phi_had_w_had_t_b = results.del_phi_had_w_had_t_b;
	  del_phi_lep_neu = results.del_phi_lep_neu;
	  del_phi_lep_w_lep_t_b = results.del_phi_lep_w_lep_t_b;
       	  del_phi_had_t_lep_t = results.del_phi_had_t_lep_t;
	  
	  del_eta_w_u_w_d = results.del_eta_w_u_w_d;
          del_eta_had_w_had_t_b = results.del_eta_had_w_had_t_b;
          del_eta_lep_neu = results.del_eta_lep_neu;
          del_eta_lep_w_lep_t_b = results.del_eta_lep_w_lep_t_b;
          del_eta_had_t_lep_t = results.del_eta_had_t_lep_t;

	  del_r_w_u_w_d = results.del_r_w_u_w_d;
          del_r_had_w_had_t_b = results.del_r_had_w_had_t_b;
          del_r_lep_neu = results.del_r_lep_neu;
          del_r_lep_w_lep_t_b = results.del_r_lep_w_lep_t_b;
          del_r_had_t_lep_t = results.del_r_had_t_lep_t;

	  theta_w_u_w_d = results.theta_w_u_w_d;
          theta_had_w_had_t_b = results.theta_had_w_had_t_b;
          theta_lep_neu = results.theta_lep_neu;
          theta_lep_w_lep_t_b = results.theta_lep_w_lep_t_b;
          theta_had_t_lep_t = results.theta_had_t_lep_t;

       	  had_t_mass = results.initial_had_t_mass;
       	  had_w_mass = results.initial_had_w_mass;
      	  lep_t_mass = results.initial_lep_t_mass;
	  
      	  TLorentzVector lep_t_partial = jet_lep_t_b + muon;
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
    }//if(!IsData && run_permutation_tree)
    
  //HF Contamination
  if(!IsData && run_hf_contamination_tree)
    {
      cout << "test run_hf_contamination_tree" << endl;
      
      //if(n_b_jet<3) return; 

      Results_Container results_container = fitter_driver->Get_Results();
      if(!fitter_driver->Check_Status()) return;
                 
      chk_hf_contamination = false;

      int index[2];
      index[0] = results_container.best_index_w_u;
      index[1] = results_container.best_index_w_d;
      
      for(int i=0; i<2; i++)
	{
	  Jet jet = vec_sel_jet[index[i]];
	  
	  //int hf_flavour = jet.GenHFHadronMatcherFlavour();
	  int hf_origin = jet.GenHFHadronMatcherOrigin();

	  if(hf_origin==21) chk_hf_contamination = true;
	}

      //vec_bjet for easy handling
      vector<Jet> vec_bjet;
      for(int i=0; i<n_sel_jet; i++)
	{
	  if(vec_btag[i]==true) vec_bjet.push_back(vec_sel_jet[i]);
	}

      //vec_cjet for easy handling
      vector<Jet> vec_cjet;
      for(int i=0; i<n_sel_jet; i++)
	{
	  if(vec_ctag[i]==true) vec_cjet.push_back(vec_sel_jet[i]);
	}

      best_mva_score = results_container.best_mva_score;
      del_phi_had_t_lep_t = results_container.best_del_phi_had_t_lep_t;
      
      theta_b_b = 999;
      for(int i=0; i<n_b_jet; i++)
        {
          Jet b0 = vec_bjet[i];

          for(int j=i+1; j<n_b_jet; j++)
            {
              Jet b1 = vec_bjet[j];

              float theta_b_b_new = b0.Angle(b1.Vect());

              if(theta_b_b_new<theta_b_b) theta_b_b = theta_b_b_new;
            }
        }

      //not useful... 
      theta_c_c = 999;
      for(int i=0; i<n_c_jet; i++)
	{
	  Jet c0 = vec_cjet[i];
	  
	  for(int j=i+1; j<n_c_jet; j++)
	    {
	      Jet c1 = vec_cjet[j];

	      float theta_c_c_new = c0.Angle(c1.Vect());

	      if(theta_c_c_new<theta_c_c) theta_c_c = theta_c_c_new;
	    }
	}

      Jet jet_w_u = vec_sel_jet[index[0]];
      Jet jet_w_d = vec_sel_jet[index[1]];
     
      theta_p_had_w = results_container.best_theta_w_u_w_d*(jet_w_u+jet_w_d).P();

      for(int i=0; i<2; i++)
	{
	  Jet jet_w_candi = vec_sel_jet[index[i]];
	  
	  int index_closest_b = -1;
	  //int index_smallest_b = -1;
	  float theta_b = 999;
	  //float mass_b = 999;
	  for(int j=0; j<n_b_jet; j++)
	    {
	      Jet bjet = vec_bjet[j];
	      
	      //if two jets are same
	      if(Compare_Jet(jet_w_candi, bjet)) continue;
	      
	      float theta_b_new = jet_w_candi.Angle(bjet.Vect());
	      //float mass_b_new = (jet_w_candi+bjet).M();

	      if(theta_b_new<theta_b)
		{
		  index_closest_b = j;
		  theta_b = theta_b_new;
		}

	      //if(mass_b_new<mass_b) mass_b = mass_b_new;
	    }//loop over bjets
	  
	  if(i==0)
	    {
	      theta_w_u_b = theta_b;
	      
	      w_u_b_bscore = vec_bjet[index_closest_b].GetTaggerResult(JetTagging::DeepJet);
	      m_w_u_b = (jet_w_candi+vec_bjet[index_closest_b]).M();
	    }
	  else if(i==1)
	    {
	      theta_w_d_b = theta_b;
	      
	      w_d_b_bscore = vec_bjet[index_closest_b].GetTaggerResult(JetTagging::DeepJet);
	      m_w_d_b = (jet_w_candi+vec_bjet[index_closest_b]).M();
	    }
	}//loop over w candidate

      if(chk_hf_contamination) hf_contamination_tree_correct->Fill();
      else hf_contamination_tree_wrong->Fill();
      
      return;
    }//if(!IsData && run_hf_contamination_tree)

  //KF Evaluation
  if(!IsData && run_reco_eval)
    {
      cout << "test run_reco_eval" << endl;

      Results_Container results_container = fitter_driver->Get_Results();
      if(!fitter_driver->Check_Status())
	{
	  cout << "Test Fitter fail" << endl;

	  if(Chk_Included(index_matched_jet)==0) chk_included = true;
	  else chk_included = false;
	  
	  if(chk_included) FillHist(param.Name+"/KF_Eval_Fail", n_sel_jet, n_b_jet, 1, 10, 0, 10, 8, 0, 8);
	
	  return;
	}
      
      //chk included
      if(Chk_Included(index_matched_jet)==0) chk_included = true;
      else chk_included = false;
      cout << "Chk Included = " << chk_included << endl;

      //Fitter correct
      int index_had_t_b = results_container.best_index_had_t_b;
      int index_w_u = results_container.best_index_w_u;
      int index_w_d = results_container.best_index_w_d;
      int index_lep_t_b = results_container.best_index_lep_t_b;
      
      // cout << "run_reco_eval" << endl;
      // cout << index_matched_jet_match[1] << " " << index_matched_jet_match[2] << endl;
      // cout << index_w_u << " " << index_w_d << endl;

      int fitter_correct = 0;
      // if(-1<index_matched_jet_match[1] && -1<index_matched_jet[2])
      // 	{
      // 	  Jet jet_gen_matched[2] = {vec_sel_jet_match.at(index_matched_jet_match[1]), vec_sel_jet_match.at(index_matched_jet_match[2])};
      // 	  Jet jet_kf_matched[2] = {vec_sel_jet.at(index_w_u), vec_sel_jet.at(index_w_d)};
	  
      // 	  fitter_correct = Compare_Jet_Pair(jet_gen_matched, jet_kf_matched);
      // 	}
      // else fitter_correct = 0;

      // if(fitter_correct==0) swapped = -1;
      // else if(fitter_correct==1) swapped = 0;
      // else if(fitter_correct==2) swapped = 1;

      if(index_matched_jet[1]==index_w_u && index_matched_jet[2]==index_w_d)
      	{
      	  fitter_correct = 1;
      	  swapped = 0;
      	}
      else if(index_matched_jet[1]==index_w_d && index_matched_jet[2] ==index_w_u)
      	{
      	  fitter_correct = 1;
      	  swapped = 1;
      	}
      else 
      	{
      	  fitter_correct = 0;
      	  swapped = -1;
      	}

      // cout << "KF Correct = " << fitter_correct << endl;
      // cout << "Swapped = " << swapped << endl;
      
      best_chi2 = results_container.best_chi2;
      best_mva_score = results_container.best_mva_score;

      if(Chk_Included(index_matched_jet)==0) chk_included = true;
      else chk_included = false;

      if(chk_included)
	{
	  if(fitter_correct)
	    {
	      FillHist(param.Name+"/KF_Eval_Correct", n_sel_jet, n_b_jet, 1, 10, 0, 10, 8, 0, 8);
	      FillHist(param.Name+"/Chi2", best_chi2, 1, 1, 200, 0, 200, 2, 0, 2);
	      FillHist(param.Name+"/MVA_Score", best_mva_score, 1, 1, 200, -1, 1, 2, 0, 2);
	    }
	  else
	    {
	      FillHist(param.Name+"/KF_Eval_Wrong", n_sel_jet, n_b_jet, 1, 10, 0, 10, 8, 0, 8);
	      FillHist(param.Name+"/Chi2", best_chi2, 0, 1, 200, 0, 200, 2, 0, 2);
	      FillHist(param.Name+"/MVA_Score", best_mva_score, 0, 1, 200, -1, 1, 2, 0, 2);
	    }
	}
      
      //mt of W 
      met_px = met.Px();
      met_py = met.Py();
      mt = Calculate_Mt(muon, met_px, met_py);
            
      //HF contamination
      cout << "test HF contamination" << endl;
      Jet jet_w[2];
      jet_w[0] = vec_sel_jet[index_w_u];
      jet_w[1] = vec_sel_jet[index_w_d];
      
      chk_hf_contamination = false;
      for(int i=0; i<2; ++i)
	{
	  int hf_flavour = jet_w[i].GenHFHadronMatcherFlavour();
	  int hf_origin = jet_w[i].GenHFHadronMatcherOrigin();

	  cout << "test hf_flavour = " << hf_flavour << ", hf_origin = " << hf_origin << endl; 

	  if(hf_origin==21) chk_hf_contamination = true;
	}
      cout << "HF Contamination = " << chk_hf_contamination << endl;
      cout << endl;

      Jet jet_had_t_b = vec_sel_jet.at(index_had_t_b);
      Jet jet_w_u = vec_sel_jet.at(index_w_u);
      Jet jet_w_d = vec_sel_jet.at(index_w_d);
      Jet jet_lep_t_b = vec_sel_jet.at(index_lep_t_b);

      bvsc_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet);
      cvsb_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);
      had_t_b_pt = jet_had_t_b.Pt();

      bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);
      w_u_pt = jet_w_u.Pt();
      m_w_u = jet_w_u.GetM();      

      bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);
      w_d_pt = jet_w_d.Pt();
      m_w_d = jet_w_d.GetM();

      bvsc_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet);
      cvsb_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);
      lep_t_b_pt = jet_lep_t_b.Pt();

      if(fitter_correct) reco_eval_tree_correct->Fill();
      else reco_eval_tree_wrong->Fill();
      
      return;
    }//if(!IsData && run_reco_eval)

  //Template_Truth
  if(!IsData && run_template_truth)
    {
      Results_Container results_container = fitter_driver->Get_Results();
      if(!fitter_driver->Check_Status()) return;

      best_mva_score = results_container.best_mva_score;

      int index_had_t_b = results_container.best_index_had_t_b;
      int index_w_u = results_container.best_index_w_u;
      int index_w_d = results_container.best_index_w_d;
      int index_lep_t_b = results_container.best_index_lep_t_b;

      if(index_matched_jet[0]!=index_had_t_b || index_matched_jet[1]!=index_w_u ||
	 index_matched_jet[2]!=index_w_d || index_matched_jet[3]!=index_lep_t_b) return;

      Jet jet_w_u = vec_sel_jet[index_w_u];
      Jet jet_w_d = vec_sel_jet[index_w_d];

      bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);

      bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);

      int index = -1;
      if(decay_mode==21) index = 0;
      else if(decay_mode==23) index = 1;
      else if(decay_mode==41) index = 2;
      else if(decay_mode==43) index = 3;
      else if(decay_mode==45) index = 4;
      if(index!=-1) template_truth_tree[index]->Fill(); 
      
      return;
    }//if(!IsData && run_template_truth)

  bool chk_fitter_status = fitter_driver->Check_Status();
  if(chk_fitter_status)
    {
      FillHist(param.Name+Form("/Cut_Flow_%d", decay_mode), Cut_Flow::KF_Pass, weight, n_cut_flow, 0, n_cut_flow);
      
      Results_Container results_container = fitter_driver->Get_Results();
      
      best_chi2 = results_container.best_chi2;
      best_mva_score = results_container.best_mva_score;
      mt = Calculate_Mt(muon, met_px, met_py);

      int index_had_t_b = results_container.best_index_had_t_b;
      int index_w_u = results_container.best_index_w_u;
      int index_w_d = results_container.best_index_w_d;
      int index_lep_t_b = results_container.best_index_lep_t_b;
      
      Jet jet_had_t_b = vec_sel_jet.at(index_had_t_b);
      Jet jet_w_u = vec_sel_jet.at(index_w_u);
      Jet jet_w_d = vec_sel_jet.at(index_w_d);
      Jet jet_lep_t_b = vec_sel_jet.at(index_lep_t_b);
            
      /* Template */
      bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);
      
      bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
      cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
      cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);

      jet_mass_w_u = jet_w_u.M();
      jet_mass_w_d = jet_w_d.M();

      result_tree->Fill();

      FillHist(param.Name+"/B_vs_C_W_U", bvsc_w_u, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_B_W_U", cvsb_w_u, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_L_W_U", cvsl_w_u, weight, 10, 0, 1);
    
      FillHist(param.Name+"/B_vs_C_W_D", bvsc_w_d, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_B_W_D", cvsb_w_d, weight, 10, 0, 1);
      FillHist(param.Name+"/C_vs_L_W_D", cvsl_w_d, weight, 10, 0, 1);

      //n b-tagged
      float n_b_jet_w = 0;
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < bvsc_w_u) n_b_jet_w++;
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < bvsc_w_d) n_b_jet_w++;
            
      FillHist(param.Name+"/N_B_Jet_W_Candidate", n_b_jet_w, weight, 3, 0, 3);

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
    }     
 
  
  return;
}//Vcb_Mu::executeEventFromParameter(AnalyzerParameter param)

//////////

float Vcb_Mu::Calculate_Mt(const Particle& lepton, const float& neu_px, const float& neu_py)
{
  //float px = lepton.Px() + neu_px;
  //float py = lepton.Py() + neu_py;
  
  //float mt = TMath::Sqrt(W_MASS*W_MASS + px*px + py*py);

  float lepton_px = lepton.Px();
  float lepton_py = lepton.Py();
  float lepton_et = Sqrt(lepton_px*lepton_px + lepton_py*lepton_py);
  
  float neu_et = Sqrt(neu_px*neu_px + neu_py*neu_py);

  TVector3 lepton3(lepton_px, lepton_py, 0);
  TVector3 neu3(neu_px, neu_py, 0);

  float angle = lepton3.Angle(neu3);
  
  float mt = Sqrt(lepton_et*neu_et*(1-Cos(angle)));

  return mt;
}//float Vcb_Mu::Calculate_Mt(const Particle& lepton, const Particle& met)

//////////

int Vcb_Mu::Chk_Included(const int index_matched_jet[4])
{
  int result = 0;
  unsigned int tmp = 1;
  for(int i=0; i<4; ++i)
    {     
      if(index_matched_jet[i]<0) result += tmp;

      tmp = tmp << 1;
    }//for(int i=0; i<4; i++)

  //0: every gen jet is matched to loose jets wo/ lepton veto and the loose jets passes baseline selection  
  return result;
}//int Vcb_Mu::Chk_Included(const int index_matched_jet[4])

//////////

bool Vcb_Mu::Compare_Jet(const Jet& jet0, const Jet& jet1)
{
  float pt0 = jet0.Pt();
  float pt1 = jet1.Pt();

  if(Abs(pt0-pt1)<1e-8) return true;
  else return false;
}//bool Vcb_Mu::Compare_Jet(const Jet& jet0, const Jet& jet1)

//////////

int Vcb_Mu::Compare_Jet_Pair(const Jet jet0[2], const Jet jet1[2])
{
  Double_t pt0[2] = {jet0[0].Pt(), jet0[1].Pt()};
  // float eta0[2] = {jet0[0].Eta(), jet0[1].Eta()};
  // float phi0[2] = {jet0[0].Phi(), jet0[1].Phi()};

  Double_t pt1[2] = {jet1[0].Pt(), jet1[1].Pt()};
  // float eta1[2] = {jet1[0].Eta(), jet1[1].Eta()};
  // float phi1[2] = {jet1[0].Phi(), jet1[1].Phi()};

  //Same jet pair
  if(Abs(pt0[0]-pt1[0])<1e-8 && Abs(pt0[1]-pt1[1])<1e-8) return 1;
  
  //Same jet pair but, swapped
  else if(Abs(pt0[0]-pt1[1])<1e-8 && Abs(pt0[1]-pt1[0])<1e-8) return 2; 
  
  //different jet pair
  return 0;
}//bool Vcb_Mu::Compare_Jet_Pair(const Jet& jet0, const Jet& jet1)

//////////

void Vcb_Mu::Gen_Match_Residual(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, const vector<float>& vec_jer, int index_gen[4], int index_matched_jet[4], bool surely_matched[4], float dr_return[4])
{
  //For gen particles which are not matched surely, this method will try to match those to jeco jet using DR matching
  
  const float rel_pt_diff_cut = 1.5;
  const float dr_cut = 0.4;

  int possible_matched_index[4] = {-999, -999, -999, -999};
  float dr_smallest[4] = {999, 999, 999, 999};
  for(int i=0; i<4; ++i)
    {
      //cout << "gen " << i << " " << index_gen[i] << " " << index_matched_jet[i] << endl;
      
      if(index_matched_jet[i]!=-999) continue;
      
      Gen gen = vec_gen[index_gen[i]];
      
      float gen_pt = gen.Pt();
      float gen_eta = gen.Eta();
      float gen_phi = gen.Phi();
     
      //cout << "gen pt = " << gen_pt << " " << gen_eta << " " << gen_phi << endl;

      //acceptance test. If gen particle is out of acceptance, matching is not tried 
      // if(JET_ETA_MATCH<Abs(gen_eta))
      //  	{
      // 	  //cout << "Out of acceptance" << endl;
      // 	  index_matched_jet[i] = -1;//gen is out of acceptance
      //  	  continue;
      //  	}

      for(unsigned int j=0; j<vec_jet.size(); ++j)
	{
	  //if jet is aleady matched surely with other gen particles, continue
	  bool chk_allocated = false;
	  for(int k=0; k<4; ++k)
	    {
	      if(j==index_matched_jet[k])
		{
		  chk_allocated = true;
		  break;
		}
	    }
	  if(chk_allocated) continue;
	  
	  //if origin is tagged by GenHFHadronMatcher, continue
	  //int jet_flavour = vec_hf_flavour.at(j);
          int jet_origin = vec_hf_origin.at(j);
	  	  
	  if(jet_origin!=-999) continue;
	  
	  Jet jet = vec_jet.at(j);

	  float jet_pt = jet.Pt();
	  float jet_eta = jet.Eta();
	  float jet_phi = jet.Phi();
	  
	  float del_pt = gen_pt - jet_pt;
	  float del_eta = gen_eta - jet_eta;
	  float del_phi = gen_phi - jet_phi;

	  float rel_pt_diff = Abs(del_pt/vec_jer[i]/jet_pt);
	  
	  float dr = Sqrt(del_eta*del_eta + del_phi*del_phi);
	  
	  if(rel_pt_diff<rel_pt_diff_cut && dr<dr_smallest[i]) 
	    {
	      dr_smallest[i] = dr;
	      possible_matched_index[i] = j;
	    }

	  //cout << j << " " << jet_flavour << " " << jet_origin << " " << rel_pt_diff << " " << dr << " " << possible_matched_index[i] << endl;
	  
	}//loop over jet
      
      //if no jet is allocated to gen particle
      if(possible_matched_index[i]==-999 || dr_cut<dr_smallest[i])
	{
	  //cout << "impossible to allocate jet to gen " << i << endl;
	  index_matched_jet[i] = -2;//impossible to allocate any jet to gen particle
	}

      //cout << endl;
    }//loop over gen 
  
  //if two or more gen particles are matched to same jet, gen particle with smallest dr is chosen
  //then recursively call this method to allcate jet to the second ranked gen 
  for(int i=0; i<3; ++i)
    {
      if(index_matched_jet[i]!=-999) continue;
      
      bool chk_overlap[4] = {false, false, false, false};
      chk_overlap[i] = true;    
      
      for(int j=i+1; j<4; ++j)
	{
	  if(index_matched_jet[j]!=-999) continue;

	  //same jet allocation found
	  if(possible_matched_index[i]==possible_matched_index[j]) chk_overlap[j] = true;
	}
      
      //find gen particle with smallest dr among overlaped gen
      int gen_index_best_among_overlap = -999; 
      float dr_smallest_among_overlap = 999;
      for(int j=i; j<4; ++j)
	{
	  if(chk_overlap[j]!=true) continue;
		  
	  if(dr_smallest[j]<dr_smallest_among_overlap)
	    {
	      gen_index_best_among_overlap = j;
	      dr_smallest_among_overlap = dr_smallest[j];
	    }
	}
      
      index_matched_jet[gen_index_best_among_overlap] = possible_matched_index[gen_index_best_among_overlap];
      dr_return[gen_index_best_among_overlap] = dr_smallest[gen_index_best_among_overlap];

      Gen_Match_Residual(vec_jet, vec_gen, vec_hf_flavour, vec_hf_origin, vec_jer, index_gen, index_matched_jet, surely_matched, dr_return);
    }
  
  return;
}//void Vcb_Mu::Gen_Match_Residual(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, int index_gen[4], int index_matched_jet[4], bool surely_matched[4])

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
  
  // if(run_debug)
  //   {
  //     cout << endl;
  //     cout << endl;
  //     if(selected_w == 1) cout << "W+ decayed hadronically" << endl;
  //     else cout << "W- decayed hadronically" << endl;
      
  //     for(unsigned int i=0; i<4; i++)
  //     	{
  //     	  if(index_gen[i]!=-999)
  //     	    {
  //     	      Gen gen = vec_gen.at(index_gen[i]);
  //     	      int pid = gen.PID();
	      
  //     	      Gen gen_mother = vec_gen.at(gen.MotherIndex());
  //     	      int m_pid = gen_mother.PID();

  //     	      float gen_eta = gen.Eta();
  //     	      float  gen_phi = gen.Phi();
  //     	      float gen_pt = gen.Pt();
	      
  //     	      cout << "(" << index_gen[i] << ", " << pid << ", " << m_pid << ", " << gen_pt << ", " << gen_eta << ", " << gen_phi << "), ";
	      
  //     	    }
  //     	  else 
  //     	    {
  //     	      cout << "(), ";
  //     	    }
  //     	}
  //     cout << endl;
           
  //     for(unsigned int i=0; i<vec_jet.size(); i++)
  //     	{
  //     	  int jet_flavour = vec_hf_flavour.at(i);
  //         int jet_origin = vec_hf_origin.at(i);
	  
  //     	  Jet jet = vec_jet.at(i);
	  
  //     	  float jet_eta = jet.Eta();
  //     	  float jet_phi = jet.Phi();
  //     	  float jet_pt =  jet.Pt();
	  
  //     	  float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
  //     	  bool b_tagged = false;
  //     	  if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score) b_tagged = true;
	  
  //     	  cout << "(" << jet_flavour << ", " << jet_origin << ", " << jet_pt << ", " << jet_eta << ", " << jet_phi << ", " << ", " << b_tagged << ")" << endl;
  //     	}
  //     cout << endl;
      
  //     for(int i=0; i<4; i++)
  // 	{
  // 	  cout << "Matched jet index = " << index_matched_jet[i] << ", surely matched = " << surely_matched[i] << ", d_r = " << dr_return[i] << endl;
  // 	}
  //     cout << endl;
  //   }
  
  Gen_Match_Residual(vec_jet, vec_gen, vec_hf_flavour, vec_hf_origin, vec_jer, index_gen, index_matched_jet, surely_matched, dr_return);

  return;
}//void Vcb_Mu::Gen_Match_TT(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, int index_matched_jet[4], bool surely_matched[4], float dr_return[4])

//////////

int Vcb_Mu::Gen_Match_W(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, const vector<float>& vec_jer, int index_gen[2], int index_matched_jet[2], bool surely_matched[2], float dr_return[2])
{
  //this method is designed to find gen particles from W and to try to match to reco jet via GenHFHadronMatcher results
 
  const int target_pdg_id_positive = 24; 
  const int target_pdg_id_negative = -24;
  
  //For CHtoCB decay
  //const int target_pdg_id_positive = 37;
  //const int target_pdg_id_negative = -37;

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
      if(pid==target_pdg_id_positive) index_last_w = i;
      if(pid==target_pdg_id_negative) index_last_aw = i;
      
      //find decay products of W+
      if(m_index==index_last_w && pid!=target_pdg_id_positive)
	{
	  if(index_d0_w==-999)  index_d0_w = i;
	  else index_d1_w = i;
	}
      
      //find decay products of W-
      if(m_index==index_last_aw && pid!=target_pdg_id_negative)
	{
	  if(index_d0_aw==-999)  index_d0_aw = i;
          else index_d1_aw = i;
	}
    }
  
  //if input sample is not TT, so both of W couldn't be found 
  if(index_last_w==-999 || index_last_aw==-999)
    {
      //if(run_debug) cout << "Can't find both of W" << endl;
      return -999;
    }
  
  //W+ decay hadronically
  if(Abs(vec_gen.at(index_d0_w).PID())<10)
    {
      if(vec_gen.at(index_d0_w).PID()%2==0)
	{
	  index_gen[0] = index_d0_w;
	  index_gen[1] = index_d1_w;
	}
      else 
	{
	  index_gen[0] = index_d1_w;
	  index_gen[1] = index_d0_w;
	}
      
      selected_w = 1;//mean W+ decayed hadronically
    }
  //W- decay hadronically
  else if(Abs(vec_gen.at(index_d0_aw).PID()<10))
    {
      if(vec_gen.at(index_d0_aw).PID()%2==0)
	{
	  index_gen[0] = index_d0_aw;
	  index_gen[1] = index_d1_aw;
	}
      else
	{
	  index_gen[0] = index_d1_aw;
          index_gen[1] = index_d0_aw;
	}

      selected_w = -1;//mean W- decayed hadronically
    }
  //no W decays hadronically
  else 
    {
      if(run_debug) cout << "No W decays hadronically" << endl;

      return -999;
    }
   
  //matching with GenHFHadronMatcher information
  for(unsigned int i=0; i<vec_jet.size(); i++)
    {
      int hf_flavour = vec_hf_flavour.at(i);
      int hf_origin = vec_hf_origin.at(i); 
      
      Jet jet = vec_jet.at(i);
      
      float jet_eta = jet.Eta();
      float jet_phi = jet.Phi();
      
      //w is tagged by GenHFHadronMatcher
      if(Abs(hf_origin)==target_pdg_id_positive)
	{
	  for(int j=0; j<2; j++)
	    {
	      if(Abs(hf_flavour)==Abs(vec_gen.at(index_gen[j]).PID()))
		{
		  index_matched_jet[j] = i;
		  surely_matched[j] = true;

		  Gen gen = vec_gen.at(index_gen[j]);

		  float gen_eta = gen.Eta();
		  float gen_phi = gen.Phi();
		  
		  float del_eta = jet_eta - gen_eta;
		  float del_phi = jet_phi - gen_phi; 

		  float del_r = Sqrt(Power(del_eta, 2.)+Power(del_phi, 2.));
		  
		  dr_return[j] = del_r;
		} 
	    }
	}
    }
 
  // cout << endl;
  // cout << "Daughters from W"; 
  // if(selected_w>0) cout << "+" << endl;
  // else cout << "-" << endl;
  // cout << "(" << index_gen[0] << ", " << vec_gen.at(index_gen[0]).PID() <<  "), (" << index_gen[1] << ", " <<  vec_gen.at(index_gen[1]).PID() << ")" << endl;
  
  // for(unsigned int i=0; i<vec_jet.size(); i++)
  //   {
  //     int jet_flavour = vec_hf_flavour.at(i);
  //     int jet_origin = vec_hf_origin.at(i);
      
  //     cout << "(" << jet_flavour << "," <<  jet_origin << "), ";
  //   }
  // cout << endl;
  // for(int i=0; i<2; i++) cout << "Matched jet index = " << index_matched_jet[i] << ", surely matched = " << surely_matched[i] << ", d_r = " << dr_return[i] << endl;
  
  return selected_w;
}//int Vcb_Mu::Gen_Match_W(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour , const vector<int>& vec_hf_origin, int index_matched_jet[2])

//////////

int Vcb_Mu::Get_W_Decay_Mode(const vector<Gen>& vec_gen)
{
  int index_last_w = -999;
  int index_last_aw = -999;
  int index_d0_w = -999;
  int index_d1_w = -999;
  int index_d0_aw = -999;
  int index_d1_aw = -999;

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
      //if(run_debug) cout << "Can't find both of W" << endl;
      return 999;
    }

  int index_had_w[2];
  //W+ decay hadronically
  if(Abs(vec_gen.at(index_d0_w).PID())<10)
    {
      index_had_w[0] = index_d0_w;
      index_had_w[1] = index_d1_w;
    }
  //W- decay hadronically
  else if(Abs(vec_gen.at(index_d0_aw).PID()<10))
    {
      index_had_w[0] = index_d0_aw;
      index_had_w[1] = index_d1_aw;
    }
  //no W decays hadronically
  else
    {
      //if(run_debug) cout << "No W decays hadronically" << endl;

      return 999;
    }

  int decay_mode = 0;
  for(unsigned int i=0; i<2; i++)
    {
      Gen gen = vec_gen.at(index_had_w[i]);

      int pid = gen.PID();
    
      if(pid%2==0) decay_mode += Abs(pid)*10;
      else decay_mode += Abs(pid);
    }
 
  return decay_mode;
}//int Vcb_Mu::Get_W_Decay_Mode(const vector<Gen>& vec_gen)

//////////

void Vcb_Mu::Index_Converter(const vector<Jet>& vec_sel_jet, const vector<Jet>& vec_sel_jet_match, const int index_matched_jet_match[4], int index_matched_jet[4])
{
  for(int i=0; i<4; i++)
    {
      index_matched_jet[i] = -999;

      if(index_matched_jet_match[i]<0) continue;
      
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

void Vcb_Mu::Sol_Neutrino_Pz(const Particle& lepton, const Particle& met, float neutrino_pz_sol[2])
{
  float lepton_mass = lepton.M();

  float met_px = met.Px();
  float met_py = met.Py();
  float met_pt = met.Pt();
  
  double k = TMath::Power(W_MASS, 2.)/2.0 - lepton_mass*lepton_mass/2.0 + lepton.Px()*met_px + lepton.Py()*met_py;
  double a = TMath::Power(lepton.Px(), 2.0) + TMath::Power(lepton.Py(), 2.0);
  double b = -2*k*lepton.Pz();
  double c = TMath::Power(lepton.Pt(), 2.0)*TMath::Power(met_pt, 2.0) - TMath::Power(k, 2.0);

  double determinant = TMath::Power(b, 2.0) - 4*a*c;

  //real solution
  if(determinant>=0)
    {
      neutrino_pz_sol[0] = (-b + TMath::Sqrt(determinant))/(2*a);
      neutrino_pz_sol[1] = (-b - TMath::Sqrt(determinant))/(2*a);
    }
  //complex solution
  else
    {
      neutrino_pz_sol[0] = -b/(2*a);
      //Resol_Neutrino_Pt();
    }

  return;
}//void Vcb_Mu::Sol_Neutrino_Pz(const Particle& lepton, const Particle& met)

//////////
