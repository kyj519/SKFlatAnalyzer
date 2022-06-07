#ifndef __Vcb_Mu_h__
#define __Vcb_Mu_h__

#include "AnalyzerCore.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TKinFitterDriver.h"
#include "Vcb_Def.h"
#include "Results_Container.h"

using namespace std;
using namespace TMath;

class Vcb_Mu : public AnalyzerCore 
{
 public:
  Vcb_Mu();
  ~Vcb_Mu();

  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter(AnalyzerParameter param);
    
 protected:
  bool run_debug;
  bool run_kf_eval;
  bool run_permutation_tree;
  bool run_syst;
  bool rm_wm_constraint;

  typedef enum cut_flow {No_Cut, Met_Filter, Trigger, TT_Selection, KF_Pass} Cut_Flow;
  int n_cut_flow = Cut_Flow::KF_Pass+1;

  vector<TString> vec_mu_id;
  vector<TString> vec_mu_id_sf_key;
  vector<TString> vec_mu_iso_sf_key;
  vector<TString> vec_mu_trig;
  TString mu_trig;
  float trig_safe_pt_cut;

  vector<JetTagging::Parameters> vec_jet_tagging_para;

  JME::JetResolution jet_resolution;
  JME::JetResolutionScaleFactor jet_resolution_sf;

  TKinFitterDriver* fitter_driver;

  vector<Muon> vec_muon;
  vector<Electron> vec_electron;
  vector<Jet> vec_jet;
  
  int n_sel_jet;
  int n_jet_bin;
  int n_bjet_bin;

  float weight_prefire;

  float weight;
    
  bool chk_matched_jets_only;
  bool chk_kf_correct;

  float had_t_b_pt;
  float w_u_pt;
  float w_d_pt;
  float lep_t_b_pt;
  
  float had_t_b_bscore;
  float lep_t_b_bscore;
  
  float del_phi_w_u_w_d;
  float del_phi_had_w_had_t_b;
  float del_phi_lep_neu;
  float del_phi_lep_w_lep_t_b;
  float del_phi_had_t_lep_t;
  
  float del_eta_w_u_w_d;
  float del_eta_had_w_had_t_b;
  float del_eta_lep_neu;
  float del_eta_lep_w_lep_t_b;
  float del_eta_had_t_lep_t;

  float del_r_w_u_w_d;
  float del_r_had_w_had_t_b;
  float del_r_lep_neu;
  float del_r_lep_w_lep_t_b;
  float del_r_had_t_lep_t;

  float theta_w_u_w_d;
  float theta_had_w_had_t_b;
  float theta_lep_neu;
  float theta_lep_w_lep_t_b;
  float theta_had_t_lep_t;
  
  float had_t_mass;
  float had_w_mass;
  float lep_t_mass;
  float lep_t_partial_mass;

  float chi2_jet_had_t_b;
  float chi2_jet_w_u;
  float chi2_jet_w_d;
  float chi2_jet_lep_t_b;

  float chi2_jet_extra;

  float chi2_constraint_had_t;
  float chi2_constraint_had_w;
  float chi2_constraint_lep_t;
  float chi2_constraint_lep_w;
 
  float chi2;

  float gen_neutrino_px;
  float gen_neutrino_py;
  float gen_neutrino_pz;
  float met_px;
  float met_py;
  float met_rebalance_px;
  float met_rebalance_py;
  float neutrino_pz_sol;
  float neutrino_pz_sol_unrebal;
  float mt_gen;
  float mt_met;
  float mt_met_rebalance;

  bool chk_jet_permutation_match;
  bool chk_real_neu_pz;
  float nu_pz_sol_0;
  float nu_pz_sol_1;

  float bvsc_w_u;
  float cvsb_w_u;
  float cvsl_w_u;

  float bvsc_w_d;
  float cvsb_w_d;
  float cvsl_w_d;

  TTree* result_tree;
  TTree* event_tree;
  TTree* permutation_tree_correct;
  TTree* permutation_tree_wrong;
 
  float Calculate_Mt(const Particle& lepton, const float& neu_px, const float& neu_py);
  int Chk_Included(const int index_matched_jet[4]);
  bool Compare_Jet_Pair(const Jet jet0[2], const Jet jet1[2]);
  void Gen_Match_TT(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, const vector<float>& vec_jer, int index_gen[4], int index_matched_jet[4], bool surely_matched[2], float dr_return[2]);
  int Gen_Match_W(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, const vector<float>& vec_jer, int index_gen[2], int index_matched_jet[2], bool surely_matched[2], float dr_return[2]);
  float Get_CvsB(const Jet& jet);
  float Get_CvsL(const Jet& jet);
  void Index_Converter(const vector<Jet>& vec_sel_jet, const vector<Jet>& vec_sel_jet_match, const int index_matched_jet_match[4], int index_matched_jet[4]);
  //void Index_Restorer(int& index_had_t_b, int& index_w_u, int& index_w_d, int& index_lep_t_b);
  Gen Neutrino(const vector<Gen>& vec_gen);
  void Sol_Neutrino_Pz(const Particle& lepton, const Particle& met, float neutrino_pz_sol[2]);
};

#endif /* __Vcb_Mu_h__ */

