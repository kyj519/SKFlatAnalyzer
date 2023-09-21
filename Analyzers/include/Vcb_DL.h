#ifndef __Vcb_DL_h__
#define __Vcb_DL_h__

#include <TMVA/Reader.h>

#include "AnalyzerCore.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TKinFitterDriver.h"
#include "Vcb_Def.h"
#include "Results_Container.h"

#include "XYMETCorrection_withUL17andUL18andUL16.h"

using namespace std;
using namespace TMath;

class Vcb_DL : public AnalyzerCore
{
public:
  Vcb_DL();
  ~Vcb_DL();

  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter(AnalyzerParameter param);

protected:
  vector<AnalyzerParameter::Syst> vec_syst_type;

  bool run_debug;
  bool run_ee;
  bool run_me;
  bool run_mm;
  bool run_syst;

  vector<TString> vec_mu_id;
  vector<TString> vec_mu_id_sf_key;
  vector<TString> vec_mu_iso_sf_key;

  // vector<TString> vec_el_id;
  // vector<TString> vec_el_id_sf_key;

  vector<TString> vec_el_trig;
  vector<TString> vec_mu_trig;
  vector<TString> vec_sl_trig;

  TString el_trig;
  TString mu_trig;

  float el_trig_safe_pt_cut;
  float mu_trig_safe_pt_cut;

  vector<JetTagging::Parameters> vec_jet_tagging_para;

  vector<Electron> vec_electron;
  vector<Muon> vec_muon;
  vector<Jet> vec_jet;
  vector<Gen> vec_gen;

  vector<Electron> vec_this_electron;
  vector<Muon> vec_this_muon;
  vector<Jet> vec_this_jet;
  Particle met;

  vector<Electron> vec_sel_electron;
  vector<Muon> vec_sel_muon;
  vector<Jet> vec_sel_jet;

  Lepton lepton[2];

  AnalyzerParameter param;

  float leading_lepton_eta;
  float leading_lepton_pt;

  float subleading_lepton_eta;
  float subleading_lepton_pt;

  float dilepton_mass;

  int n_sel_jet;
  int n_b_jet;
  float ht;

  float leading_jet_bvsc;
  float leading_jet_cvsb;
  float leading_jet_cvsl;

  float subleading_jet_bvsc;
  float subleading_jet_cvsb;
  float subleading_jet_cvsl;

  float met_pt;
  float met_phi;

  float bvsc_third;
  float bvsc_fourth;

  int decay_mode;

  vector<int> vec_gen_hf_flavour;
  vector<int> vec_gen_hf_origin;

  vector<int> vec_sel_gen_hf_flavour;
  vector<int> vec_sel_gen_hf_origin;

  float weight;

  float weight_b_tag;
  float weight_b_tag_down_hf;
  float weight_b_tag_up_hf;
  float weight_b_tag_down_jes;
  float weight_b_tag_up_jes;
  float weight_b_tag_down_lfstats1;
  float weight_b_tag_up_lfstats1;
  float weight_b_tag_down_lfstats2;
  float weight_b_tag_up_lfstats2;
  float weight_b_tag_down_cferr1;
  float weight_b_tag_up_cferr1;
  float weight_b_tag_down_cferr2;
  float weight_b_tag_up_cferr2;
  float weight_b_tag_down_hfstats1;
  float weight_b_tag_up_hfstats1;
  float weight_b_tag_down_hfstats2;
  float weight_b_tag_up_hfstats2;

  float weight_c_tag;
  float weight_c_tag_down_extrap;
  float weight_c_tag_up_extrap;
  float weight_c_tag_down_interp;
  float weight_c_tag_up_interp;
  float weight_c_tag_down_lhe_scale_muf;
  float weight_c_tag_up_lhe_scale_muf;
  float weight_c_tag_down_lhe_scale_mur;
  float weight_c_tag_up_lhe_scale_mur;
  float weight_c_tag_down_ps_fsr_fixed;
  float weight_c_tag_up_ps_fsr_fixed;
  float weight_c_tag_down_ps_isr_fixed;
  float weight_c_tag_up_ps_isr_fixed;
  float weight_c_tag_down_pu;
  float weight_c_tag_up_pu;
  float weight_c_tag_down_stat;
  float weight_c_tag_up_stat;
  float weight_c_tag_down_xsec_brunc_dyjets_b;
  float weight_c_tag_up_xsec_brunc_dyjets_b;
  float weight_c_tag_down_xsec_brunc_dyjets_c;
  float weight_c_tag_up_xsec_brunc_dyjets_c;
  float weight_c_tag_down_xsec_brunc_wjets_c;
  float weight_c_tag_up_xsec_brunc_wjets_c;
  float weight_c_tag_down_jer;
  float weight_c_tag_up_jer;
  float weight_c_tag_down_jes_total;
  float weight_c_tag_up_jes_total;

  float weight_el_id;
  float weight_el_id_down;
  float weight_el_id_up;

  float weight_el_reco;
  float weight_el_reco_down;
  float weight_el_reco_up;

  float weight_hem_veto;
  float weight_lumi;
  float weight_mc;

  float weight_mu_id;
  float weight_mu_id_down;
  float weight_mu_id_up;

  float weight_mu_iso;
  float weight_mu_iso_down;
  float weight_mu_iso_up;

  float weight_pdf_alternative;
  float weight_pdf_error_set[100];
  float weight_pdf_as_down;
  float weight_pdf_as_up;

  float weight_pileup;
  float weight_pileup_down;
  float weight_pileup_up;

  float weight_prefire;
  float weight_prefire_down;
  float weight_prefire_up;

  float weight_ps[4];

  float weight_pujet_veto;
  float weight_pujet_veto_down;
  float weight_pujet_veto_up;

  float weight_scale_variation_1;
  float weight_scale_variation_2;
  float weight_scale_variation_3;
  float weight_scale_variation_4;
  float weight_scale_variation_6;
  float weight_scale_variation_8;

  float weight_sl_trig;
  float weight_sl_trig_el_down;
  float weight_sl_trig_el_up;
  float weight_sl_trig_mu_down;
  float weight_sl_trig_mu_up;

  float weight_top_pt;

  map<AnalyzerParameter::Syst, TTree *> map_result_tree;

  XYMETCorrection_withUL17andUL18andUL16 xy_met_correction;

  float Calculate_HT(const vector<Jet> &vec_jet);
  int Check_Process(const vector<Gen> &vec_gen);
  void Clear();
  void Make_Result_Tree(const AnalyzerParameter &param);
  Particle Rebalance_Met();
  void Set_Result_Tree();
};

#endif /* __Vcb_DL_h__ */
