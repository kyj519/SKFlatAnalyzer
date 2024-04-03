#ifndef __Vcb_Tagging_RF_h__
#define __Vcb_Tagging_RF_h__

#include <TMVA/Reader.h>

#include "AnalyzerCore.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TKinFitterDriver.h"
#include "Vcb_Def.h"
#include "Results_Container.h"

#include "XYMETCorrection_withUL17andUL18andUL16.h"

using namespace std;
using namespace TMath;

class Vcb_Tagging_RF : public AnalyzerCore
{
public:
  Vcb_Tagging_RF();
  ~Vcb_Tagging_RF();

  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter(AnalyzerParameter param);

protected:
  vector<AnalyzerParameter::Syst> vec_syst_type;

  bool run_mu_ch;
  bool run_el_ch;
  bool run_debug;

  vector<TString> vec_mu_id;
  vector<TString> vec_mu_id_sf_key;
  vector<TString> vec_mu_iso_sf_key;
  vector<TString> vec_mu_trig;
  TString mu_trig;
  float mu_trig_safe_pt_cut;

  vector<TString> vec_el_id;
  vector<TString> vec_el_id_sf_key;
  vector<TString> vec_el_trig;
  TString el_trig;
  float el_trig_safe_pt_cut;

  vector<TString> vec_sl_trig; // single lepton
  TString sl_trig;
  float sl_trig_safe_pt_cut;

  vector<JetTagging::Parameters> vec_jet_tagging_para;

  AnalyzerParameter param;

  vector<Muon> vec_muon;
  vector<Electron> vec_electron;
  vector<Jet> vec_jet;

  vector<Muon> vec_this_muon;
  vector<Electron> vec_this_electron;
  vector<Jet> vec_this_jet;

  vector<Gen> vec_gen;

  Muon muon;
  Electron electron;
  Lepton lepton;

  vector<Jet> vec_sel_jet;
  vector<bool> vec_btag;
  vector<bool> vec_ctag;

  Particle met;
  float met_pt;
  float met_phi;

  XYMETCorrection_withUL17andUL18andUL16 xy_met_correction;

  int n_pv;

  float lepton_pt;
  float lepton_eta;

  int n_sel_jet;
  int n_b_jet;
  int n_c_jet;

  float ht;

  float leading_jet_bvsc;
  float leading_jet_cvsb;
  float leading_jet_cvsl;
  float leading_jet_eta;
  float leading_jet_pt;

  float subleading_jet_bvsc;
  float subleading_jet_cvsb;
  float subleading_jet_cvsl;
  float subleading_jet_eta;
  float subleading_jet_pt;

  int decay_mode;

  vector<int> vec_gen_hf_flavour;
  vector<int> vec_gen_hf_origin;

  vector<int> vec_sel_gen_hf_flavour;
  vector<int> vec_sel_gen_hf_origin;

  float weight;

  float weight_hem_veto;
  float weight_lumi;
  float weight_mc;

  float weight_pileup;
  float weight_pileup_down;
  float weight_pileup_up;

  float weight_ps[4];

  float weight_prefire;
  float weight_top_pt;
  float weight_pujet_veto;

  float weight_scale_variation_1;
  float weight_scale_variation_2;
  float weight_scale_variation_3;
  float weight_scale_variation_4;
  float weight_scale_variation_6;
  float weight_scale_variation_8;

  float weight_b_tag;
  float weight_b_tag_down_hf;
  float weight_b_tag_up_hf;
  float weight_b_tag_down_lf;
  float weight_b_tag_up_lf;
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

  float sf_mu_trig;
  float sf_mu_id;
  float sf_mu_iso;

  float sf_el_id;
  float sf_el_reco;

  float sf_sl_trig;

  map<AnalyzerParameter::Syst, TTree *> map_result_tree;

  float Calculate_HT(const vector<Jet> &jet);
  void Clear();
  Particle Rebalance_Met();
  void Set_Result_Tree();
};

#endif /* __Vcb_Tagging_RF_h__ */
