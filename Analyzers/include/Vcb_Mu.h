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
  bool run_syst;
  bool rm_wm_constraint;

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
  
  float weight_prefire;

  float weight;

  float chi2;
  
  float bvsc_w_u;
  float cvsb_w_u;
  float cvsl_w_u;

  float bvsc_w_d;
  float cvsb_w_d;
  float cvsl_w_d;

  TTree* result_tree;

  int Chk_Included(const vector<Jet>& vec_sel_jet, const vector<Jet>& vec_sel_jet_matched, const int index_matched_jet[4]);
  bool Compare_Jet_Pair(const Jet jet0[2], const Jet jet1[2]);
  void Gen_Match_TT(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, const vector<float>& vec_jer, int index_gen[4], int index_matched_jet[4], bool surely_matched[2], float dr_return[2]);
  int Gen_Match_W(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, const vector<float>& vec_jer, int index_gen[2], int index_matched_jet[2], bool surely_matched[2], float dr_return[2]);
  float Get_CvsB(const Jet& jet);
  float Get_CvsL(const Jet& jet);
};

#endif /* __Vcb_Mu_h__ */

