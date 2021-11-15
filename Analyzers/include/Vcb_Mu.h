#ifndef Vcb_Mu_h
#define Vcb_Mu_h

#include "AnalyzerCore.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TKinFitterDriver.h"
#include "Vcb_Def.h"
#include "Results_Container.h"

using namespace std;

class Vcb_Mu : public AnalyzerCore 
{
 public:
  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();
  
  Vcb_Mu();
  ~Vcb_Mu();
  
 protected:
  bool run_syst;

  vector<TString> vec_mu_id;
  vector<TString> vec_mu_id_sf_key;
  vector<TString> vec_mu_iso_sf_key;
  vector<TString> iso_mu_trig_name;
  float trig_safe_pt_cut;

  vector<JetTagging::Parameters> vec_jet_tagging_para;

  vector<Muon> vec_muon;
  vector<Jet> vec_jet;
  vector<Electron> vec_electron;

  JME::JetResolution jet_resolution;
  JME::JetResolutionScaleFactor jet_resolution_sf;

  float weight_prefire;

  TKinFitterDriver* fitter_driver;

  float Get_CvsB(const Jet& jet);
  float Get_CvsL(const Jet& jet);
};



#endif

