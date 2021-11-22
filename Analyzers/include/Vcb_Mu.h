#ifndef __Vcb_Mu_h__
#define __Vcb_Mu_h__

#include "AnalyzerCore.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TKinFitterDriver.h"
#include "Vcb_Def.h"
#include "Results_Container.h"

using namespace std;

class Vcb_Mu : public AnalyzerCore 
{
 public:
  Vcb_Mu();
  ~Vcb_Mu();

  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter(AnalyzerParameter param);
    
 protected:
  bool run_syst;

  vector<TString> vec_mu_id;
  vector<TString> vec_mu_id_sf_key;
  vector<TString> vec_mu_iso_sf_key;
  vector<TString> vec_mu_trig;
  float trig_safe_pt_cut;

  vector<JetTagging::Parameters> vec_jet_tagging_para;

  JME::JetResolution jet_resolution;
  JME::JetResolutionScaleFactor jet_resolution_sf;

  TKinFitterDriver* fitter_driver;

  vector<Muon> vec_muon;
  vector<Electron> vec_electron;
  vector<Jet> vec_jet;
  
  float weight_prefire;

  float Get_CvsB(const Jet& jet);
  float Get_CvsL(const Jet& jet);
};

#endif /* __Vcb_Mu_h__ */

