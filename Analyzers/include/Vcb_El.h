#ifndef __Vcb_El_h__
#define __Vcb_El_h__

#include "AnalyzerCore.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TKinFitterDriver.h"
#include "Vcb_Def.h"
#include "Results_Container.h"

using namespace std;

class Vcb_El : public AnalyzerCore
{
 public:
  Vcb_El();
  ~Vcb_El();

  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter(AnalyzerParameter param);
  
 protected:
  bool run_syst;

  vector<TString> vec_el_id;
  vector<TString> vec_el_id_sf_key;
  vector<TString> vec_el_trig;
  float trig_safe_pt_cut;

  vector<JetTagging::Parameters> vec_jet_tagging_para;

  JME::JetResolution jet_resolution;
  JME::JetResolutionScaleFactor jet_resolution_sf;

  TKinFitterDriver* fitter_driver;

  vector<Electron> vec_electron;
  vector<Muon> vec_muon;
  vector<Jet> vec_jet;
};

#endif /* __Vcb_El_h__ */
