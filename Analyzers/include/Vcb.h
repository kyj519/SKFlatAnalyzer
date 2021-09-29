#ifndef Vcb_h
#define Vcb_h

#include "AnalyzerCore.h"

#include "Define_Def.h"

using namespace std;

class Vcb : public AnalyzerCore {
  
 public:
  
  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();
  
  Vcb();
  ~Vcb();
  
 protected:
  vector<TString> vec_muid, vec_muid_sf_key;
  TString iso_mu_trig_name;
  float trig_safe_pt_cut;

  vector<Muon> vec_muon;
  vector<Jet> vec_jet;

};



#endif

