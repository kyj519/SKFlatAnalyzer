#ifndef MeasureJetTaggingEfficiency_h
#define MeasureJetTaggingEfficiency_h

#include "AnalyzerCore.h"


class MeasureJetTaggingEfficiency : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEvent();

  MeasureJetTaggingEfficiency();
  ~MeasureJetTaggingEfficiency();

 protected:
  vector<double> vec_etabins;
  vector<double> vec_ptbins;

  double PtMax;
  
  int NEtaBin;
  int NPtBin;

  double* etabins;
  double* ptbins;
    
  //==== Read what to measrue from data/Run2Legacy_v4/<Year>/BTag/taggermap.txt
  vector<string> Taggers;
  vector<string> WPs;
  vector<pair<double, double>> cut_values;
};



#endif

