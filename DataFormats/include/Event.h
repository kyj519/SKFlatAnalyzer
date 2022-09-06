#ifndef Event_h
#define Event_h

#include "Particle.h"
#include "TString.h"
#include "TObject.h"

class Event : public TObject {
public:

  Event();
  ~Event();

  void SetnPV(double n);
  inline int nPV() const { return j_nPV; }

  void SetTrigger(vector<string> HLT_TriggerName);
  bool PassTrigger(TString trig);
  bool PassTrigger(std::vector<TString> trigs);
  double GetTriggerLumi(TString trig);
  bool IsPDForTrigger(TString trig, TString PD);

  void SetMET(const TString& type, double pt, double phi);
  Particle GetMETVector(const TString& type="PF") const;

  void SetEra(TString era){
    j_DataEra=era;
    j_DataYear=TString(era(0,4)).Atoi();
  }
  TString GetEra() const { return j_DataEra; }
  int GetYear() const { return j_DataYear; }


private:
  int j_nPV;
  vector<string> j_HLT_TriggerName;
  Particle j_METVector_pf;
  Particle j_METVector_puppi;
  int j_DataYear;
  TString j_DataEra;

  ClassDef(Event,1)
};

#endif
