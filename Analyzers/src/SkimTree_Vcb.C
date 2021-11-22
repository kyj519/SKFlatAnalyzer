#include "SkimTree_Vcb.h"

//////////

SkimTree_Vcb::SkimTree_Vcb()
{
  newtree = NULL;
}//SkimTree_Vcb::SkimTree_Vcb()

//////////

SkimTree_Vcb::~SkimTree_Vcb()
{
}//SkimTree_Vcb::~SkimTree_Vcb()

//////////

void SkimTree_Vcb::initializeAnalyzer()
{
  outfile->cd();
  cout << "[SkimTree_Vcb::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);
  
  vec_single_muon_trigger.clear();
  vec_single_electron_trigger.clear();

  if(DataYear==2016)
    {
      vec_single_muon_trigger = {"HLT_IsoMu24_v"};
    }
  else if(DataYear==2017)
    {
      vec_single_muon_trigger = {"HLT_IsoMu27_v"};
    }
  else if(DataYear==2018)
    {
      vec_single_muon_trigger = {"HLT_IsoMu24_v"};
    }
  else std::runtime_error("No trigger configuration for year");
    
  cout << "[SkimTree_Vcb::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<vec_single_muon_trigger.size(); i++)  cout << "[SkimTree_Vcb::initializeAnalyzer] " << vec_single_muon_trigger.at(i) << endl;
  for(unsigned int i=0; i<vec_single_electron_trigger.size(); i++)  cout << "[SkimTree_Vcb::initializeAnalyzer] " << vec_single_electron_trigger.at(i) << endl;

  return;
}//void SkimTree_Vcb::initializeAnalyzer()

//////////

void SkimTree_Vcb::executeEvent()
{
  //met filter
  if(!PassMETFilter()) return;
  
  Event ev = GetEvent();
  Particle met = ev.GetMETVector();

  //trigger
  bool chk_mu_trig = ev.PassTrigger(vec_single_muon_trigger);
  bool chk_el_trig = ev.PassTrigger(vec_single_electron_trigger);
  
  if(chk_mu_trig==false&&chk_el_trig==false) return;
  
  vec_muon = GetAllMuons();
  vec_electron = GetAllElectrons();
  vec_jet = GetAllJets();
 
  vector<Muon> vec_sel_muon = SelectMuons(vec_muon, "POGTight", MUON_PT_SKIM, MUON_ETA_SKIM);
  vector<Electron> vec_sel_electron = SelectElectrons(vec_electron, "passTightID", ELECTRON_PT_SKIM, ELECTRON_ETA_SKIM);
  vector<Jet> vec_sel_jet = SelectJets(vec_jet, "tight", JET_PT_SKIM, JET_ETA_SKIM);

  //n of btag
  int nbtag = 0;
  for(auto& jet : vec_sel_jet)
    {
      float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
      if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Loose) < tagging_score) nbtag++;
    }
  
  //pre selection
  if(met.Pt()<MET_PT_SKIM) return;
  if(1<vec_sel_muon.size()) return;
  if(2<vec_sel_electron.size()) return;
  if(vec_sel_jet.size()<4) return;
  if(nbtag<2) return;

  //fill
  if(newtree->Fill()<0) exit(EIO);

  return;
}//void SkimTree_Vcb::executeEvent()

//////////

void SkimTree_Vcb::executeEventFromParameter(AnalyzerParameter param)
{
  return;
}//void SkimTree_Vcb::executeEvent(AnalyzerParameter param)

//////////

void SkimTree_Vcb::WriteHist()
{
  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

  return;
}//void SkimTree_Vcb::WriteHist()

//////////
