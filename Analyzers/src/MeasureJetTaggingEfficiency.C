#include "MeasureJetTaggingEfficiency.h"

//////////

MeasureJetTaggingEfficiency::MeasureJetTaggingEfficiency()
{
}//MeasureJetTaggingEfficiency::MeasureJetTaggingEfficiency

//////////

MeasureJetTaggingEfficiency::~MeasureJetTaggingEfficiency()
{
}//MeasureJetTaggingEfficiency::~MeasureJetTaggingEfficiency()

//////////

void MeasureJetTaggingEfficiency::initializeAnalyzer()
{
  TString datapath = getenv("DATA_DIR");
  TString btagpath = datapath+"/"+DataEra+"/BTag/";

  Taggers.clear();
  WPs.clear();
  cut_values.clear();

  ifstream in_tagger(btagpath+"/CutValues.txt");
  string btaggerline;
  while(getline(in_tagger,btaggerline))
    {
      std::istringstream is_tag( btaggerline );
      TString tstring_taggerline = btaggerline;
      if(tstring_taggerline.Contains("#")) continue;
      TString a;
      string b,c;
      float d;
      float e;
      
      is_tag >> a; // ERA
      is_tag >> b; // TAGGER
      is_tag >> c; // WP
      is_tag >> d; // cut value
      is_tag >> e;// cut value for CvsB, dummy for B-tagging

      if(a == DataEra)
	{
	  Taggers.push_back(b);
	  WPs.push_back(c);
	  cut_values.push_back(make_pair(d,e));
	}
    }// end of taggermap loop
  
  cout << "[MeasureJetTaggingEfficiency::initializeAnalyzer()] What to measure :" << endl;
  cout << "[MeasureJetTaggingEfficiency::initializeAnalyzer()] Tagger\tWP\tCutValue" << endl;
  for(unsigned i_m=0; i_m<Taggers.size(); i_m++)
    {
      string Tagger = Taggers.at(i_m);
      string WP = WPs.at(i_m);
      auto cut_pair = cut_values.at(i_m);
 
      cout << "[MeasureJetTaggingEfficiency::initializeAnalyzer()] " << Tagger << "\t" << WP << "\t" << cut_pair.first << "\t" << cut_pair.second << endl;
    }

  vec_etabins = {0.0, 0.8, 1.6, 2., 2.5};
  vec_ptbins = {20., 30., 50., 70., 100., 140., 200., 300., 600., 1000.};//PT bins used in POG SF measurements
  //for average users, this binning will be sufficient.
  //but eta-dependence of efficiency can be larger for |eta|>~2, where track & muon detector information of jet constituents starts to get lost, which is critical in tagging.
  //precision analysis with high-eta b may use finer binnings there, but beware of small number of b-jets in high-eta, high-pt bins if you use ttbar sample; proper optimization of bin size should be studied.
  
  PtMax = vec_ptbins.at( vec_ptbins.size()-1 );
  NEtaBin = vec_etabins.size()-1;
  NPtBin = vec_ptbins.size()-1;

  etabins = new double[NEtaBin+1];
  for(int i=0; i<NEtaBin+1; i++) etabins[i] = vec_etabins.at(i);
  ptbins = new double[NPtBin+1];
  for(int i=0; i<NPtBin+1; i++) ptbins[i] = vec_ptbins.at(i);
}//void MeasureJetTaggingEfficiency::initializeAnalyzer

//////////

void MeasureJetTaggingEfficiency::executeEvent()
{
  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
  Particle METv = ev.GetMETVector();

  vector<Jet> jets = GetJets("tightLepVeto", 20., 2.5);
  float weight = 1.;
  float w_Gen  = MCweight();
  float w_Norm = ev.GetTriggerLumi("Full");
  float w_PU   = GetPileUpWeight(nPileUp, 0);
  weight *= w_Gen*w_Norm*w_PU; 
  //tagging performance depends on PU, so it is better reweight to proper PU profile
  
  //==== code to measure btag efficiencies in TT MC
  //==== Reference : https://github.com/rappoccio/usercode/blob/Dev_53x/EDSHyFT/plugins/BTaggingEffAnalyzer.cc
  for(unsigned int ij = 0 ; ij < jets.size(); ij++)
    {
      TString flav= "B";
      if(fabs(jets.at(ij).hadronFlavour()) == 4) flav= "C";
      if(fabs(jets.at(ij).hadronFlavour()) == 0) flav= "Light";
      
      double this_Eta = fabs(jets.at(ij).Eta());//POG recommendation is to use |eta|
      double this_Pt = jets.at(ij).Pt()<PtMax ? jets.at(ij).Pt() : PtMax-1; // put overflows in the last bin

      //==== First, fill the denominator
      FillHist("Jet_"+DataEra+"_eff_"+flav+"_denom", this_Eta, this_Pt, weight, NEtaBin, etabins, NPtBin, ptbins);

      //==== Now looping over (tagger,working point)
      for(unsigned i_m=0; i_m<Taggers.size(); i_m++)
	{ 
	  
	  string Tagger = Taggers.at(i_m);
	  string WP = WPs.at(i_m);
	  auto cut_pair = cut_values.at(i_m);

	  if(Tagger.find("_C")==string::npos)
	    {
	      double this_taggerresult = jets.at(ij).GetTaggerResult(JetTagging::StringToTagger(Tagger));
	  
	      if(cut_pair.first<this_taggerresult) FillHist("Jet_"+DataEra+"_"+Tagger+"_"+WP+"_eff_"+flav+"_num", this_Eta, this_Pt, weight, NEtaBin, etabins, NPtBin, ptbins);
	    }
	  else 
	    {
	      double this_cvsl = jets.at(ij).GetTaggerResult(JetTagging::StringToTagger(Tagger, 0));
	      double this_cvsb = jets.at(ij).GetTaggerResult(JetTagging::StringToTagger(Tagger, 1));

	      if(cut_pair.first<this_cvsl && cut_pair.second<this_cvsb) FillHist("Jet_"+DataEra+"_"+Tagger+"_"+WP+"_eff_"+flav+"_num", this_Eta, this_Pt, weight, NEtaBin, etabins, NPtBin, ptbins);
	    }
	} // END Loop (tagger,working point)
    } // END Loop jet

  return;
}//void MeasureJetTaggingEfficiency::executeEvent
