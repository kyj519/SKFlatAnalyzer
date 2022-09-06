#ifndef JetTaggingParameters_h
#define JetTaggingParameters_h

#include <iostream>

using namespace std;

namespace JetTagging{

  enum Tagger{
    DeepCSV, DeepCSV_C, DeepCSV_CvsL, DeepCSV_CvsB,
    DeepJet, DeepJet_C, DeepJet_CvsL, DeepJet_CvsB,
    particleNet_TvsQCD, particleNet_WvsQCD, particleNet_ZvsQCD, 
    particleNet_HbbvsQCD, particleNet_HccvsQCD, particleNet_H4qvsQCD, particleNet_QCD, 
    particleNetMD_Xbb, particleNetMD_Xcc, particleNetMD_Xqq, particleNetMD_QCD,
  };
  inline string TaggerToString(Tagger t){
    if(t == DeepCSV) return "DeepCSV";
    if(t == DeepCSV_C) return "DeepCSV_C";
    if(t == DeepCSV_CvsL) return "DeepCSV_CvsL";
    if(t == DeepCSV_CvsB) return "DeepCSV_CvsB";
    if(t == DeepJet) return "DeepJet";
    if(t == DeepJet_C) return "DeepJet_C";
    if(t == DeepJet_CvsL) return "DeepJet_CvsL";
    if(t == DeepJet_CvsB) return "DeepJet_CvsB";
    if(t == particleNet_TvsQCD) return "particleNet_TvsQCD";
    if(t == particleNet_WvsQCD) return "particleNet_WvsQCD";
    if(t == particleNet_ZvsQCD) return "particleNet_ZvsQCD";
    if(t == particleNet_HbbvsQCD) return "particleNet_HbbvsQCD";
    if(t == particleNet_HccvsQCD) return "particleNet_HccvsQCD";
    if(t == particleNet_H4qvsQCD) return "particleNet_H4qvsQCD";
    if(t == particleNet_QCD) return "particleNet_QCD";
    if(t == particleNetMD_Xbb) return "particleNetMD_Xbb";
    if(t == particleNetMD_Xcc) return "particleNetMD_Xcc";
    if(t == particleNetMD_Xqq) return "particleNetMD_Xqq";
    if(t == particleNetMD_QCD) return "particleNetMD_QCD";
    
    cerr << "[JetTaggingParameters::TaggerToString] Wrong tagger enum : " << t << endl;
    exit(ENODATA);

    return "None";
  }
  inline Tagger StringToTagger(string s){
    if(s == "DeepCSV") return DeepCSV;
    if(s == "DeepCSV_CvsL") return DeepCSV_CvsL;
    if(s == "DeepCSV_CvsB") return DeepCSV_CvsB;
    if(s == "DeepJet") return DeepJet;
    if(s == "DeepJet_CvsL") return DeepJet_CvsL;
    if(s == "DeepJet_CvsB") return DeepJet_CvsB;
    if(s == "particleNet_TvsQCD") return particleNet_TvsQCD;
    if(s == "particleNet_WvsQCD") return particleNet_WvsQCD;
    if(s == "particleNet_ZvsQCD") return particleNet_ZvsQCD;
    if(s == "particleNet_HbbvsQCD") return particleNet_HbbvsQCD;
    if(s == "particleNet_HccvsQCD") return particleNet_HccvsQCD;
    if(s == "particleNet_H4qvsQCD") return particleNet_H4qvsQCD;
    if(s == "particleNet_QCD") return particleNet_QCD;
    if(s == "particleNetMD_Xbb") return particleNetMD_Xbb;
    if(s == "particleNetMD_Xcc") return particleNetMD_Xcc;
    if(s == "particleNetMD_Xqq") return particleNetMD_Xqq;
    if(s == "particleNetMD_QCD") return particleNetMD_QCD;

    cerr << "[JetTaggingParameters:: StringToTagger] Wrong tagger string : " << s << endl;
    exit(ENODATA);

    return DeepCSV;
  }

  inline Tagger StringToTagger(string s, int i){
    if(s == "DeepCSV_C"){
      if(i==0) return DeepCSV_CvsL;
      else if(i==1) return DeepCSV_CvsB;
    }
    else if(s == "DeepJet_C"){
      if(i==0) return DeepJet_CvsL;
      else if(i==1) return DeepJet_CvsB;
    }
    
    cerr << "[JetTaggingParameters:: StringToTagger] Wrong tagger string : " << s << ", or wrong index : " << i << endl;
    exit(ENODATA);
    
    return DeepCSV;
  }
  
  enum WP{
    Loose, Medium, Tight
  };

  inline string WPToString(WP wp){
    if(wp == Loose) return "Loose";
    if(wp == Medium) return "Medium";
    if(wp == Tight) return "Tight";

    cerr << "[JetTaggingParameters::WPToString] Wrong WP enum : " << wp << endl;
    exit(ENODATA);

    return "None";
  }

  enum MeasurmentType{
    incl, mujets, comb, wcharm, iterativefit 
  };

  inline string MeasurmentTypeToString(MeasurmentType mt){
    if(mt == incl) return "incl";
    if(mt == mujets) return "mujets";
    if(mt == comb) return "comb";
    if(mt == wcharm) return "wcharm";
    if(mt == iterativefit) return "iterativefit";

    cerr << "[JetTaggingParameters::MeasurmentTypeToString] Wrong MeasurmentType enum : " << mt << endl;
    exit(ENODATA);

    return "None";

  }

  class Parameters{
  public:
    inline Parameters(){

    }
    inline Parameters(Tagger t, WP wp, MeasurmentType mt_l, MeasurmentType mt_h){
      j_Tagger = t;
      j_WP = wp;
      j_MeasurmentType_Light = mt_l;
      j_MeasurmentType_Heavy = mt_h;
    }

    Tagger j_Tagger;
    WP j_WP;
    MeasurmentType j_MeasurmentType_Light, j_MeasurmentType_Heavy;

    inline void Print(){
      printf("(%s,%s,%s,%s)\n",TaggerToString(j_Tagger).c_str(),WPToString(j_WP).c_str(),MeasurmentTypeToString(j_MeasurmentType_Light).c_str(),MeasurmentTypeToString(j_MeasurmentType_Heavy).c_str());
    }

  };

} // namespace

#endif
