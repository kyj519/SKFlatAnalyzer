#include "Jet.h"

ClassImp(Jet)

Jet::Jet() : Particle() {
  j_area=-999.;
  j_partonFlavour=-999;
  j_hadronFlavour=-999;
  j_CSVv2=-999.;
  j_DeepCSV=-999.;
  j_DeepFlavour_b=-999;
  j_DeepFlavour_bb=-999;
  j_DeepFlavour_lepb=-999;
  j_DeepFlavour_c=-999;
  j_DeepFlavour_uds=-999;
  j_DeepFlavour_g=-999;
  j_CvsL=-999.;
  j_CvsB=-999.;
  j_DeepCvsL=-999.;
  j_DeepCvsB=-999.;
  j_chargedHadronEnergyFraction=-999.;
  j_neutralHadronEnergyFraction=-999.;
  j_neutralEmEnergyFraction=-999.;
  j_chargedEmEnergyFraction=-999.;
  j_muonEnergyFraction=-999.;
  j_chargedMultiplicity=-999;
  j_neutralMultiplicity=-999;
  j_PileupJetId=-999.;
  j_En_up=1.;
  j_En_down=1.;;
  j_Res_up = 1.;
  j_Res_down = 1.;
  j_tightJetID=false;
  j_tightLepVetoJetID=false;
}

Jet::~Jet(){

}

void Jet::SetArea(double area){
  j_area = area;
}
void Jet::SetGenFlavours(double pf, double hf){
  j_partonFlavour = pf;
  j_hadronFlavour = hf;
}
void Jet::SetTaggerResults(std::vector<double> ds){
  j_CSVv2             = ds.at(0);
  j_DeepCSV           = ds.at(1);
  j_DeepCvsL          = ds.at(2);
  j_DeepCvsB          = ds.at(3);
  j_DeepFlavour_b     = ds.at(4);
  j_DeepFlavour_bb    = ds.at(5);
  j_DeepFlavour_lepb  = ds.at(6);
  j_DeepFlavour_c     = ds.at(7);
  j_DeepFlavour_uds   = ds.at(8);
  j_DeepFlavour_g     = ds.at(9);
  j_CvsL              = ds.at(10);
  j_CvsB              = ds.at(11);
}
void Jet::SetEnergyFractions(double cH, double nH, double nEM, double cEM, double muE){
  j_chargedHadronEnergyFraction = cH;
  j_neutralHadronEnergyFraction = nH;
  j_neutralEmEnergyFraction = nEM;
  j_chargedEmEnergyFraction = cEM;
  j_muonEnergyFraction = muE;
}
void Jet::SetMultiplicities(double cM, double nM){
  j_chargedMultiplicity = cM;
  j_neutralMultiplicity = nM;
}
void Jet::SetPileupJetId(double v){
  j_PileupJetId = v;
}

void Jet::SetEnShift(double en_up, double en_down){
  j_En_up = en_up;
  j_En_down = en_down;
}

void Jet::SetResShift(double res_up, double res_down){
  j_Res_up = res_up;
  j_Res_down = res_down;
}

void Jet::SetTightJetID(double b){
  j_tightJetID = b;
}
void Jet::SetTightLepVetoJetID(double b){
  j_tightLepVetoJetID = b;
}

bool Jet::PassID(TString ID) const {

  if(ID=="tight") return Pass_tightJetID();
  if(ID=="tightLepVeto") return Pass_tightLepVetoJetID();

  cout << "[Jet::PassID] No id : " << ID << endl;
  exit(ENODATA);

  return false;

}

double Jet::GetTaggerResult(JetTagging::Tagger tg) const {

  if(tg==JetTagging::CSVv2) return j_CSVv2;
  else if(tg==JetTagging::DeepCSV) return j_DeepCSV;
  else if(tg==JetTagging::DeepJet) return j_DeepFlavour_b+j_DeepFlavour_bb+j_DeepFlavour_lepb;
  else if(tg==JetTagging::DeepFlavour_b) return j_DeepFlavour_b;
  else if(tg==JetTagging::DeepFlavour_bb) return j_DeepFlavour_bb;
  else if(tg==JetTagging::DeepFlavour_lepb) return j_DeepFlavour_lepb;
  else if(tg==JetTagging::DeepFlavour_c) return j_DeepFlavour_c;
  else if(tg==JetTagging::DeepFlavour_uds) return j_DeepFlavour_uds;
  else if(tg==JetTagging::DeepFlavour_g) return j_DeepFlavour_g;
  else if(tg==JetTagging::CvsL) return j_CvsL;
  else if(tg==JetTagging::CvsB) return j_CvsB;
  else if(tg==JetTagging::DeepCvsL) return j_DeepCvsL;
  else if(tg==JetTagging::DeepCvsB) return j_DeepCvsB;
  else{
    cout << "[Jet::GetTaggerResult] ERROR; Wrong tagger : " << tg << endl;
    return -999;
  }
}

int Jet::GetPUID(double Pt, double eta) const{
  unsigned PtNum = 0;
  unsigned etaNum = 0;
  if(10.<=Pt && Pt < 20.) PtNum = 1;
  else if(20.<=Pt && Pt < 30.) PtNum = 5;
  else if(30.<=Pt && Pt < 40.) PtNum = 9;
  else if(40.<=Pt && Pt < 50.) PtNum = 13;
  if( fabs(eta) < 2.5) etaNum = 0;
  else if(2.5<=fabs(eta) && fabs(eta) < 2.75) etaNum = 1;
  else if(2.75<=fabs(eta) && fabs(eta) < 3.0) etaNum = 2;
  else if(3.0<=fabs(eta) && fabs(eta) < 5.0) etaNum = 3; 
  unsigned int binID = PtNum + etaNum;
  
  double wp[16][3];
  double temp[48] = {0.77,0.38,-0.31,-0.21,0.90,0.60,-0.12,-0.13,0.96,0.82,0.20,0.09,0.98,0.92,0.47,0.29
                    ,0.26,-0.33,-0.54,-0.37,0.68,-0.04,-0.43,-0.30,0.90,0.36,-0.16,-0.09,0.96,0.61,0.14,0.12
                    ,-0.95,-0.72,-0.68,-0.47,-0.88,-0.55,-0.60,-0.43,-0.63,-0.18,-0.43,-0.24,-0.19,0.22,-0.13,-0.03};
  for(unsigned int i = 0; i < 3; i++){
    for(unsigned int j = 0; j < 16; j++){
      wp[j][i] = temp[i*16+j];
    }
  }
  
  if(j_PileupJetId>=wp[binID][0]) return 0b111;
  else if(j_PileupJetId >= wp[binID][1]) return 0b11;
  else if(j_PileupJetId >= wp[binID][2]) return 0b1;
  else return 0b0;
}

