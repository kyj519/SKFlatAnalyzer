#include "Jet.h"

ClassImp(Jet)

    Jet::Jet() : Particle()
{
  j_area = -999.;
  j_partonFlavour = -999;
  j_hadronFlavour = -999;
  j_GenHFHadronMatcher_flavour = -999;
  j_GenHFHadronMatcher_origin = -999;
  j_DeepCSV = -999.;
  j_DeepCSV_CvsL = -999.;
  j_DeepCSV_CvsB = -999.;
  j_DeepJet = -999;
  j_DeepJet_CvsL = -999;
  j_DeepJet_CvsB = -999;
  j_chargedHadronEnergyFraction = -999.;
  j_neutralHadronEnergyFraction = -999.;
  j_neutralEmEnergyFraction = -999.;
  j_chargedEmEnergyFraction = -999.;
  j_muonEnergyFraction = -999.;
  j_chargedMultiplicity = -999;
  j_neutralMultiplicity = -999;
  j_PileupJetId = -999.;
  j_En_up = 1.;
  j_En_down = 1.;
  ;
  j_Res_up = 1.;
  j_Res_down = 1.;
  j_bJetNN_corr = 1.;
  j_bJetNN_res = -999.;
  j_cJetNN_corr = 1.;
  j_cJetNN_res = -999.;

  j_tightJetID = false;
  j_tightLepVetoJetID = false;
}

Jet::~Jet()
{
}

void Jet::SetArea(double area)
{
  j_area = area;
}
void Jet::SetGenFlavours(int pf, int hf)
{
  j_partonFlavour = pf;
  j_hadronFlavour = hf;
}
void Jet::SetGenHFHadronMatcher(int flavour, int origin)
{
  j_GenHFHadronMatcher_flavour = flavour;
  j_GenHFHadronMatcher_origin = origin;
}
void Jet::SetTaggerResults(std::vector<double> ds)
{
  j_DeepCSV = ds.at(0);
  j_DeepCSV_CvsL = ds.at(1);
  j_DeepCSV_CvsB = ds.at(2);
  j_DeepJet = ds.at(3);
  j_DeepJet_CvsL = ds.at(4);
  j_DeepJet_CvsB = ds.at(5);
}
void Jet::SetEnergyFractions(double cH, double nH, double nEM, double cEM, double muE)
{
  j_chargedHadronEnergyFraction = cH;
  j_neutralHadronEnergyFraction = nH;
  j_neutralEmEnergyFraction = nEM;
  j_chargedEmEnergyFraction = cEM;
  j_muonEnergyFraction = muE;
}
void Jet::SetMultiplicities(double cM, double nM)
{
  j_chargedMultiplicity = cM;
  j_neutralMultiplicity = nM;
}
void Jet::SetPileupJetId(double v)
{
  j_PileupJetId = v;
}

void Jet::SetEnShift(double en_up, double en_down)
{
  j_En_up = en_up;
  j_En_down = en_down;
}

void Jet::SetResShift(double res_up, double res_down)
{
  j_Res_up = res_up;
  j_Res_down = res_down;
}
void Jet::SetBJetNNCorrection(double bJetNN_corr, double bJetNN_res)
{
  j_bJetNN_corr = bJetNN_corr;
  j_bJetNN_res = bJetNN_res;
}
void Jet::SetCJetNNCorrection(double cJetNN_corr, double cJetNN_res)
{
  j_cJetNN_corr = cJetNN_corr;
  j_cJetNN_res = cJetNN_res;
}
void Jet::SetTightJetID(double b)
{
  j_tightJetID = b;
}
void Jet::SetTightLepVetoJetID(double b)
{
  j_tightLepVetoJetID = b;
}

//////////

bool Jet::Pass_PileupJetVeto(const TString &wp) const
{
  float pt = this->Pt();
  float eta = this->Eta();

  // PU Jet veto ID should be applied for 10<pt<50
  if (10. > pt || pt >= 50.)
  {
    return true;
  }

  unsigned int pt_bin = 0;
  if (10. <= pt && pt < 20.)
    pt_bin = 0;
  else if (20. <= pt && pt < 30.)
    pt_bin = 4;
  else if (30. <= pt && pt < 40.)
    pt_bin = 8;
  else if (40. <= pt && pt < 50.)
    pt_bin = 12;

  unsigned int eta_bin = 0;
  if (fabs(eta) < 2.5)
    eta_bin = 0;
  else if (2.5 <= fabs(eta) && fabs(eta) < 2.75)
    eta_bin = 1;
  else if (2.75 <= fabs(eta) && fabs(eta) < 3.0)
    eta_bin = 2;
  else if (3.0 <= fabs(eta) && fabs(eta) < 5.0)
    eta_bin = 3;

  unsigned int bin = pt_bin + eta_bin;

  double wp_cut[16][3] = {{0.77, 0.26, -0.95},
                          {0.38, -0.33, -0.72},
                          {-0.31, -0.54, -0.68},
                          {-0.21, -0.37, -0.47},
                          {0.90, 0.68, -0.88},
                          {0.60, -0.04, -0.55},
                          {-0.12, -0.43, -0.60},
                          {-0.13, -0.30, -0.43},
                          {0.96, 0.90, -0.63},
                          {0.82, 0.36, -0.18},
                          {0.20, -0.16, -0.43},
                          {0.09, -0.09, -0.24},
                          {0.98, 0.96, -0.19},
                          {0.92, 0.61, 0.22},
                          {0.47, 0.14, -0.13},
                          {0.29, 0.12, -0.03}};

  if (wp.Contains("Tight"))
  {
    if (wp_cut[bin][0] < j_PileupJetId)
      return true;
    else
      return false;
  }
  else if (wp.Contains("Medium"))
  {
    if (wp_cut[bin][1] < j_PileupJetId)
      return true;
    else
      return false;
  }
  else if (wp.Contains("Loose"))
  {
    if (wp_cut[bin][2] < j_PileupJetId)
      return true;
    else
      return false;
  }
  else
  {
    cout << "Wrong Jet::Pass_PileupJetVeto wrong wp_type: " << wp << endl;
    exit(ENODATA);
  }

  return false;
} // bool Jet::Pass_PileupJetVeto(const TString& wp) const

//////////

bool Jet::PassID(TString ID) const
{

  if (ID == "tight")
    return Pass_tightJetID();
  if (ID == "tightLepVeto")
    return Pass_tightLepVetoJetID();
  if (ID == "LoosePileupJetVeto")
    return Pass_PileupJetVeto("Loose");
  if (ID == "MediumPileupJetVeto")
    return Pass_PileupJetVeto("Medium");
  if (ID == "TightPileupJetVeto")
    return Pass_PileupJetVeto("Tight");

  cout << "[Jet::PassID] No id : " << ID << endl;
  exit(ENODATA);

  return false;
}

double Jet::GetTaggerResult(JetTagging::Tagger tg) const
{

  if (tg == JetTagging::DeepCSV)
    return j_DeepCSV;
  else if (tg == JetTagging::DeepCSV_CvsL)
    return j_DeepCSV_CvsL;
  else if (tg == JetTagging::DeepCSV_CvsB)
    return j_DeepCSV_CvsB;
  else if (tg == JetTagging::DeepJet)
    return j_DeepJet;
  else if (tg == JetTagging::DeepJet_CvsL)
    return j_DeepJet_CvsL;
  else if (tg == JetTagging::DeepJet_CvsB)
    return j_DeepJet_CvsB;
  else
  {
    cout << "[Jet::GetTaggerResult] ERROR; Wrong tagger : " << tg << endl;
    return -999;
  }
}
