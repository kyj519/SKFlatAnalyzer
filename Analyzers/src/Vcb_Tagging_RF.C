#include "Vcb_Tagging_RF.h"

//////////

Vcb_Tagging_RF::Vcb_Tagging_RF()
{
  result_tree = NULL;
} // Vcb_Tagging_RF::Vcb_Tagging_RF()

//////////

Vcb_Tagging_RF::~Vcb_Tagging_RF()
{
  outfile->cd();
  result_tree->Write();

} // Vcb_Tagging_RF::~Vcb_Tagging_RF()

//////////

void Vcb_Tagging_RF::initializeAnalyzer()
{
  run_mu_ch = HasFlag("RunMu");
  cout << "[Vcb_Tagging_RF::initializeAnalyzer] RunMu = " << run_mu_ch << endl;

  run_el_ch = HasFlag("RunEl");
  cout << "[Vcb_Tagging_RF::initializeAnalyzer] RunEl = " << run_el_ch << endl;

  if (run_mu_ch == run_el_ch)
  {
    cout << "One of RunMu or RunEl should be set." << endl;
    exit(1);
  }

  run_debug = HasFlag("RunDebug");
  cout << "[Vcb_Tagging_RF::initializeAnalyzer] RunDebug = " << run_debug << endl;

  // set single muon object
  vec_mu_id = {"POGTightWithTightIso"};
  vec_mu_id_sf_key = {"NUM_TightID_DEN_TrackerMuons"};
  vec_mu_iso_sf_key = {"NUM_TightRelIso_DEN_TightIDandIPCut"};

  // set single electron
  vec_el_id = {"passTightID"};
  vec_el_id_sf_key = {"ID_SF_passTightID"};

  // set sigle lepton trigger
  if (DataYear == 2016)
  {
    if (run_mu_ch)
    {
      vec_sl_trig.push_back("HLT_IsoMu24_v");
      sl_trig = "IsoMu24";
      sl_trig_safe_pt_cut = 26.;
    }
    else if (run_el_ch)
    {
      vec_sl_trig.push_back("HLT_Ele27_WPTight_Gsf_v");
      sl_trig = "Ele27";
      sl_trig_safe_pt_cut = 30.;
    }
  } // if (DataYear == 2016)
  else if (DataYear == 2017)
  {
    if (run_mu_ch)
    {
      vec_sl_trig.push_back("HLT_IsoMu27_v");
      sl_trig = "IsoMu27";
      sl_trig_safe_pt_cut = 30.;
    }
    else if (run_el_ch)
    {
      vec_sl_trig.push_back("HLT_Ele35_WPTight_Gsf_v");
      sl_trig = "Ele35";
      sl_trig_safe_pt_cut = 37.;
    }
  } // else if (DataYear == 2017)
  else if (DataYear == 2018)
  {
    if (run_mu_ch)
    {
      vec_sl_trig.push_back("HLT_IsoMu24_v");
      sl_trig = "IsoMu24";
      sl_trig_safe_pt_cut = 26.;
    }
    else if (run_el_ch)
    {
      vec_sl_trig.push_back("HLT_Ele32_WPTight_Gsf_v");
      sl_trig = "Ele32";
      sl_trig_safe_pt_cut = 35.;
    }
  } // else if (DataYear == 2018)
  else
    std::runtime_error("No trigger configuration for year");

  for (auto &trigger_name : vec_sl_trig)
    cout << "[Vcb::initializeAnalyzer] Single Lepton Trigger Name = " << trigger_name << endl;
  cout << "[Vcb::initializeAnalyzer] Single Lepton Trigger Safe Pt Cut = " << sl_trig_safe_pt_cut << endl;

  // Jet Tagging Parameters
  if (run_debug)
  {
    vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
    vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet_C, JetTagging::Medium, JetTagging::incl, JetTagging::wcharm));
  }
  else
  {
    vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::iterativefit, JetTagging::iterativefit));
    vec_jet_tagging_para.push_back(JetTagging::Parameters(JetTagging::DeepJet_C, JetTagging::Medium, JetTagging::iterativefit, JetTagging::iterativefit));
  }
  mcCorr->SetJetTaggingParameters(vec_jet_tagging_para);

  Set_Result_Tree();

  return;
} // void Vcb_Tagging_RF::initializeAnalyzer()

//////////

void Vcb_Tagging_RF::executeEvent()
{
  // init and clear
  vec_muon.clear();
  vec_electron.clear();
  vec_jet.clear();

  vec_muon = GetAllMuons();
  vec_electron = GetAllElectrons();
  vec_jet = GetAllJets();

  param.Clear();

  param.syst_ = AnalyzerParameter::Central;

  if (run_mu_ch)
  {
    param.Muon_Tight_ID = "POGTightWithTightIso";
    param.Muon_ID_SF_Key = "NUM_TightID_DEN_TrackerMuons";
    param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";

    param.Electron_Tight_ID = "passTightID";
    param.Electron_Loose_ID = "passLooseID";

    param.Name = param.Muon_Tight_ID + "_" + param.GetSystType();
  }
  else if (run_el_ch)
  {
    param.Electron_Tight_ID = "passTightID";
    param.Electron_Loose_ID = "passLooseID";
    param.Electron_ID_SF_Key = "ID_SF_passTightID";

    param.Muon_Tight_ID = "POGTightWithTightIso";

    param.Name = param.Electron_Tight_ID + "_" + param.GetSystType();
  }

  param.Jet_ID = "tight";
  param.PUJet_Veto_ID = "MediumPileupJetVeto";

  executeEventFromParameter(param);

  return;
} // void Vcb_Tagging_RF::executeEvent()

//////////

void Vcb_Tagging_RF::executeEventFromParameter(AnalyzerParameter param)
{
  Clear();

  Event ev = GetEvent();

  if (!IsData)
  {
    // lumi
    weight_lumi = ev.GetTriggerLumi("Full");
    weight *= weight_lumi;

    // MCweight +1 or -1
    weight_mc = MCweight();
    weight *= weight_mc;

    // pileup reweight
    weight_pileup = mcCorr->GetPileUpWeight(nPileUp, 0);
    weight *= weight_pileup;

    // L1 prefire
    weight_prefire = GetPrefireWeight(0);
    weight *= weight_prefire;

    // Top Pt reweight
    vec_gen = GetGens();
    weight_top_pt = mcCorr->GetTopPtReweight(vec_gen);
    weight *= weight_top_pt;
  }

  // met filter
  if (!PassMETFilter())
    return;

  // set objects
  vec_this_muon = vec_muon;
  vec_this_electron = vec_electron;
  vec_this_jet = vec_jet;

  met = ev.GetMETVector("PUPPI");
  pair<double, double> met_corr = xy_met_correction.METXYCorr_Met_MetPhi(met.Pt(), met.Phi(), run, to_string(DataYear), !IsData, nPV, true, true);

  // Particle met = ev.GetMETVector("PF");
  // pair<double, double> met_corr = xy_met_correction.METXYCorr_Met_MetPhi(met.Pt(), met.Phi(), run, to_string(DataYear), !IsData, nPV, true, false);
  // cout << met.Pt() << " " << met_corr.first << " " << met.Phi() << " " << met_corr.second << endl;

  met.SetPtEtaPhiE(met_corr.first, 0, met_corr.second, met_corr.first);

  /////////////////
  /*Setup Objects*/
  /////////////////

  // for
  vector<Muon> vec_sel_muon = SelectMuons(vec_this_muon, param.Muon_Tight_ID, sl_trig_safe_pt_cut, MUON_ETA);
  vector<Electron> vec_sel_electron = SelectElectrons(vec_this_electron, param.Electron_Tight_ID, sl_trig_safe_pt_cut, ELECTRON_ETA);

  // for lepton veto
  vector<Muon> vec_muon_veto = SelectMuons(vec_this_muon, param.Muon_Tight_ID, MUON_PT_VETO, MUON_ETA);
  vector<Electron> vec_electron_veto = SelectElectrons(vec_this_electron, param.Electron_Loose_ID, ELECTRON_PT_VETO, ELECTRON_ETA);

  float jet_eta_cut = 999;
  if (DataYear == 2016)
    jet_eta_cut = JET_ETA_2016;
  else if (DataYear == 2017 || DataYear == 2018)
    jet_eta_cut = JET_ETA;

  // Jet selection
  vec_sel_jet = SelectJets(vec_this_jet, param.Jet_ID, JET_PT, jet_eta_cut);
  vec_sel_jet = JetsVetoLeptonInside(vec_sel_jet, vec_electron_veto, vec_muon_veto, DR_LEPTON_VETO);
  n_sel_jet = vec_sel_jet.size();

  // sort jet as pt ordering
  sort(vec_sel_jet.begin(), vec_sel_jet.end(), PtComparing);

  // single lepton trigger
  if (!ev.PassTrigger(vec_sl_trig))
    return;

  if (!IsData)
  {
    // SF for muon trigger effi
    if (run_mu_ch)
      sf_sl_trig = mcCorr->MuonTrigger_SF("POGTight", sl_trig, vec_sel_muon, 0);

    // SF for electron trigger effi
    else if (run_el_ch)
      sf_sl_trig = mcCorr->ElectronTrigger_SF(param.Electron_Tight_ID, sl_trig, vec_sel_electron, 0);

    weight *= sf_sl_trig;
  }

  // cut on lepton
  // veto additional lepton
  if (run_mu_ch)
  {
    if (vec_sel_muon.size() != 1)
      return;
    if (vec_muon_veto.size() != 1)
      return;
    if (vec_electron_veto.size() != 0)
      return;

    muon = vec_sel_muon.at(0);
    lepton = muon;
  }
  else if (run_el_ch)
  {
    if (vec_sel_electron.size() != 1)
      return;
    if (vec_electron_veto.size() != 1)
      return;
    if (vec_muon_veto.size() != 0)
      return;

    electron = vec_sel_electron.at(0);
    lepton = electron;
  }

  lepton_pt = lepton.Pt();
  lepton_eta = lepton.Eta();

  if (lepton_pt <= sl_trig_safe_pt_cut)
    return;

  if (!IsData)
  {
    if (run_mu_ch)
    {
      // SF for muon id
      sf_mu_id = mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muon.Eta(), muon.MiniAODPt(), 0);
      weight *= sf_mu_id;

      // SF for muon iso
      sf_mu_iso = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muon.Eta(), muon.MiniAODPt(), 0);
      weight *= sf_mu_iso;

      // if(run_debug) cout << "M ISO SF = " << muon.Eta() << "\t" << muon.MiniAODPt() << "\t" << sf_mu_iso_effi << endl;
    }
    else if (run_el_ch)
    {
      // SF for electron id
      sf_el_id = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electron.scEta(), electron.UncorrPt(), 0);
      weight *= sf_el_id;

      // SF for electron Reco eff
      sf_el_reco = mcCorr->ElectronReco_SF(electron.scEta(), electron.UncorrPt(), 0);
      weight *= sf_el_reco;
    }
  }

  // cut on jet
  if (n_sel_jet < 4)
    return;

  if (!IsData)
  {
    // SF for PUJet Veto
    weight_pujet_veto = mcCorr->PileupJetVeto_Reweight(vec_sel_jet, param.PUJet_Veto_ID, 0);
    // weight *= weight_pujet_veto;
  }

  // n of btag
  n_b_jet = 0;
  for (auto &jet : vec_sel_jet)
  {
    float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
    if (mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score)
    {
      n_b_jet++;
      vec_btag.push_back(true);
    }
    else
      vec_btag.push_back(false);
  }

  // n of ctag
  n_c_jet = 0;
  for (auto &jet : vec_sel_jet)
  {
    float cvsb = jet.GetTaggerResult(JetTagging::DeepJet_CvsB);
    float cvsl = jet.GetTaggerResult(JetTagging::DeepJet_CvsL);

    float cvsl_wp = -1;
    float cvsb_wp = -1;

    if (DataYear == 2017)
    {
      cvsl_wp = cvsl_2017_m;
      cvsb_wp = cvsb_2017_m;
    }
    else if (DataYear == 2018)
    {
      cvsl_wp = cvsl_2018_m;
      cvsb_wp = cvsb_2018_m;
    }

    if (cvsl_wp < cvsl && cvsb_wp < cvsb)
    {
      n_c_jet++;
      vec_ctag.push_back(true);
    }
    else
      vec_ctag.push_back(false);
  }
  // cout << "test n_of_ctag" << n_c_jet << endl;

  if (!IsData)
  {
    // SF for b-tagging
    weight_b_tag = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "central");

    weight_b_tag_down_hf = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_hf");
    weight_b_tag_up_hf = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_hf");

    weight_b_tag_down_jes = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_jes");
    weight_b_tag_up_jes = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_jes");

    weight_b_tag_down_lfstats1 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_lfstats1");
    weight_b_tag_up_lfstats1 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_lfstats1");

    weight_b_tag_down_lfstats2 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_lfstats2");
    weight_b_tag_up_lfstats2 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_lfstats2");

    weight_b_tag_down_cferr1 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_cferr1");
    weight_b_tag_up_cferr1 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_cferr1");

    weight_b_tag_down_cferr2 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_cferr2");
    weight_b_tag_up_cferr2 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_cferr2");

    weight_b_tag_down_hfstats1 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_hfstats1");
    weight_b_tag_up_hfstats1 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_hfstats1");

    weight_b_tag_down_hfstats2 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_hfstats2");
    weight_b_tag_up_hfstats2 = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_hfstats2");

    weight *= weight_b_tag;
  }

  if (!IsData)
  {
    // SF for c-tagging
    weight_c_tag = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "central");

    weight_c_tag_down_extrap = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "Extrap_Down");
    weight_c_tag_up_extrap = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "Extrap_Up");

    weight_c_tag_down_interp = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "Interp_Down");
    weight_c_tag_up_interp = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "Interp_Up");

    weight_c_tag_down_lhe_scale_muf = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "LHEScaleWeight_muF_Down");
    weight_c_tag_up_lhe_scale_muf = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "LHEScaleWeight_muF_Up");

    weight_c_tag_down_lhe_scale_mur = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "LHEScaleWeight_muR_Down");
    weight_c_tag_up_lhe_scale_mur = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "LHEScaleWeight_muR_Up");

    weight_c_tag_down_ps_fsr_fixed = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "PSWeightFSRFixed_Down");
    weight_c_tag_up_ps_fsr_fixed = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "PSWeightFSRFixed_Up");

    weight_c_tag_down_ps_isr_fixed = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "PSWeightISRFixed_Down");
    weight_c_tag_up_ps_isr_fixed = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "PSWeightISRFixed_Up");

    weight_c_tag_down_pu = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "PUWeight_Down");
    weight_c_tag_up_pu = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "PUWeight_Up");

    weight_c_tag_down_stat = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "Stat_Down");
    weight_c_tag_up_stat = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "Stat_Up");

    weight_c_tag_down_xsec_brunc_dyjets_b = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "XSec_BRUnc_DYJets_b_Down");
    weight_c_tag_up_xsec_brunc_dyjets_b = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "XSec_BRUnc_DYJets_b_Up");

    weight_c_tag_down_xsec_brunc_dyjets_c = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "XSec_BRUnc_DYJets_c_Down");
    weight_c_tag_up_xsec_brunc_dyjets_c = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "XSec_BRUnc_DYJets_c_Up");

    weight_c_tag_down_xsec_brunc_wjets_c = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "XSec_BRUnc_WJets_c_Down");
    weight_c_tag_up_xsec_brunc_wjets_c = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "XSec_BRUnc_WJets_c_Up");

    weight_c_tag_down_jer = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jer_Down");
    weight_c_tag_up_jer = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jer_Up");

    weight_c_tag_down_jes_total = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jesTotal_Down");
    weight_c_tag_up_jes_total = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jesTotal_Up");

    weight *= weight_c_tag;
  }

  // no b-tag is required
  // if (n_b_jet < 2)
  //  return;

  // cut on MET
  met_pt = met.Pt();
  // met_phi = met.Phi();

  if (met_pt < MET_PT)
    return;

  leading_jet_bvsc = vec_sel_jet[0].GetTaggerResult(JetTagging::DeepJet);
  leading_jet_cvsb = vec_sel_jet[0].GetTaggerResult(JetTagging::DeepJet_CvsB);
  leading_jet_cvsl = vec_sel_jet[0].GetTaggerResult(JetTagging::DeepJet_CvsL);
  leading_jet_eta = vec_sel_jet[0].Eta();
  leading_jet_pt = vec_sel_jet[0].Pt();

  subleading_jet_bvsc = vec_sel_jet[1].GetTaggerResult(JetTagging::DeepJet);
  subleading_jet_cvsb = vec_sel_jet[1].GetTaggerResult(JetTagging::DeepJet_CvsB);
  subleading_jet_cvsl = vec_sel_jet[1].GetTaggerResult(JetTagging::DeepJet_CvsL);
  subleading_jet_eta = vec_sel_jet[1].Eta();
  subleading_jet_pt = vec_sel_jet[1].Pt();

  ht = Calculate_HT(vec_sel_jet);

  result_tree->Fill();

  return;
} // Vcb_Tagging_RF::executeEventFromParameter(AnalyzerParameter param)

//////////

float Vcb_Tagging_RF::Calculate_HT(const vector<Jet> &jet)
{
  float ht = 0;
  for (unsigned int i = 0; i < jet.size(); i++)
  {
    ht += jet[i].Pt();
  }

  return ht;
} // float Vcb_Tagging_RF::Calculate_HT(const vector<Jet>& jet)

//////////

void Vcb_Tagging_RF::Clear()
{
  weight = 1;

  weight_b_tag = 1;
  weight_c_tag = 1;
  weight_lumi = 1;
  weight_mc = 1;
  weight_pileup = 1;
  weight_pujet_veto = 1;
  weight_prefire = 1;
  weight_top_pt = 1;

  sf_mu_id = 1;
  sf_mu_iso = 1;
  sf_mu_trig = 1;

  sf_el_id = 1;
  sf_el_reco = 1;

  sf_sl_trig = 1;

  vec_this_muon.clear();
  vec_this_electron.clear();
  vec_this_jet.clear();

  vec_sel_jet.clear();

  vec_btag.clear();
  vec_ctag.clear();

  return;
} // void Vcb_Tagging_RF::Clear()

//////////

void Vcb_Tagging_RF::Set_Result_Tree()
{
  result_tree = new TTree("Result_Tree", "Result_Tree");

  result_tree->Branch("sf_mu_id", &sf_mu_id);
  result_tree->Branch("sf_mu_iso", &sf_mu_iso);
  result_tree->Branch("sf_mu_trig", &sf_mu_trig);

  result_tree->Branch("sf_el_id", &sf_el_id);
  result_tree->Branch("sf_el_reco", &sf_el_reco);

  result_tree->Branch("sf_sl_trig", &sf_sl_trig);

  result_tree->Branch("weight_b_tag", &weight_b_tag);
  result_tree->Branch("weight_b_tag_down_hf", &weight_b_tag_down_hf);
  result_tree->Branch("weight_b_tag_up_hf", &weight_b_tag_up_hf);
  result_tree->Branch("weight_b_tag_down_jes", &weight_b_tag_down_jes);
  result_tree->Branch("weight_b_tag_up_jes", &weight_b_tag_up_jes);
  result_tree->Branch("weight_b_tag_down_lfstats1", &weight_b_tag_down_lfstats1);
  result_tree->Branch("weight_b_tag_up_lfstats1", &weight_b_tag_up_lfstats1);
  result_tree->Branch("weight_b_tag_down_lfstats2", &weight_b_tag_down_lfstats2);
  result_tree->Branch("weight_b_tag_up_lfstats2", &weight_b_tag_up_lfstats2);
  result_tree->Branch("weight_b_tag_down_cferr1", &weight_b_tag_down_cferr1);
  result_tree->Branch("weight_b_tag_up_cferr1", &weight_b_tag_up_cferr1);
  result_tree->Branch("weight_b_tag_down_cferr2", &weight_b_tag_down_cferr2);
  result_tree->Branch("weight_b_tag_up_cferr2", &weight_b_tag_up_cferr2);
  result_tree->Branch("weight_b_tag_down_hfstats1", &weight_b_tag_down_hfstats1);
  result_tree->Branch("weight_b_tag_up_hfstats1", &weight_b_tag_up_hfstats1);
  result_tree->Branch("weight_b_tag_down_hfstats2", &weight_b_tag_down_hfstats2);
  result_tree->Branch("weight_b_tag_up_hfstats2", &weight_b_tag_up_hfstats2);

  result_tree->Branch("weight_c_tag", &weight_c_tag);
  result_tree->Branch("weight_c_tag_down_extrap", &weight_c_tag_down_extrap);
  result_tree->Branch("weight_c_tag_up_extrap", &weight_c_tag_up_extrap);
  result_tree->Branch("weight_c_tag_down_interp", &weight_c_tag_down_interp);
  result_tree->Branch("weight_c_tag_up_interp", &weight_c_tag_up_interp);
  result_tree->Branch("weight_c_tag_down_lhe_scale_muf", &weight_c_tag_down_lhe_scale_muf);
  result_tree->Branch("weight_c_tag_up_lhe_scale_muf", &weight_c_tag_up_lhe_scale_muf);
  result_tree->Branch("weight_c_tag_down_lhe_scale_mur", &weight_c_tag_down_lhe_scale_mur);
  result_tree->Branch("weight_c_tag_up_lhe_scale_mur", &weight_c_tag_up_lhe_scale_mur);
  result_tree->Branch("weight_c_tag_down_ps_fsr_fixed", &weight_c_tag_down_ps_fsr_fixed);
  result_tree->Branch("weight_c_tag_up_ps_fsr_fixed", &weight_c_tag_up_ps_fsr_fixed);
  result_tree->Branch("weight_c_tag_down_ps_isr_fixed", &weight_c_tag_down_ps_isr_fixed);
  result_tree->Branch("weight_c_tag_up_ps_isr_fixed", &weight_c_tag_up_ps_isr_fixed);
  result_tree->Branch("weight_c_tag_down_pu", &weight_c_tag_down_pu);
  result_tree->Branch("weight_c_tag_up_pu", &weight_c_tag_up_pu);
  result_tree->Branch("weight_c_tag_down_stat", &weight_c_tag_down_stat);
  result_tree->Branch("weight_c_tag_up_stat", &weight_c_tag_up_stat);
  result_tree->Branch("weight_c_tag_down_xsec_brunc_dyjets_b", &weight_c_tag_down_xsec_brunc_dyjets_b);
  result_tree->Branch("weight_c_tag_up_xsec_brunc_dyjets_b", &weight_c_tag_up_xsec_brunc_dyjets_b);
  result_tree->Branch("weight_c_tag_down_xsec_brunc_dyjets_c", &weight_c_tag_down_xsec_brunc_dyjets_c);
  result_tree->Branch("weight_c_tag_up_xsec_brunc_dyjets_c", &weight_c_tag_up_xsec_brunc_dyjets_c);
  result_tree->Branch("weight_c_tag_down_xsec_brunc_wjets_c", &weight_c_tag_down_xsec_brunc_wjets_c);
  result_tree->Branch("weight_c_tag_up_xsec_brunc_wjets_c", &weight_c_tag_up_xsec_brunc_wjets_c);
  result_tree->Branch("weight_c_tag_down_jer", &weight_c_tag_down_jer);
  result_tree->Branch("weight_c_tag_up_jer", &weight_c_tag_up_jer);
  result_tree->Branch("weight_c_tag_down_jes_total", &weight_c_tag_down_jes_total);
  result_tree->Branch("weight_c_tag_up_jes_total", &weight_c_tag_up_jes_total);

  result_tree->Branch("weight_lumi", &weight_lumi);
  result_tree->Branch("weight_mc", &weight_mc);
  result_tree->Branch("weight_pileup", &weight_pileup);
  result_tree->Branch("weight_prefire", &weight_prefire);
  result_tree->Branch("weight_pujet_veto", &weight_pujet_veto);
  result_tree->Branch("weight_top_pt", &weight_top_pt);

  result_tree->Branch("n_vertex", &nPV);

  // result_tree->Branch("lepton_pt", &lepton_pt);
  // result_tree->Branch("lepton_eta", &lepton_eta);

  result_tree->Branch("n_jets", &n_sel_jet);
  result_tree->Branch("n_bjets", &n_b_jet);
  result_tree->Branch("n_cjets", &n_c_jet);

  result_tree->Branch("ht", &ht);

  result_tree->Branch("leading_jet_bvsc", &leading_jet_bvsc);
  result_tree->Branch("leading_jet_cvsb", &leading_jet_cvsb);
  result_tree->Branch("leading_jet_cvsl", &leading_jet_cvsl);
  result_tree->Branch("leading_jet_eta", &leading_jet_eta);
  result_tree->Branch("leading_jet_pt", &leading_jet_pt);

  result_tree->Branch("subleading_jet_bvsc", &subleading_jet_bvsc);
  result_tree->Branch("subleading_jet_cvsb", &subleading_jet_cvsb);
  result_tree->Branch("subleading_jet_cvsl", &subleading_jet_cvsl);
  result_tree->Branch("subleading_jet_eta", &subleading_jet_eta);
  result_tree->Branch("subleading_jet_pt", &subleading_jet_pt);

  return;
} // void Vcb_Tagging_RF::Set_Result_Tree()

//////////