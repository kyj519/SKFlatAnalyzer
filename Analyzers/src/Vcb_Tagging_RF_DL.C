#include "Vcb_Tagging_RF_DL.h"

//////////

Vcb_Tagging_RF_DL::Vcb_Tagging_RF_DL()
{
} // Vcb_Tagging_RF_DL::Vcb_Tagging_RF_DL()

//////////

Vcb_Tagging_RF_DL::~Vcb_Tagging_RF_DL()
{
  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    param.syst_ = vec_syst_type.at(i);

    outfile->cd(param.GetSystType());
    map_result_tree[param.syst_]->Write();
  }
} // Vcb_Tagging_RF_DL::~Vcb_Tagging_RF_DL()

//////////

void Vcb_Tagging_RF_DL::initializeAnalyzer()
{
  run_mm_ch = HasFlag("RunMM");
  cout << "[Vcb_Tagging_RF_DL::initializeAnalyzer] RunMM = " << run_mm_ch << endl;

  run_ee_ch = HasFlag("RunEE");
  cout << "[Vcb_Tagging_RF_DL::initializeAnalyzer] RunEE = " << run_ee_ch << endl;

  run_me_ch = HasFlag("RunME");
  cout << "[Vcb_Tagging_RF_DL::initializeAnalyzer] RunME = " << run_me_ch << endl;

  int chk_ch = (int)run_mm_ch + (int)run_ee_ch + (int)run_me_ch;
  if (chk_ch != 1)
  {
    cout << "One of RunMM, RunEE or RunME should be set." << endl;
    exit(1);
  }

  run_debug = HasFlag("RunDebug");
  cout << "[Vcb_Tagging_RF_DL::initializeAnalyzer] RunDebug = " << run_debug << endl;

  // set single muon object
  // vec_mu_id = {"POGTightWithTightIso"};
  // vec_mu_id_sf_key = {"NUM_TightID_DEN_TrackerMuons"};
  // vec_mu_iso_sf_key = {"NUM_TightRelIso_DEN_TightIDandIPCut"};

  // set single electron
  // vec_el_id = {"passTightID"};
  // vec_el_id_sf_key = {"ID_SF_passTightID"};

  // set sigle lepton trigger
  if (DataYear == 2016)
  {
    vec_mu_trig.push_back("HLT_IsoMu24_v");
    mu_trig = "IsoMu24";
    mu_trig_safe_pt_cut = 26.;

    vec_el_trig.push_back("HLT_Ele27_WPTight_Gsf_v");
    el_trig = "Ele27";
    el_trig_safe_pt_cut = 30.;
  } // if (DataYear == 2016)
  else if (DataYear == 2017)
  {
    vec_mu_trig.push_back("HLT_IsoMu27_v");
    mu_trig = "IsoMu27";
    mu_trig_safe_pt_cut = 30.;

    // vec_el_trig.push_back("HLT_Ele35_WPTight_Gsf_v");
    // el_trig = "Ele35";
    // el_trig_safe_pt_cut = 37.;

    vec_el_trig.push_back("HLT_Ele32_WPTight_Gsf_L1DoubleEG_");
    el_trig = "Ele32";
    el_trig_safe_pt_cut = 35.;
  } // else if (DataYear == 2017)
  else if (DataYear == 2018)
  {
    vec_mu_trig.push_back("HLT_IsoMu24_v");
    mu_trig = "IsoMu24";
    mu_trig_safe_pt_cut = 26.;

    vec_el_trig.push_back("HLT_Ele32_WPTight_Gsf_v");
    el_trig = "Ele32";
    el_trig_safe_pt_cut = 35.;
  } // else if (DataYear == 2018)
  else
    std::runtime_error("No trigger configuration for year");

  for (auto &trigger_name : vec_mu_trig)
    vec_sl_trig.push_back(trigger_name);
  for (auto &trigger_name : vec_el_trig)
    vec_sl_trig.push_back(trigger_name);

  for (auto &trigger_name : vec_sl_trig)
    cout << "[Vcb::initializeAnalyzer] Single Lepton Trigger Name = " << trigger_name << endl;
  cout << "[Vcb::initializeAnalyzer] Single Muon Trigger Safe Pt Cut = " << mu_trig_safe_pt_cut << endl;
  cout << "[Vcb::initializeAnalyzer] Single Electron Trigger Safe Pt Cut = " << el_trig_safe_pt_cut << endl;

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

  if (MCSample.Contains("mtop17") || MCSample.Contains("CP5") || MCSample.Contains("hdamp"))
    vec_syst_type = {AnalyzerParameter::Central};
  else
    vec_syst_type = {AnalyzerParameter::Central,
                     AnalyzerParameter::JetEnDown,
                     AnalyzerParameter::JetEnUp,
                     AnalyzerParameter::JetResDown,
                     AnalyzerParameter::JetResUp};

  // to make output dir
  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    param.syst_ = vec_syst_type.at(i);
    FillHist(param.GetSystType() + "/Dummy", 0, weight, 1, 0, 1);
  }

  Set_Result_Tree();

  return;
} // void Vcb_Tagging_RF_DL::initializeAnalyzer()

//////////

void Vcb_Tagging_RF_DL::executeEvent()
{
  // init and clear
  vec_muon.clear();
  vec_electron.clear();
  vec_jet.clear();

  vec_muon = GetAllMuons();
  vec_electron = GetAllElectrons();
  vec_jet = GetAllJets();

  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    param.Clear();

    param.Muon_Tight_ID = "POGTightWithTightIso";
    param.Muon_Loose_ID = "POGLoose";
    param.Muon_ID_SF_Key = "NUM_TightID_DEN_TrackerMuons";
    param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";

    // param.Electron_Tight_ID = "passTightID";
    // param.Electron_Loose_ID = "passLooseID";
    // param.Electron_ID_SF_Key = "ID_SF_passTightID";

    param.Electron_Tight_ID = "passMVAID_iso_WP80";
    param.Electron_Loose_ID = "passMVAID_iso_WP90";
    // param.Electron_ID_SF_Key = "ID_SF_passTightID";

    param.Jet_ID = "tight";
    param.PUJet_Veto_ID = "LoosePileupJetVeto";

    param.syst_ = vec_syst_type.at(i);

    param.Name = param.GetSystType();

    vec_gen_hf_flavour.clear();
    vec_gen_hf_origin.clear();

    vec_sel_gen_hf_flavour.clear();
    vec_sel_gen_hf_origin.clear();

    executeEventFromParameter(param);
  }

  return;
} // void Vcb_Tagging_RF_DL::executeEvent()

//////////

void Vcb_Tagging_RF_DL::executeEventFromParameter(AnalyzerParameter param)
{
  Clear();

  Event ev = GetEvent();

  if (!IsData)
  {
    vec_gen = GetGens();
    decay_mode = Get_W_Decay_Mode(vec_gen);

    for (unsigned int i = 0; i < vec_jet.size(); i++)
    {
      Jet jet = vec_jet[i];

      // acceptance cut
      if (jet.Pt() < JET_PT)
        continue;

      float jet_eta_cut = JET_ETA;
      if (DataYear == 2016)
        jet_eta_cut = JET_ETA_2016;

      if (jet_eta_cut < abs(jet.Eta()))
        continue;

      vec_gen_hf_flavour.push_back(jet.GenHFHadronMatcherFlavour());
      vec_gen_hf_origin.push_back(jet.GenHFHadronMatcherOrigin());
    }

    // lumi
    weight_lumi = ev.GetTriggerLumi("Full");
    weight *= weight_lumi;

    // MCweight +1 or -1
    weight_mc = MCweight();
    weight *= weight_mc;

    // pileup reweight
    weight_pileup = mcCorr->GetPileUpWeight(nPileUp, 0);
    if (param.syst_ == AnalyzerParameter::Central)
    {
      weight_pileup_down = mcCorr->GetPileUpWeight(nPileUp, -1);
      weight_pileup_up = mcCorr->GetPileUpWeight(nPileUp, +1);
    }
    weight *= weight_pileup;

    // L1 prefire
    weight_prefire = GetPrefireWeight(0);
    weight *= weight_prefire;

    // Top Pt reweight
    weight_top_pt = mcCorr->GetTopPtReweight(vec_gen);
    weight *= weight_top_pt;

    // Scale Variation
    if (param.syst_ == AnalyzerParameter::Central)
    {
      weight_scale_variation_1 = GetScaleVariation(1);
      weight_scale_variation_2 = GetScaleVariation(2);
      weight_scale_variation_3 = GetScaleVariation(3);
      weight_scale_variation_4 = GetScaleVariation(4);
      weight_scale_variation_6 = GetScaleVariation(6);
      weight_scale_variation_8 = GetScaleVariation(8);
    }

    // PS Reweight
    Get_Reweight_PS(weight_ps);
  }

  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::No_Cut, weight, n_cut_flow, 0, n_cut_flow);

  // met filter
  if (!PassMETFilter())
    return;

  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Met_Filter, weight, n_cut_flow, 0, n_cut_flow);

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

  //////////////////////
  /* syst for objects */
  //////////////////////

  if (param.syst_ == AnalyzerParameter::JetEnDown)
  {
    vec_this_jet = ScaleJets(vec_jet, -1);
    met = Rebalance_Met();
  }
  else if (param.syst_ == AnalyzerParameter::JetEnUp)
  {
    vec_this_jet = ScaleJets(vec_jet, +1);
    met = Rebalance_Met();
  }

  if (param.syst_ == AnalyzerParameter::JetResDown)
  {
    vec_this_jet = SmearJets(vec_jet, -1);
    met = Rebalance_Met();
  }
  else if (param.syst_ == AnalyzerParameter::JetResUp)
  {
    vec_this_jet = SmearJets(vec_jet, +1);
    met = Rebalance_Met();
  }

  /////////////////
  /*Setup Objects*/
  /////////////////

  // for lepton
  vector<Muon> vec_sel_muon = SelectMuons(vec_this_muon, param.Muon_Tight_ID, sl_trig_safe_pt_cut, MUON_ETA);
  vector<Electron> vec_sel_electron = SelectElectrons(vec_this_electron, param.Electron_Tight_ID, sl_trig_safe_pt_cut, ELECTRON_ETA);

  // for lepton veto
  vector<Muon> vec_muon_veto = SelectMuons(vec_this_muon, param.Muon_Loose_ID, MUON_PT_VETO, MUON_ETA);
  vector<Electron> vec_electron_veto = SelectElectrons(vec_this_electron, param.Electron_Loose_ID, ELECTRON_PT_VETO, ELECTRON_ETA);

  float jet_eta_cut = 999;
  if (DataYear == 2016)
    jet_eta_cut = JET_ETA_2016;
  else if (DataYear == 2017 || DataYear == 2018)
    jet_eta_cut = JET_ETA;

  // Jet selection
  vec_sel_jet = SelectJets(vec_this_jet, param.Jet_ID, JET_PT, jet_eta_cut);
  vec_sel_jet = SelectJets(vec_sel_jet, param.PUJet_Veto_ID, JET_PT, jet_eta_cut);
  vec_sel_jet = JetsVetoLeptonInside(vec_sel_jet, vec_electron_veto, vec_muon_veto, DR_LEPTON_VETO);
  n_sel_jet = vec_sel_jet.size();

  // sort jet as pt ordering
  sort(vec_sel_jet.begin(), vec_sel_jet.end(), PtComparing);

  if (!IsData)
  {
    // SF for PUJet Veto
    weight_pujet_veto = mcCorr->PileupJetVeto_Reweight(vec_sel_jet, param.PUJet_Veto_ID, 0);
    weight *= weight_pujet_veto;
  }

  for (unsigned int i = 0; i < vec_sel_jet.size(); i++)
  {
    Jet jet = vec_sel_jet[i];

    vec_sel_gen_hf_flavour.push_back(jet.GenHFHadronMatcherFlavour());
    vec_sel_gen_hf_origin.push_back(jet.GenHFHadronMatcherOrigin());
  }

  weight_hem_veto = Weight_HEM_Veto(vec_sel_jet);
  weight *= weight_hem_veto;

  // single lepton trigger
  if (!ev.PassTrigger(vec_mu_trig) && (!ev.PassTrigger(vec_el_trig) || !HLT_SE_Filter_2017(vec_sel_electron)))
    return;

  if (!IsData)
  {
    weight_sl_trig = mcCorr->SingleLepton_Trigger_SF("POGTight", mu_trig, vec_sel_muon, 0, "passTightID", el_trig, vec_sel_electron, 0);
    weight *= weight_sl_trig;

    // muon reweight
    for (unsigned int i = 0; i < vec_sel_muon.size(); i++)
    {
      weight_mu_id *= mcCorr->MuonID_SF(param.Muon_ID_SF_Key, vec_sel_muon[i].Eta(), vec_sel_muon[i].MiniAODPt(), 0);
      weight_mu_iso *= mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, vec_sel_muon[i].Eta(), vec_sel_muon[i].MiniAODPt(), 0);
    }
    weight *= weight_mu_id;
    weight *= weight_mu_iso;

    // electron reweight
    for (unsigned int i = 0; i < vec_sel_electron.size(); i++)
    {
      weight_el_id *= mcCorr->ElectronID_SF(param.Electron_Tight_ID, vec_sel_electron[i].scEta(), vec_sel_electron[i].UncorrPt(), 0);
      weight_el_reco *= mcCorr->ElectronReco_SF(vec_sel_electron[i].scEta(), vec_sel_electron[i].UncorrPt(), 0);
    }
    weight *= weight_el_id;
    weight *= weight_el_reco;
  } // if (!IsData)

  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Trigger, weight, n_cut_flow, 0, n_cut_flow);

  // cut on double lepton and veto additional loose electron
  if (run_mm_ch)
  {
    if (vec_sel_muon.size() != 2 || vec_electron_veto.size() != 0)
    //if (vec_sel_muon.size() != 2)
      return;

    lepton[0] = vec_sel_muon[0];
    lepton[1] = vec_sel_muon[1];
  }
  else if (run_me_ch)
  {
    if (vec_sel_muon.size() != 1 || vec_sel_electron.size() != 1 || vec_electron_veto.size() != 1)
      return;

    lepton[0] = vec_sel_muon[0];
    lepton[1] = vec_sel_electron[0];
  }
  else if (run_ee_ch)
  {
    if (vec_sel_muon.size() != 0 || vec_sel_electron.size() != 2 || vec_electron_veto.size() != 2)
      return;

    lepton[0] = vec_sel_electron[0];
    lepton[1] = vec_sel_electron[1];
  }

  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Two_Lepton, weight, n_cut_flow, 0, n_cut_flow);

  // lepton pt sorting
  if (lepton[0].Pt() < lepton[1].Pt())
  {
    Lepton temp = lepton[1];
    lepton[1] = lepton[0];
    lepton[0] = temp;
  }

  // cut on lepton charge
  if (lepton[0].Charge() + lepton[1].Charge() != 0)
    return;

  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Lepton_Charge_Sum, weight, n_cut_flow, 0, n_cut_flow);

  // cut on dilepton mass to veto low mass resonance
  TLorentzVector dilepton = lepton[0];
  dilepton += lepton[1];
  dilepton_mass = dilepton.M();
  if (dilepton_mass < 15)
    return;

  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Lepton_Low_Mass, weight, n_cut_flow, 0, n_cut_flow);

  // cut on dilepton mass to veto Z
  if (abs(dilepton_mass - Z_MASS) < 15)
    return;

  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Z_Mass, weight, n_cut_flow, 0, n_cut_flow);

  // cut on jet
  if (n_sel_jet < 4)
    return;

  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::N_Sel_Jet, weight, n_cut_flow, 0, n_cut_flow);

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
    if (param.syst_ == AnalyzerParameter::Central)
    {
      weight_b_tag = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "central");

      weight_b_tag_down_hf = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_hf");
      weight_b_tag_up_hf = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_hf");

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
    }
    else if (param.syst_ == AnalyzerParameter::JetEnDown)
      weight_b_tag_down_jes = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_jes");
    else if (param.syst_ == AnalyzerParameter::JetEnUp)
      weight_b_tag_up_jes = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_jes");
  }

  if (!IsData)
  {
    // SF for c-tagging
    if (param.syst_ == AnalyzerParameter::Central)
    {
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
    } // central
    else if (param.syst_ == AnalyzerParameter::JetEnDown)
      weight_c_tag_down_jes_total = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jesTotal_Down");
    else if (param.syst_ == AnalyzerParameter::JetEnUp)
      weight_c_tag_up_jes_total = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jesTotal_Up");
    else if (param.syst_ == AnalyzerParameter::JetResDown)
      weight_c_tag_down_jer = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jer_Down");
    else if (param.syst_ == AnalyzerParameter::JetResUp)
      weight_c_tag_up_jer = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jer_Up");
  } // if (!IsData)

  // no b-tag is required
  // if (n_b_jet < 2)
  //  return;

  // cut on MET
  met_pt = met.Pt();
  // met_phi = met.Phi();

  if (!run_me_ch && met_pt < MET_PT_DL)
    return;

  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Met, weight, n_cut_flow, 0, n_cut_flow);

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

  map_result_tree[param.syst_]->Fill();

  return;
} // Vcb_Tagging_RF_DL::executeEventFromParameter(AnalyzerParameter param)

//////////

float Vcb_Tagging_RF_DL::Calculate_HT(const vector<Jet> &jet)
{
  float ht = 0;
  for (unsigned int i = 0; i < jet.size(); i++)
  {
    ht += jet[i].Pt();
  }

  return ht;
} // float Vcb_Tagging_RF_DL::Calculate_HT(const vector<Jet>& jet)

//////////

int Vcb_Tagging_RF_DL::Check_Process(const vector<Gen> &vec_gen)
{
  // PrintGen(vec_gen);

  int index_last_t = -999;
  int index_last_at = -999;
  // int index_first_b_from_t = -999;
  // int index_first_ab_from_at = -999;
  int index_first_b_from_g = -999;
  int index_first_ab_from_g = -999;
  int index_first_c_from_g = -999;
  int index_first_ac_from_g = -999;
  for (unsigned int i = 0; i < vec_gen.size(); i++)
  {
    Gen gen = vec_gen.at(i);

    int pid = gen.PID();
    int m_index = gen.MotherIndex();
    int m_pid = vec_gen[m_index].PID();

    // find last index of t and tbar
    if (pid == 6)
      index_last_t = i;
    if (pid == -6)
      index_last_at = i;

    // // find b from t decay
    // if (m_index == index_last_t && pid == 5)
    //   index_first_b_from_t = i;
    // if (m_index == index_last_at && pid == -5)
    //   index_first_ab_from_at = i;

    // find b from g radiation
    if (m_pid == 21 && pid == 5)
      index_first_b_from_g = i;
    if (m_pid == 21 && pid == -5)
      index_first_ab_from_g = i;

    // find c from g radiation
    if (m_pid == 21 && pid == 4)
      index_first_c_from_g = i;
    if (m_pid == 21 && pid == -4)
      index_first_ac_from_g = i;
  }

  // cout << endl;
  // cout << "test" << index_last_t << " " << index_last_at << " " << index_first_b_from_g << " " << index_first_ab_from_g << " " << index_first_c_from_g << " " << index_first_ac_from_g << endl;
  // if (index_first_b_from_g != -999 || index_first_c_from_g != -999)
  //   cout << "found" << endl;

  bool chk_ttbar_found = false;
  if (index_last_t != -999 && index_last_at != -999)
    chk_ttbar_found = true;

  bool chk_bbbar_from_gluon = false;
  if ((index_first_b_from_g != -999 && index_first_ab_from_g != -999) && (vec_gen[index_first_b_from_g].MotherIndex() == vec_gen[index_first_ab_from_g].MotherIndex()))
    chk_bbbar_from_gluon = true;

  bool chk_ccbar_from_gluon = false;
  if ((index_first_c_from_g != -999 && index_first_ac_from_g != -999) && (vec_gen[index_first_c_from_g].MotherIndex() == vec_gen[index_first_ac_from_g].MotherIndex()))
    chk_ccbar_from_gluon = true;

  if (chk_ttbar_found == false)
    return -1;
  else
  {
    if (chk_bbbar_from_gluon)
      return 2;
    else if (chk_ccbar_from_gluon)
      return 1;
    else
      return 0;
  }

  return -999;
} // int Check_Process(cont vector<Gen> &vec_gen)

//////////

void Vcb_Tagging_RF_DL::Clear()
{
  vec_gen_hf_flavour.clear();
  vec_gen_hf_origin.clear();

  weight = 1;

  weight_b_tag = 1;
  weight_c_tag = 1;

  weight_hem_veto = 1;
  weight_lumi = 1;
  weight_mc = 1;

  weight_pileup = 1;
  weight_pujet_veto = 1;
  weight_prefire = 1;
  weight_top_pt = 1;

  weight_mu_id = 1;
  weight_mu_iso = 1;

  weight_el_id = 1;
  weight_el_reco = 1;

  weight_sl_trig = 1;

  vec_this_muon.clear();
  vec_this_electron.clear();
  vec_this_jet.clear();

  vec_sel_jet.clear();

  vec_btag.clear();
  vec_ctag.clear();

  return;
} // void Vcb_Tagging_RF_DL::Clear()

//////////

Particle Vcb_Tagging_RF_DL::Rebalance_Met()
{
  Particle met_rebal;
  met_rebal += met;

  // before
  for (auto &muon : vec_muon)
    met_rebal += muon;
  for (auto &electron : vec_electron)
    met_rebal += electron;
  for (auto &jet : vec_jet)
    met_rebal += jet;

  // after
  for (auto &muon : vec_this_muon)
    met_rebal -= muon;
  for (auto &electron : vec_this_electron)
    met_rebal -= electron;
  for (auto &jet : vec_this_jet)
    met_rebal -= jet;

  return met_rebal;
} // Particle Vcb_Tagging_RF_DL::Rebalance_Met()

//////////

void Vcb_Tagging_RF_DL::Set_Result_Tree()
{
  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    AnalyzerParameter::Syst syst_type = vec_syst_type.at(i);

    TTree *result_tree = new TTree("Result_Tree", "Result_Tree");

    if (syst_type == AnalyzerParameter::Central)
    {
      result_tree->Branch("weight_b_tag", &weight_b_tag);
      result_tree->Branch("weight_b_tag_down_hf", &weight_b_tag_down_hf);
      result_tree->Branch("weight_b_tag_up_hf", &weight_b_tag_up_hf);
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
    }
    else if (syst_type == AnalyzerParameter::JetEnDown)
      result_tree->Branch("weight_b_tag_down_jes", &weight_b_tag_down_jes);
    else if (syst_type == AnalyzerParameter::JetEnUp)
      result_tree->Branch("weight_b_tag_up_jes", &weight_b_tag_up_jes);

    if (syst_type == AnalyzerParameter::Central)
    {
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
    }
    else if (syst_type == AnalyzerParameter::JetEnDown)
      result_tree->Branch("weight_c_tag_down_jes_total", &weight_c_tag_down_jes_total);
    else if (syst_type == AnalyzerParameter::JetEnUp)
      result_tree->Branch("weight_c_tag_up_jes_total", &weight_c_tag_up_jes_total);
    else if (syst_type == AnalyzerParameter::JetResDown)
      result_tree->Branch("weight_c_tag_down_jer", &weight_c_tag_down_jer);
    else if (syst_type == AnalyzerParameter::JetResUp)
      result_tree->Branch("weight_c_tag_up_jer", &weight_c_tag_up_jer);

    result_tree->Branch("weight_el_id", &weight_el_id);
    result_tree->Branch("weight_el_reco", &weight_el_reco);

    result_tree->Branch("weight_mu_id", &weight_mu_id);
    result_tree->Branch("weight_mu_iso", &weight_mu_iso);

    result_tree->Branch("weight_sl_trig", &weight_sl_trig);

    result_tree->Branch("weight_hem_veto", &weight_hem_veto);
    result_tree->Branch("weight_lumi", &weight_lumi);
    result_tree->Branch("weight_mc", &weight_mc);
    result_tree->Branch("weight_pileup", &weight_pileup);
    if (syst_type == AnalyzerParameter::Central)
    {
      result_tree->Branch("weight_pileup_down", &weight_pileup_down);
      result_tree->Branch("weight_pileup_up", &weight_pileup_up);
    }

    if (syst_type == AnalyzerParameter::Central)
      result_tree->Branch("weight_ps", weight_ps, "weight_ps[4]");

    if (syst_type == AnalyzerParameter::Central)
    {
      result_tree->Branch("weight_scale_variation_1", &weight_scale_variation_1);
      result_tree->Branch("weight_scale_variation_2", &weight_scale_variation_2);
      result_tree->Branch("weight_scale_variation_3", &weight_scale_variation_3);
      result_tree->Branch("weight_scale_variation_4", &weight_scale_variation_4);
      result_tree->Branch("weight_scale_variation_6", &weight_scale_variation_6);
      result_tree->Branch("weight_scale_variation_8", &weight_scale_variation_8);
    }

    result_tree->Branch("weight_prefire", &weight_prefire);
    result_tree->Branch("weight_pujet_veto", &weight_pujet_veto);
    result_tree->Branch("weight_top_pt", &weight_top_pt);

    result_tree->Branch("n_vertex", &nPV);

    result_tree->Branch("n_jets", &n_sel_jet);
    result_tree->Branch("n_bjets", &n_b_jet);
    result_tree->Branch("n_cjets", &n_c_jet);

    result_tree->Branch("decay_mode", &decay_mode);

    result_tree->Branch("Gen_HF_Flavour", &vec_gen_hf_flavour);
    result_tree->Branch("Gen_HF_Origin", &vec_gen_hf_origin);

    result_tree->Branch("Sel_Gen_HF_Flavour", &vec_sel_gen_hf_flavour);
    result_tree->Branch("Sel_Gen_HF_Origin", &vec_sel_gen_hf_origin);

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

    map_result_tree.insert({syst_type, result_tree});
  } // loop over syst

  return;
} // void Vcb_Tagging_RF_DL::Set_Result_Tree()

//////////
