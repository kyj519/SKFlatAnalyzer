#include "Vcb_Tagging_RF.h"

//////////

Vcb_Tagging_RF::Vcb_Tagging_RF()
{
} // Vcb_Tagging_RF::Vcb_Tagging_RF()

//////////

Vcb_Tagging_RF::~Vcb_Tagging_RF()
{
  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    param.syst_ = vec_syst_type.at(i);

    outfile->cd(param.GetSystType());
    map_result_tree[param.syst_]->Write();
  }
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

  // set sigle lepton trigger
  if (DataYear == 2016)
  {
    if (run_mu_ch)
    {
      vec_mu_trig.push_back("HLT_IsoMu24_v");
      sl_trig = "IsoMu24";
      sl_trig_safe_pt_cut = 26.;
    }
    else if (run_el_ch)
    {
      vec_el_trig.push_back("HLT_Ele27_WPTight_Gsf_v");
      sl_trig = "Ele27";
      sl_trig_safe_pt_cut = 30.;
    }
  } // if (DataYear == 2016)
  else if (DataYear == 2017)
  {
    if (run_mu_ch)
    {
      vec_mu_trig.push_back("HLT_IsoMu27_v");
      sl_trig = "IsoMu27";
      sl_trig_safe_pt_cut = 30.;
    }
    else if (run_el_ch)
    {
      // vec_sl_trig.push_back("HLT_Ele35_WPTight_Gsf_v");
      // sl_trig = "Ele35";
      // sl_trig_safe_pt_cut = 37.;

      vec_el_trig.push_back("HLT_Ele32_WPTight_Gsf_L1DoubleEG_");
      sl_trig = "Ele32";
      sl_trig_safe_pt_cut = 35.;
    }
  } // else if (DataYear == 2017)
  else if (DataYear == 2018)
  {
    if (run_mu_ch)
    {
      vec_mu_trig.push_back("HLT_IsoMu24_v");
      sl_trig = "IsoMu24";
      sl_trig_safe_pt_cut = 26.;
    }
    else if (run_el_ch)
    {
      vec_el_trig.push_back("HLT_Ele32_WPTight_Gsf_v");
      sl_trig = "Ele32";
      sl_trig_safe_pt_cut = 35.;
    }
  } // else if (DataYear == 2018)
  else
    std::runtime_error("No trigger configuration for year");

  for (auto &trigger_name : vec_mu_trig)
    vec_sl_trig.push_back(trigger_name);
  for (auto &trigger_name : vec_el_trig)
    vec_sl_trig.push_back(trigger_name);

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
} // void Vcb_Tagging_RF::initializeAnalyzer()

//////////

void Vcb_Tagging_RF::executeEvent()
{
  // Apply Jet Veto Map
  if (IsEventJetMapVetoed())
    return;

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

    // param.Muon_Tight_ID = "POGTightWithTightIso";
    param.Muon_Tight_ID = "POGTight";
    param.Muon_Loose_ID = "POGLoose";

    param.Muon_ID_SF_Key = "NUM_TightID_DEN_TrackerMuons";
    param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";

    // param.Electron_Tight_ID = "passTightID";
    // param.Electron_Loose_ID = "passLooseID";
    // param.Electron_Tight_ID = "passMVAID_iso_WP80";
    // param.Electron_Loose_ID = "passMVAID_iso_WP90";
    param.Electron_Tight_ID = "passMVAID_noIso_WP80";
    param.Electron_Loose_ID = "passMVAID_noIso_WP90";

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
} // void Vcb_Tagging_RF::executeEvent()

//////////

void Vcb_Tagging_RF::executeEventFromParameter(AnalyzerParameter param)
{
  Clear();

  Event ev = GetEvent();

  if (!IsData)
  {
    vec_gen = GetGens();
    decay_mode = Get_W_Decay_Mode(vec_gen);

    // TTBB or TTCC check using vec_jet
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

  for (unsigned int i = 0; i < vec_sel_jet.size(); i++)
  {
    Jet jet = vec_sel_jet[i];

    vec_sel_gen_hf_flavour.push_back(jet.GenHFHadronMatcherFlavour());
    vec_sel_gen_hf_origin.push_back(jet.GenHFHadronMatcherOrigin());
  }

  weight_hem_veto = Weight_HEM_Veto(vec_sel_jet);

  // single lepton trigger
  if (run_mu_ch)
  {
    if (!ev.PassTrigger(vec_mu_trig))
      return;
  }
  else if (run_el_ch)
  {
    if (!ev.PassTrigger(vec_el_trig) || !HLT_SE_Filter_2017(vec_sel_electron))
      return;
  }

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
  if (run_el_ch)
    lepton_pt_uncorr = electron.UncorrPt();
  lepton_eta = lepton.Eta();
  lepton_rel_iso = lepton.RelIso();

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

    if (DataEra == "2016preVFP")
    {
      cvsl_wp = cvsl_2016a_m;
      cvsb_wp = cvsb_2016a_m;
    }
    else if (DataEra == "2016postVFP")
    {
      cvsl_wp = cvsl_2016b_m;
      cvsb_wp = cvsb_2016b_m;
    }
    else if (DataEra == "2017")
    {
      cvsl_wp = cvsl_2017_m;
      cvsb_wp = cvsb_2017_m;
    }
    else if (DataEra == "2018")
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

      weight_b_tag_down_lf = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "down_lf");
      weight_b_tag_up_lf = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "up_lf");

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
    else
      weight_b_tag = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "central");

    weight *= weight_b_tag;
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
    }
    else if (param.syst_ == AnalyzerParameter::JetEnDown)
      weight_c_tag_down_jes_total = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jesTotal_Down");
    else if (param.syst_ == AnalyzerParameter::JetEnUp)
      weight_c_tag_up_jes_total = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jesTotal_Up");
    else if (param.syst_ == AnalyzerParameter::JetResDown)
      weight_c_tag_down_jer = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jer_Down");
    else if (param.syst_ == AnalyzerParameter::JetResUp)
      weight_c_tag_up_jer = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "jer_Up");
    else
      weight_c_tag = mcCorr->GetCTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(1), "central");

    weight *= weight_c_tag;
  }

  // no b-tag is required
  // if (n_b_jet < 2)
  //  return;

  // cut on MET
  met_pt = met.Pt();
  // met_phi = met.Phi();

  // if (met_pt < MET_PT)
  //   return;

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

Particle Vcb_Tagging_RF::Rebalance_Met()
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
} // Particle Vcb_Tagging_RF::Rebalance_Met()

//////////

void Vcb_Tagging_RF::Set_Result_Tree()
{
  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    AnalyzerParameter::Syst syst_type = vec_syst_type.at(i);

    TTree *result_tree = new TTree("Result_Tree", "Result_Tree");

    result_tree->Branch("sf_mu_id", &sf_mu_id);
    result_tree->Branch("sf_mu_iso", &sf_mu_iso);
    result_tree->Branch("sf_mu_trig", &sf_mu_trig);

    result_tree->Branch("sf_el_id", &sf_el_id);
    result_tree->Branch("sf_el_reco", &sf_el_reco);

    result_tree->Branch("sf_sl_trig", &sf_sl_trig);

    if (syst_type == AnalyzerParameter::Central)
    {
      result_tree->Branch("weight_b_tag", &weight_b_tag);
      result_tree->Branch("weight_b_tag_down_hf", &weight_b_tag_down_hf);
      result_tree->Branch("weight_b_tag_up_hf", &weight_b_tag_up_hf);
      result_tree->Branch("weight_b_tag_down_lf", &weight_b_tag_down_lf);
      result_tree->Branch("weight_b_tag_up_lf", &weight_b_tag_up_lf);
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
    else
      result_tree->Branch("weight_b_tag", &weight_b_tag);

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
    else
      result_tree->Branch("weight_c_tag", &weight_c_tag);

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

    result_tree->Branch("lepton_pt", &lepton_pt);
    result_tree->Branch("lepton_eta", &lepton_eta);
    result_tree->Branch("lepton_pt_uncorr", &lepton_pt_uncorr);
    result_tree->Branch("lepton_rel_iso", &lepton_rel_iso);

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

    result_tree->Branch("met_pt", &met_pt);

    result_tree->Branch("decay_mode", &decay_mode);

    result_tree->Branch("Gen_HF_Flavour", &vec_gen_hf_flavour);
    result_tree->Branch("Gen_HF_Origin", &vec_gen_hf_origin);

    result_tree->Branch("Sel_Gen_HF_Flavour", &vec_sel_gen_hf_flavour);
    result_tree->Branch("Sel_Gen_HF_Origin", &vec_sel_gen_hf_origin);

    map_result_tree.insert({syst_type, result_tree});
  }

  return;
} // void Vcb_Tagging_RF::Set_Result_Tree()

//////////