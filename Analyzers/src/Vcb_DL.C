#include "Vcb_DL.h"

//////////

Vcb_DL::Vcb_DL()
{
} // Vcb_DL::Vcb_DL()

//////////

Vcb_DL::~Vcb_DL()
{
  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    param.syst_ = vec_syst_type.at(i);

    outfile->cd(param.GetSystType());
    map_result_tree[param.syst_]->Write();
  }

} // Vcb_DL::~Vcb_DL()

//////////

void Vcb_DL::initializeAnalyzer()
{
  run_debug = HasFlag("RunDebug");
  cout << "[Vcb_DL::initializeAnalyzer] RunDebug = " << run_debug << endl;

  run_mm = HasFlag("RunMM");
  cout << "[Vcb_DL::initializeAnalyzer] RunMM = " << run_mm << endl;

  run_me = HasFlag("RunME");
  cout << "[Vcb_DL::initializeAnalyzer] RunME = " << run_me << endl;

  run_ee = HasFlag("RunEE");
  cout << "[Vcb_DL::initializeAnalyzer] RunEE = " << run_ee << endl;

  int chk_ch = run_mm + run_me + run_ee;
  if (chk_ch != 1)
  {
    cout << "One of RunMM, RunME, or RunEE should be set." << endl;
    exit(1);
  }

  run_syst = HasFlag("RunSyst");
  cout << "[Vcb_DL::initializeAnalyzer] RunSyst = " << run_syst << endl;

  // set muon object
  vec_mu_id = {"POGTightWithTightIso"};
  vec_mu_id_sf_key = {"NUM_TightID_DEN_TrackerMuons"};
  vec_mu_iso_sf_key = {"NUM_TightRelIso_DEN_TightIDandIPCut"};

  // set electron  object
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

  // PDF syst
  if (run_syst)
  {
    LHAPDFHandler LHAPDFHandler_Prod;
    // LHAPDFHandler_Prod.CentralPDFName = "NNPDF31_nnlo_hessian_pdfas";
    LHAPDFHandler_Prod.CentralPDFName = "NNPDF31_nnlo_as_0118_mc_hessian_pdfas";
    LHAPDFHandler_Prod.init();

    LHAPDFHandler LHAPDFHandler_New;
    LHAPDFHandler_New.CentralPDFName = "NNPDF31_nlo_hessian_pdfas";
    LHAPDFHandler_New.ErrorSetMember_Start = 1;
    LHAPDFHandler_New.ErrorSetMember_End = 100;
    LHAPDFHandler_New.AlphaSMember_Down = 101;
    LHAPDFHandler_New.AlphaSMember_Up = 102;
    LHAPDFHandler_New.init();

    pdfReweight->SetProdPDF(LHAPDFHandler_Prod.PDFCentral);
    pdfReweight->SetNewPDF(LHAPDFHandler_New.PDFCentral);
    pdfReweight->SetNewPDFErrorSet(LHAPDFHandler_New.PDFErrorSet);
    pdfReweight->SetNewPDFAlphaS(LHAPDFHandler_New.PDFAlphaSDown, LHAPDFHandler_New.PDFAlphaSUp);
  }

  if (!IsData && run_syst)
  {
    vec_syst_type = {AnalyzerParameter::Central,
                     AnalyzerParameter::JetEnDown,
                     AnalyzerParameter::JetEnUp,
                     AnalyzerParameter::JetResDown,
                     AnalyzerParameter::JetResUp,
                     AnalyzerParameter::UnclusteredEnergyDown,
                     AnalyzerParameter::UnclusteredEnergyUp};
    if (run_me || run_ee)
    {
      vec_syst_type.push_back(AnalyzerParameter::ElectronEnDown);
      vec_syst_type.push_back(AnalyzerParameter::ElectronEnUp);
      vec_syst_type.push_back(AnalyzerParameter::ElectronResDown);
      vec_syst_type.push_back(AnalyzerParameter::ElectronResUp);
    }
  }
  else
    vec_syst_type = {AnalyzerParameter::Central};

  // to make output dir
  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    param.syst_ = vec_syst_type.at(i);
    FillHist(param.GetSystType() + "/Dummy", 0, weight, 1, 0, 1);
  }

  Set_Result_Tree();

  return;
} // void Vcb_DL::initializeAnalyzer()

//////////

void Vcb_DL::executeEvent()
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
    // param setting
    param.Clear();

    param.Muon_Tight_ID = "POGTightWithTightIso";
    param.Muon_Loose_ID = "POGLoose";
    param.Muon_ID_SF_Key = "NUM_TightID_DEN_TrackerMuons";
    param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";

    // param.Electron_Tight_ID = "passTightID";
    // param.Electron_Loose_ID = "passLooseID";

    param.Electron_Tight_ID = "passMVAID_iso_WP80";
    param.Electron_Loose_ID = "passMVAID_iso_WP90";

    // param.Electron_ID_SF_Key = "ID_SF_passTightID";

    param.Jet_ID = "tight";
    param.PUJet_Veto_ID = "LoosePileupJetVeto";

    param.syst_ = vec_syst_type.at(i);

    executeEventFromParameter(param);
  }

  return;
} // void Vcb_DL::executeEvent()

//////////

void Vcb_DL::executeEventFromParameter(AnalyzerParameter param)
{
  Clear();

  Event ev = GetEvent();

  if (!IsData)
  {
    vec_gen = GetGens();

    process = Check_Process(vec_gen);
    // cout << "test process = " << process << endl;

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

    // PDF
    if (param.syst_ == AnalyzerParameter::Central && run_syst)
    {
      weight_pdf_alternative = GetPDFReweight();
      for (int i = 0; i < 100; i++)
        weight_pdf_error_set[i] = GetPDFReweight(i);
      weight_pdf_as_up = GetPDFReweight("As_Up");
      weight_pdf_as_down = GetPDFReweight("As_Down");
    }

    // pileup reweight
    weight_pileup = mcCorr->GetPileUpWeight(nPileUp, 0);
    if (param.syst_ == AnalyzerParameter::Central && run_syst)
    {
      weight_pileup_down = mcCorr->GetPileUpWeight(nPileUp, -1);
      weight_pileup_up = mcCorr->GetPileUpWeight(nPileUp, +1);
    }
    weight *= weight_pileup;

    // L1 prefire
    weight_prefire = GetPrefireWeight(0);
    if (param.syst_ == AnalyzerParameter::Central && run_syst)
    {
      weight_prefire_down = GetPrefireWeight(-1);
      weight_prefire_up = GetPrefireWeight(+1);
    }
    weight *= weight_prefire;

    // Top Pt reweight
    weight_top_pt = mcCorr->GetTopPtReweight(vec_gen);
    weight *= weight_top_pt;

    // Scale Variation
    if (param.syst_ == AnalyzerParameter::Central && run_syst)
    {
      weight_scale_variation_1 = GetScaleVariation(1);
      weight_scale_variation_2 = GetScaleVariation(2);
      weight_scale_variation_3 = GetScaleVariation(3);
      weight_scale_variation_4 = GetScaleVariation(4);
      weight_scale_variation_6 = GetScaleVariation(6);
      weight_scale_variation_8 = GetScaleVariation(8);

      // PS Reweight
      Get_Reweight_PS(weight_ps);
    }
  } //  if (!IsData)

  // met filter
  if (!PassMETFilter())
    return;

  // set objects
  vec_this_muon = vec_muon;
  vec_this_electron = vec_electron;
  vec_this_jet = vec_jet;

  if (param.syst_ == AnalyzerParameter::UnclusteredEnergyDown)
    Met_Syst_Unclustered(met, -1);
  else if (param.syst_ == AnalyzerParameter::UnclusteredEnergyUp)
    Met_Syst_Unclustered(met, 1);
  else
    met = ev.GetMETVector("PUPPI");

  pair<double, double> met_corr = xy_met_correction.METXYCorr_Met_MetPhi(met.Pt(), met.Phi(), run, to_string(DataYear), !IsData, nPV, true, true);
  met.SetPtEtaPhiE(met_corr.first, 0, met_corr.second, met_corr.first);

  // syst for objects
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
  // if(param.syst_ == AnalyzerParameter::MuonEnUp) vec_this_muon = ScaleMuons(vec_this_muon, +1);
  // if(param.syst_ == AnalyzerParameter::MuonEnDown) vec_this_muon = ScaleMuons(vec_this_muon, -1);

  if (param.syst_ == AnalyzerParameter::ElectronEnDown)
  {
    vec_this_electron = ScaleElectrons(vec_this_electron, -1);
    met = Rebalance_Met();
  }
  else if (param.syst_ == AnalyzerParameter::ElectronEnUp)
  {
    vec_this_electron = ScaleElectrons(vec_this_electron, +1);
    met = Rebalance_Met();
  }

  if (param.syst_ == AnalyzerParameter::ElectronResDown)
  {
    vec_this_electron = SmearElectrons(vec_this_electron, -1);
    met = Rebalance_Met();
  }
  else if (param.syst_ == AnalyzerParameter::ElectronResUp)
  {
    vec_this_electron = SmearElectrons(vec_this_electron, +1);
    met = Rebalance_Met();
  }

  /////////////////
  /*Setup Objects*/
  /////////////////

  // for lepton selection
  vec_sel_muon = SelectMuons(vec_this_muon, param.Muon_Tight_ID, mu_trig_safe_pt_cut, MUON_ETA);
  vec_sel_electron = SelectElectrons(vec_this_electron, param.Electron_Tight_ID, el_trig_safe_pt_cut, ELECTRON_ETA);

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

  // sort jet as pt ordering
  sort(vec_sel_jet.begin(), vec_sel_jet.end(), PtComparing);

  weight_hem_veto = Weight_HEM_Veto(vec_sel_jet);

  // single lepton trigger
  if (IsData && run_me && (DataStream == "EGamma" || DataStream == "SingleElectron"))
  {
    if (ev.PassTrigger(vec_mu_trig) || (!ev.PassTrigger(vec_el_trig) || !HLT_SE_Filter_2017(vec_sel_electron)))
      return;
  }
  else
  {
    if (!ev.PassTrigger(vec_mu_trig) && (!ev.PassTrigger(vec_el_trig) || !HLT_SE_Filter_2017(vec_sel_electron)))
      return;
  }

  if (!IsData)
  {
    weight_sl_trig = mcCorr->SingleLepton_Trigger_SF("POGTight", mu_trig, vec_sel_muon, 0, "passTightID", el_trig, vec_sel_electron, 0);
    if (param.syst_ == AnalyzerParameter::Central && run_syst)
    {
      weight_sl_trig_el_down = mcCorr->SingleLepton_Trigger_SF("POGTight", mu_trig, vec_sel_muon, 0, "passTightID", el_trig, vec_sel_electron, -1);
      weight_sl_trig_el_up = mcCorr->SingleLepton_Trigger_SF("POGTight", mu_trig, vec_sel_muon, 0, "passTightID", el_trig, vec_sel_electron, 1);
      weight_sl_trig_mu_down = mcCorr->SingleLepton_Trigger_SF("POGTight", mu_trig, vec_sel_muon, -1, "passTightID", el_trig, vec_sel_electron, 0);
      weight_sl_trig_mu_up = mcCorr->SingleLepton_Trigger_SF("POGTight", mu_trig, vec_sel_muon, 1, "passTightID", el_trig, vec_sel_electron, 0);
    }

    weight *= weight_sl_trig;

    // muon reweight
    for (unsigned int i = 0; i < vec_sel_muon.size(); i++)
    {
      weight_mu_id *= mcCorr->MuonID_SF(param.Muon_ID_SF_Key, vec_sel_muon[i].Eta(), vec_sel_muon[i].MiniAODPt(), 0);
      if (param.syst_ == AnalyzerParameter::Central && run_syst)
      {
        weight_mu_id_down *= mcCorr->MuonID_SF(param.Muon_ID_SF_Key, vec_sel_muon[i].Eta(), vec_sel_muon[i].MiniAODPt(), -1);
        weight_mu_id_up *= mcCorr->MuonID_SF(param.Muon_ID_SF_Key, vec_sel_muon[i].Eta(), vec_sel_muon[i].MiniAODPt(), 1);
      }

      weight_mu_iso *= mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, vec_sel_muon[i].Eta(), vec_sel_muon[i].MiniAODPt(), 0);
      if (param.syst_ == AnalyzerParameter::Central && run_syst)
      {
        weight_mu_iso_down *= mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, vec_sel_muon[i].Eta(), vec_sel_muon[i].MiniAODPt(), -1);
        weight_mu_iso_up *= mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, vec_sel_muon[i].Eta(), vec_sel_muon[i].MiniAODPt(), 1);
      }
    }
    weight *= weight_mu_id;
    weight *= weight_mu_iso;

    // electron reweight
    for (unsigned int i = 0; i < vec_sel_electron.size(); i++)
    {
      weight_el_id *= mcCorr->ElectronID_SF(param.Electron_Tight_ID, vec_sel_electron[i].scEta(), vec_sel_electron[i].UncorrPt(), 0);
      if (param.syst_ == AnalyzerParameter::Central && run_syst)
      {
        weight_el_id_down *= mcCorr->ElectronID_SF(param.Electron_Tight_ID, vec_sel_electron[i].scEta(), vec_sel_electron[i].UncorrPt(), -1);
        weight_el_id_up *= mcCorr->ElectronID_SF(param.Electron_Tight_ID, vec_sel_electron[i].scEta(), vec_sel_electron[i].UncorrPt(), 1);
      }

      weight_el_reco *= mcCorr->ElectronReco_SF(vec_sel_electron[i].scEta(), vec_sel_electron[i].UncorrPt(), 0);
      if (param.syst_ == AnalyzerParameter::Central && run_syst)
      {
        weight_el_reco_down *= mcCorr->ElectronReco_SF(vec_sel_electron[i].scEta(), vec_sel_electron[i].UncorrPt(), -1);
        weight_el_reco_up *= mcCorr->ElectronReco_SF(vec_sel_electron[i].scEta(), vec_sel_electron[i].UncorrPt(), 1);
      }
    }
    weight *= weight_el_id;
    weight *= weight_el_reco;
  }

  // cut on double lepton and veto additional loose electron
  if (run_mm)
  {
    if (vec_sel_muon.size() != 2 || vec_muon_veto.size() != 2 || vec_electron_veto.size() != 0)
      return;

    lepton[0] = vec_sel_muon[0];
    lepton[1] = vec_sel_muon[1];
  }
  else if (run_me)
  {
    if (vec_sel_muon.size() != 1 || vec_muon_veto.size() != 1 || vec_sel_electron.size() != 1 || vec_electron_veto.size() != 1)
      return;

    lepton[0] = vec_sel_muon[0];
    lepton[1] = vec_sel_electron[0];
  }
  else if (run_ee)
  {
    if (vec_sel_muon.size() != 0 || vec_sel_electron.size() != 2 || vec_electron_veto.size() != 2)
      return;

    lepton[0] = vec_sel_electron[0];
    lepton[1] = vec_sel_electron[1];
  }

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

  // cut on dilepton mass to veto low mass resonance
  TLorentzVector dilepton = lepton[0];
  dilepton += lepton[1];
  dilepton_mass = dilepton.M();
  if (dilepton_mass < 15)
    return;

  // cut on dilepton mass to veto Z
  // if (abs(dilepton_mass - Z_MASS) < 15)
  //  return;

  // cut on jet
  n_sel_jet = vec_sel_jet.size();
  if (n_sel_jet < 4)
    return;

  if (!IsData)
  {
    // SF for PU jet veto
    weight_pujet_veto = mcCorr->PileupJetVeto_Reweight(vec_sel_jet, param.PUJet_Veto_ID, 0);
    if (param.syst_ == AnalyzerParameter::Central && run_syst)
    {
      weight_pujet_veto_down = mcCorr->PileupJetVeto_Reweight(vec_sel_jet, param.PUJet_Veto_ID, -1);
      weight_pujet_veto_up = mcCorr->PileupJetVeto_Reweight(vec_sel_jet, param.PUJet_Veto_ID, 1);
    }
  }

  // n_b_jet
  for (auto &jet : vec_sel_jet)
  {
    float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
    if (mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score)
      n_b_jet++;
  }

  // cut on n_b_jet
  if (n_b_jet < 2)
    return;

  if (!IsData)
  {
    // SF for b-tagging
    if (param.syst_ == AnalyzerParameter::Central && run_syst)
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
    else
      weight_b_tag = mcCorr->GetBTaggingReweight_1d(vec_sel_jet, vec_jet_tagging_para.at(0), "central");

    weight *= weight_b_tag;

    // SF for c-tagging

    if (param.syst_ == AnalyzerParameter::Central && run_syst)
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
  } // if (!IsData)

  met_pt = met.Pt();
  met_phi = met.Phi();

  // cut on MET
  if (!run_me && met_pt < MET_PT_DL)
    return;

  Make_Result_Tree(param);

  return;
} // Vcb_DL::executeEventFromParameter(AnalyzerParameter param)

//////////

float Vcb_DL::Calculate_HT(const vector<Jet> &vec_jet)
{
  float ht = 0;

  for (unsigned int i = 0; i < vec_jet.size(); i++)
    ht += vec_jet[i].Pt();

  return ht;
} // float Vcb_DL::Calculate_HT(const vector<Jet> &vec_jet)

//////////

int Vcb_DL::Check_Process(const vector<Gen> &vec_gen)
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
} // int Vcb_DL::Check_Process(const vector<Gen>& vec_gen)

//////////

void Vcb_DL::Clear()
{
  n_sel_jet = 0;
  n_b_jet = 0;

  process = -999;

  vec_gen_hf_flavour.clear();
  vec_gen_hf_origin.clear();

  bvsc_third = -9999;
  bvsc_fourth = -9999;

  weight = 1;

  weight_b_tag = 1;
  weight_b_tag_down_hf = 1;
  weight_b_tag_up_hf = 1;
  weight_b_tag_down_jes = 1;
  weight_b_tag_up_jes = 1;
  weight_b_tag_down_lfstats1 = 1;
  weight_b_tag_up_lfstats1 = 1;
  weight_b_tag_down_lfstats2 = 1;
  weight_b_tag_up_lfstats2 = 1;
  weight_b_tag_down_cferr1 = 1;
  weight_b_tag_up_cferr1 = 1;
  weight_b_tag_down_cferr2 = 1;
  weight_b_tag_up_cferr2 = 1;
  weight_b_tag_down_hfstats1 = 1;
  weight_b_tag_up_hfstats1 = 1;
  weight_b_tag_down_hfstats2 = 1;
  weight_b_tag_up_hfstats2 = 1;

  weight_c_tag = 1;
  weight_c_tag_down_extrap = 1;
  weight_c_tag_up_extrap = 1;
  weight_c_tag_down_interp = 1;
  weight_c_tag_up_interp = 1;
  weight_c_tag_down_lhe_scale_muf = 1;
  weight_c_tag_up_lhe_scale_muf = 1;
  weight_c_tag_down_lhe_scale_mur = 1;
  weight_c_tag_up_lhe_scale_mur = 1;
  weight_c_tag_down_ps_fsr_fixed = 1;
  weight_c_tag_up_ps_fsr_fixed = 1;
  weight_c_tag_down_ps_isr_fixed = 1;
  weight_c_tag_up_ps_isr_fixed = 1;
  weight_c_tag_down_pu = 1;
  weight_c_tag_up_pu = 1;
  weight_c_tag_down_stat = 1;
  weight_c_tag_up_stat = 1;
  weight_c_tag_down_xsec_brunc_dyjets_b = 1;
  weight_c_tag_up_xsec_brunc_dyjets_b = 1;
  weight_c_tag_down_xsec_brunc_dyjets_c = 1;
  weight_c_tag_up_xsec_brunc_dyjets_c = 1;
  weight_c_tag_down_xsec_brunc_wjets_c = 1;
  weight_c_tag_up_xsec_brunc_wjets_c = 1;
  weight_c_tag_down_jer = 1;
  weight_c_tag_up_jer = 1;
  weight_c_tag_down_jes_total = 1;
  weight_c_tag_up_jes_total = 1;

  weight_el_id = 1;
  weight_el_id_down = 1;
  weight_el_id_up = 1;

  weight_el_reco = 1;
  weight_el_reco_down = 1;
  weight_el_reco_up = 1;

  weight_lumi = 1;
  weight_mc = 1;

  weight_mu_id = 1;
  weight_mu_id_down = 1;
  weight_mu_id_up = 1;

  weight_mu_iso = 1;
  weight_mu_iso_down = 1;
  weight_mu_iso_up = 1;

  weight_pdf_alternative = 1;
  weight_pdf_error_set[100] = {1};
  weight_pdf_as_down = 1;
  weight_pdf_as_up = 1;

  weight_pileup = 1;
  weight_pileup_down = 1;
  weight_pileup_up = 1;

  weight_prefire = 1;
  weight_prefire_down = 1;
  weight_prefire_up = 1;

  weight_pujet_veto = 1;
  weight_pujet_veto_down = 1;
  weight_pujet_veto_up = 1;

  weight_scale_variation_1 = 1;
  weight_scale_variation_2 = 1;
  weight_scale_variation_3 = 1;
  weight_scale_variation_4 = 1;
  weight_scale_variation_6 = 1;
  weight_scale_variation_8 = 1;

  weight_sl_trig = 1;
  weight_sl_trig_el_down = 1;
  weight_sl_trig_el_up = 1;
  weight_sl_trig_mu_down = 1;
  weight_sl_trig_mu_up = 1;

  weight_top_pt = 1;

  return;
} // void Vcb_DL::Clear()

//////////

void Vcb_DL::Make_Result_Tree(const AnalyzerParameter &param)
{
  leading_lepton_eta = lepton[0].Eta();
  leading_lepton_pt = lepton[0].Pt();

  subleading_lepton_eta = lepton[1].Eta();
  subleading_lepton_pt = lepton[1].Pt();

  ht = Calculate_HT(vec_sel_jet);

  leading_jet_bvsc = vec_sel_jet[0].GetTaggerResult(JetTagging::DeepJet);
  leading_jet_cvsb = vec_sel_jet[0].GetTaggerResult(JetTagging::DeepJet_CvsB);
  leading_jet_cvsl = vec_sel_jet[0].GetTaggerResult(JetTagging::DeepJet_CvsL);

  subleading_jet_bvsc = vec_sel_jet[1].GetTaggerResult(JetTagging::DeepJet);
  subleading_jet_cvsb = vec_sel_jet[1].GetTaggerResult(JetTagging::DeepJet_CvsB);
  subleading_jet_cvsl = vec_sel_jet[1].GetTaggerResult(JetTagging::DeepJet_CvsL);

  sort(vec_sel_jet.begin(), vec_sel_jet.end(), BvsCComparing);

  bvsc_third = vec_sel_jet[2].GetTaggerResult(JetTagging::DeepJet);
  bvsc_fourth = vec_sel_jet[3].GetTaggerResult(JetTagging::DeepJet);

  map_result_tree[param.syst_]->Fill();

  return;
} // void Vcb_DL::Make_Result_Tree()

/////////

Particle Vcb_DL::Rebalance_Met()
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
} // Particle Vcb_DL::Rebalance_Met()

//////////

void Vcb_DL::Set_Result_Tree()
{

  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    AnalyzerParameter::Syst syst_type = vec_syst_type.at(i);

    TTree *result_tree = new TTree("Result_Tree", "Result_Tree");

    result_tree->Branch("n_pv", &nPV);

    result_tree->Branch("leading_lepton_eta", &leading_lepton_eta);
    result_tree->Branch("leading_lepton_pt", &leading_lepton_pt);

    result_tree->Branch("subleading_lepton_eta", &subleading_lepton_eta);
    result_tree->Branch("subleading_lepton_pt", &subleading_lepton_pt);

    result_tree->Branch("dilepton_mass", &dilepton_mass);

    result_tree->Branch("n_jet", &n_sel_jet);
    result_tree->Branch("n_b_jet", &n_b_jet);
    result_tree->Branch("ht", &ht);

    result_tree->Branch("leading_jet_bvsc", &leading_jet_bvsc);
    result_tree->Branch("leading_jet_cvsb", &leading_jet_cvsb);
    result_tree->Branch("leading_jet_cvsl", &leading_jet_cvsl);

    result_tree->Branch("subleading_jet_bvsc", &subleading_jet_bvsc);
    result_tree->Branch("subleading_jet_cvsb", &subleading_jet_cvsb);
    result_tree->Branch("subleading_jet_cvsl", &subleading_jet_cvsl);

    result_tree->Branch("met_pt", &met_pt);
    result_tree->Branch("met_phi", &met_phi);

    result_tree->Branch("process", &process);

    result_tree->Branch("bvsc_third", &bvsc_third);
    result_tree->Branch("bvsc_fourth", &bvsc_fourth);

    result_tree->Branch("Gen_HF_Flavour", &vec_gen_hf_flavour);
    result_tree->Branch("Gen_HF_Origin", &vec_gen_hf_origin);

    if (syst_type == AnalyzerParameter::Central && run_syst)
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
    else
      result_tree->Branch("weight_b_tag", &weight_b_tag);

    if (syst_type == AnalyzerParameter::Central && run_syst)
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

    result_tree->Branch("weight_el_id", &weight_el_id);
    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_el_id_down", &weight_el_id_down);
      result_tree->Branch("weight_el_id_up", &weight_el_id_up);
    }

    result_tree->Branch("weight_el_reco", &weight_el_reco);
    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_el_reco_down", &weight_el_reco_down);
      result_tree->Branch("weight_el_reco_up", &weight_el_reco_up);
    }

    result_tree->Branch("weight_hem_veto", &weight_hem_veto);
    result_tree->Branch("weight_lumi", &weight_lumi);
    result_tree->Branch("weight_mc", &weight_mc);

    result_tree->Branch("weight_mu_id", &weight_mu_id);
    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_mu_id_down", &weight_mu_id_down);
      result_tree->Branch("weight_mu_id_up", &weight_mu_id_up);
    }

    result_tree->Branch("weight_mu_iso", &weight_mu_iso);
    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_mu_iso_down", &weight_mu_iso_down);
      result_tree->Branch("weight_mu_iso_up", &weight_mu_iso_up);
    }

    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_pdf_alternative", &weight_pdf_alternative);
      result_tree->Branch("weight_pdf_error_set", weight_pdf_error_set, "weight_pdf_error_set[100]/F");
      result_tree->Branch("weight_pdf_as_down", &weight_pdf_as_down);
      result_tree->Branch("weight_pdf_as_up", &weight_pdf_as_up);
    }

    result_tree->Branch("weight_pileup", &weight_pileup);
    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_pileup_down", &weight_pileup_down);
      result_tree->Branch("weight_pileup_up", &weight_pileup_up);
    }

    result_tree->Branch("weight_prefire", &weight_prefire);
    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_prefire_down", &weight_prefire_down);
      result_tree->Branch("weight_prefire_up", &weight_prefire_up);
    }

    if (syst_type == AnalyzerParameter::Central && run_syst)
      result_tree->Branch("weight_ps", weight_ps, "weight_ps[4]");

    result_tree->Branch("weight_pujet_veto", &weight_pujet_veto);
    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_pujet_veto_down", &weight_pujet_veto_down);
      result_tree->Branch("weight_pujet_veto_up", &weight_pujet_veto_up);
    }

    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_scale_variation_1", &weight_scale_variation_1);
      result_tree->Branch("weight_scale_variation_2", &weight_scale_variation_2);
      result_tree->Branch("weight_scale_variation_3", &weight_scale_variation_3);
      result_tree->Branch("weight_scale_variation_4", &weight_scale_variation_4);
      result_tree->Branch("weight_scale_variation_6", &weight_scale_variation_6);
      result_tree->Branch("weight_scale_variation_8", &weight_scale_variation_8);
    }

    result_tree->Branch("weight_sl_trig", &weight_sl_trig);
    if (syst_type == AnalyzerParameter::Central && run_syst)
    {
      result_tree->Branch("weight_sl_trig_el_down", &weight_sl_trig_el_down);
      result_tree->Branch("weight_sl_trig_el_up", &weight_sl_trig_el_up);
      result_tree->Branch("weight_sl_trig_mu_down", &weight_sl_trig_mu_down);
      result_tree->Branch("weight_sl_trig_mu_up", &weight_sl_trig_mu_up);
    }

    result_tree->Branch("weight_top_pt", &weight_top_pt);

    map_result_tree.insert({syst_type, result_tree});
  }

  return;
} // void Vcb_DL::Set_Permutation_Tree()

//////////
