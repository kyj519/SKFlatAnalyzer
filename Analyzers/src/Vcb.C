#include "Vcb.h"

//////////

Vcb::Vcb()
{
  permutation_tree_correct = NULL;
  permutation_tree_wrong = NULL;
} // Vcb::Vcb()

//////////

Vcb::~Vcb()
{
  if (run_permutation_tree)
  {
    outfile->cd();
    permutation_tree_correct->Write();
    permutation_tree_wrong->Write();
  }

  if (run_hf_contamination_tree)
  {
    outfile->cd();
    hf_contamination_tree_correct->Write();
    hf_contamination_tree_wrong->Write();
  }

  if (run_template)
  {
    outfile->cd();

    for (int i = 0; i < 5; i++)
    {
      template_tree[i]->Write();
    }
  }

  if (run_template_truth)
  {
    outfile->cd();

    for (int i = 0; i < 5; i++)
    {
      template_truth_tree[i]->Write();
    }
  }

  if (run_result)
  {
    for (unsigned int i = 0; i < vec_syst_type.size(); i++)
    {
      param.syst_ = vec_syst_type.at(i);
      outfile->cd(param.GetSystType());
      map_result_tree[param.syst_]->Write();
    }
  }

  outfile->Close();

  delete fitter_driver;
  if (!run_permutation_tree)
  {
    delete reader_swapper[0];
    delete reader_swapper[1];
  }
} // Vcb::~Vcb()

//////////

void Vcb::initializeAnalyzer()
{
  run_mu_ch = HasFlag("RunMu");
  cout << "[Vcb::initializeAnalyzer] RunMu = " << run_mu_ch << endl;

  run_el_ch = HasFlag("RunEl");
  cout << "[Vcb::initializeAnalyzer] RunEl = " << run_el_ch << endl;

  if (run_mu_ch == run_el_ch)
  {
    cout << "One of RunMu or RunEl should be set." << endl;
    exit(1);
  }

  if (run_mu_ch)
    channel = "Muon";
  else if (run_el_ch)
    channel = "Electron";

  run_debug = HasFlag("RunDebug");
  cout << "[Vcb::initializeAnalyzer] RunDebug = " << run_debug << endl;

  run_permutation_tree = HasFlag("RunPermutationTree");
  cout << "[Vcb::initializeAnalyzer] RunPermutationTree = " << run_permutation_tree << endl;

  run_hf_contamination_tree = HasFlag("RunHFContaminationTree");
  cout << "[Vcb::initializeAnalyzer] RunHFContaminationTree = " << run_hf_contamination_tree << endl;

  run_chi = HasFlag("RunChi");
  cout << "[Vcb::initializeAnalyzer] RunChi = " << run_chi << endl;

  run_template_truth = HasFlag("RunTemplateTruth");
  cout << "[Vcb::initializeAnalyzer] RuTemplateTruth = " << run_template_truth << endl;

  run_template = HasFlag("RunTemplate");
  cout << "[Vcb::initializeAnalyzer] RunTemplate = " << run_template << endl;

  run_result = HasFlag("RunResult");
  cout << "[Vcb::initializeAnalyzer] RunResult = " << run_result << endl;

  run_syst = HasFlag("RunSyst");
  cout << "[Vcb::initializeAnalyzer] RunSyst = " << run_syst << endl;

  // to check additional figure of merit using hadronic w mass constrain to suppress background
  rm_wm_constraint = HasFlag("RM_WM");
  cout << "[Vcb::initializeAnalyzer] RM_WM = " << rm_wm_constraint << endl;

  rm_bjet_energy_reg_nn = HasFlag("RM_REG_NN");
  cout << "[Vcb::initializeAnalyzer] RM_REG_NN = " << rm_bjet_energy_reg_nn << endl;

  // set single muon object
  vec_mu_id = {"POGTightWithTightIso"};
  vec_mu_id_sf_key = {"NUM_TightID_DEN_TrackerMuons"};
  vec_mu_iso_sf_key = {"NUM_TightRelIso_DEN_TightIDandIPCut"};

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

  for (unsigned int i = 0; i < vec_mu_trig.size(); i++)
    vec_sl_trig.push_back(vec_mu_trig[i]);

  for (unsigned int i = 0; i < vec_el_trig.size(); i++)
    vec_sl_trig.push_back(vec_el_trig[i]);

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

  // Retrieve POG JER
  std::string jetPtResolutionPath = "";
  std::string jetPtResolutionSFPath = "";
  std::string BASE_PATH = std::getenv("SKFlat_WD") + std::string("/data/Run2UltraLegacy_v3/");

  if (DataEra == "2016preVFP")
  {
    jetPtResolutionPath = BASE_PATH + "2016preVFP/JME/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt";
    jetPtResolutionSFPath = BASE_PATH + "2016preVFP/JME/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt";
  }
  else if (DataEra == "2016postVFP")
  {
    jetPtResolutionPath = BASE_PATH + "2016postVFP/JME/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt";
    jetPtResolutionSFPath = BASE_PATH + "2016postVFP/JME/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt";
  }
  else if (DataEra == "2017")
  {
    jetPtResolutionPath = BASE_PATH + "2017/JME/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt";
    jetPtResolutionSFPath = BASE_PATH + "2017/JME/Summer19UL17_JRV2_MC_SF_AK4PFchs.txt";
  }
  else if (DataEra == "2018")
  {
    jetPtResolutionPath = BASE_PATH + "2018/JME/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt";
    jetPtResolutionSFPath = BASE_PATH + "2018/JME/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt";
  }
  else
  {
    throw std::runtime_error("No JME configuration for year");
  }
  jet_resolution = JME::JetResolution(jetPtResolutionPath);
  jet_resolution_sf = JME::JetResolutionScaleFactor(jetPtResolutionSFPath);

  fitter_driver = new TKinFitterDriver(DataEra, channel, run_permutation_tree, run_chi, rm_wm_constraint, rm_bjet_energy_reg_nn);

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

  if (!IsDATA && run_syst)
  {
    vec_syst_type = {AnalyzerParameter::Central,
                     AnalyzerParameter::JetEnDown,
                     AnalyzerParameter::JetEnUp,
                     AnalyzerParameter::JetResDown,
                     AnalyzerParameter::JetResUp,
                     AnalyzerParameter::UnclusteredEnergyDown,
                     AnalyzerParameter::UnclusteredEnergyUp};
    // if (run_el_ch)
    // {
    //   vec_syst_type.push_back(AnalyzerParameter::ElectronEnDown);
    //   vec_syst_type.push_back(AnalyzerParameter::ElectronEnUp);
    //   vec_syst_type.push_back(AnalyzerParameter::ElectronResDown);
    //   vec_syst_type.push_back(AnalyzerParameter::ElectronResUp);
    // }
  }
  else
    vec_syst_type = {AnalyzerParameter::Central};

  if (run_permutation_tree)
  {
    chk_matched_jets_only = false;

    Set_Permutation_Tree();
  } // if(run_permutation_tree)

  if (run_hf_contamination_tree)
  {
    chk_matched_jets_only = false;

    Set_HF_Contamination_Tree();
  }

  if (run_template_truth)
  {
    Set_Template_Truth_Tree();
  } // run_template_truth

  if (run_result)
  {
    chk_matched_jets_only = false;

    Set_Result_Tree();
  }

  if (!run_permutation_tree)
    Set_Reader_Swapper();

  // Set_Reader_HF_Contamination();

  return;
} // void Vcb::initializeAnalyzer()

//////////

void Vcb::executeEvent()
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
    // setup setting
    param.Clear();

    param.Muon_Tight_ID = "POGTightWithTightIso";
    param.Muon_Loose_ID = "POGLoose";

    param.Muon_ID_SF_Key = "NUM_TightID_DEN_TrackerMuons";
    param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";

    // param.Electron_Tight_ID = "passTightID";
    // param.Electron_Loose_ID = "passLooseID";

    param.Electron_Tight_ID = "passMVAID_iso_WP80";
    param.Electron_Loose_ID = "passMVAID_iso_WP90";

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
} // void Vcb::executeEvent()

//////////

void Vcb::executeEventFromParameter(AnalyzerParameter param)
{
  // cout << param.Name << endl;

  Clear();

  Event ev = GetEvent();

  if (!IsDATA)
  {
    // separate W decay mode
    vector<Gen> vec_gen = GetGens();
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

    // if(!TMath::Finite(lumi_weight)) cout << lumi_weight << endl;
    // if(!TMath::Finite(mc_weight)) cout << mc_weight << endl;
    // if(!TMath::Finite(pileup_weight)) cout << pileup_weight << endl;
    // if(!TMath::Finite(prefire_weight)) cout << prefire_weight << endl;
  }

  // no cut
  FillHist(param.Name + "/Cut_Flow", Cut_Flow::No_Cut, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::No_Cut, weight, n_cut_flow, 0, n_cut_flow);

  // met filter
  if (!PassMETFilter())
    return;

  FillHist(param.Name + "/Cut_Flow", Cut_Flow::Met_Filter, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Met_Filter, weight, n_cut_flow, 0, n_cut_flow);

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

  // xy correction
  pair<double, double> met_corr = xy_met_correction.METXYCorr_Met_MetPhi(met.Pt(), met.Phi(), run, to_string(DataYear), !IsDATA, nPV, true, true);

  // Particle met = ev.GetMETVector("PF");
  // pair<double, double> met_corr = xy_met_correction.METXYCorr_Met_MetPhi(met.Pt(), met.Phi(), run, to_string(DataYear), !IsDATA, nPV, true, false);
  // cout << met.Pt() << " " << met_corr.first << " " << met.Phi() << " " << met_corr.second << endl;

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
  vec_sel_jet = JetsVetoLeptonInside(vec_sel_jet, vec_electron_veto, vec_muon_veto, DR_LEPTON_VETO);
  n_sel_jet = vec_sel_jet.size();

  int n_pu_jet_before = 0;
  for (unsigned int i = 0; i < vec_sel_jet.size(); i++)
  {
    Jet jet = vec_sel_jet.at(i);

    int hf_flavour = jet.GenHFHadronMatcherFlavour();
    int hf_origin = jet.GenHFHadronMatcherOrigin();

    if (hf_flavour == -999 && hf_origin == -999)
      n_pu_jet_before++;
  }
  int n_real_jet_before = n_sel_jet - n_pu_jet_before;

  vec_sel_jet = SelectJets(vec_sel_jet, param.PUJet_Veto_ID, JET_PT, jet_eta_cut);
  n_sel_jet = vec_sel_jet.size();

  int n_pu_jet_after = 0;
  for (unsigned int i = 0; i < vec_sel_jet.size(); i++)
  {
    Jet jet = vec_sel_jet.at(i);

    int hf_flavour = jet.GenHFHadronMatcherFlavour();
    int hf_origin = jet.GenHFHadronMatcherOrigin();

    if (hf_flavour == -999 && hf_origin == -999)
      n_pu_jet_after++;
  }
  int n_real_jet_after = n_sel_jet - n_pu_jet_after;

  weight_hem_veto = Weight_HEM_Veto(vec_sel_jet);

  // jets for gen matching. Not to introduce matching bias, jets for matching have loose pt cut
  vec_sel_jet_match = SelectJets(vec_this_jet, param.Jet_ID, JET_PT_MATCH, JET_ETA_MATCH);
  vec_sel_jet_match = SelectJets(vec_sel_jet_match, param.PUJet_Veto_ID, JET_PT_MATCH, JET_ETA_MATCH);

  // sort jet as pt ordering
  sort(vec_sel_jet.begin(), vec_sel_jet.end(), PtComparing);
  sort(vec_sel_jet_match.begin(), vec_sel_jet_match.end(), PtComparing);

  for (unsigned int i = 0; i < vec_sel_jet.size(); i++)
  {
    Jet jet = vec_sel_jet[i];

    vec_sel_gen_hf_flavour.push_back(jet.GenHFHadronMatcherFlavour());
    vec_sel_gen_hf_origin.push_back(jet.GenHFHadronMatcherOrigin());
  }

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

  if (!IsDATA)
  {
    if (run_mu_ch)
    {
      // SF for muon trigger effi
      weight_sl_trig = mcCorr->MuonTrigger_SF("POGTight", sl_trig, vec_sel_muon, 0);
      if (run_syst)
      {
        weight_sl_trig_down = mcCorr->MuonTrigger_SF("POGTight", sl_trig, vec_sel_muon, -1);
        weight_sl_trig_up = mcCorr->MuonTrigger_SF("POGTight", sl_trig, vec_sel_muon, +1);
      }
    }
    else if (run_el_ch)
    {
      weight_sl_trig = mcCorr->ElectronTrigger_SF(param.Electron_Tight_ID, sl_trig, vec_sel_electron, 0);
      if (run_syst)
      {
        weight_sl_trig_down = mcCorr->ElectronTrigger_SF(param.Electron_Tight_ID, sl_trig, vec_sel_electron, -1);
        weight_sl_trig_up = mcCorr->ElectronTrigger_SF(param.Electron_Tight_ID, sl_trig, vec_sel_electron, +1);
      }
    }
    weight *= weight_sl_trig;

    // cout << "n_muon = " << vec_sel_muon.size() << endl;
    // for(int i=0; i<vec_sel_muon.size(); i++)
    //   {
    //     Muon muon = vec_sel_muon.at(i);
    //     cout << muon.Eta() << "\t" << muon.MiniAODPt() << endl;
    //   }
    // cout << sf_sl_trig_effi << endl;
  }
  FillHist(param.Name + "/Cut_Flow", Cut_Flow::Trigger, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Trigger, weight, n_cut_flow, 0, n_cut_flow);

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

  if (!IsDATA)
  {
    if (run_mu_ch)
    {
      // SF for muon id
      weight_mu_id = mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muon.Eta(), muon.MiniAODPt(), 0);
      if (run_syst)
      {
        weight_mu_id_down = mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muon.Eta(), muon.MiniAODPt(), -1);
        weight_mu_id_up = mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muon.Eta(), muon.MiniAODPt(), +1);
      }
      weight *= weight_mu_id;

      // SF for muon iso
      weight_mu_iso = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muon.Eta(), muon.MiniAODPt(), 0);
      if (run_syst)
      {
        weight_mu_iso_down = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muon.Eta(), muon.MiniAODPt(), -1);
        weight_mu_iso_up = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muon.Eta(), muon.MiniAODPt(), +1);
      }
      weight *= weight_mu_iso;

      // if(run_debug) cout << "M ISO SF = " << muon.Eta() << "\t" << muon.MiniAODPt() << "\t" << sf_mu_iso_effi << endl;
    }
    else if (run_el_ch)
    {
      // SF for electron id
      weight_el_id = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electron.scEta(), electron.UncorrPt(), 0);
      if (run_syst)
      {
        weight_el_id_down = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electron.scEta(), electron.UncorrPt(), -1);
        weight_el_id_up = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electron.scEta(), electron.UncorrPt(), +1);
      }
      weight *= weight_el_id;

      // SF for electron Reco eff
      weight_el_reco = mcCorr->ElectronReco_SF(electron.scEta(), electron.UncorrPt(), 0);
      if (run_syst)
      {
        weight_el_reco_down = mcCorr->ElectronReco_SF(electron.scEta(), electron.UncorrPt(), -1);
        weight_el_reco_up = mcCorr->ElectronReco_SF(electron.scEta(), electron.UncorrPt(), +1);
      }
      weight *= weight_el_reco;
    }
  }
  FillHist(param.Name + "/Cut_Flow", Cut_Flow::Single_Lepton, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::Single_Lepton, weight, n_cut_flow, 0, n_cut_flow);

  // cut on jet
  if (n_sel_jet < 4)
    return;

  if (!IsDATA)
  {
    // SF for PUJet Veto
    weight_pujet_veto = mcCorr->PileupJetVeto_Reweight(vec_sel_jet, param.PUJet_Veto_ID, 0);
    if (run_syst)
    {
      weight_pujet_veto_down = mcCorr->PileupJetVeto_Reweight(vec_sel_jet, param.PUJet_Veto_ID, -1);
      weight_pujet_veto_up = mcCorr->PileupJetVeto_Reweight(vec_sel_jet, param.PUJet_Veto_ID, +1);
    }

    // weight *= weight_pujet_veto;
  }

  FillHist(param.Name + "/Cut_Flow", Cut_Flow::At_Least_Four_Jets, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::At_Least_Four_Jets, weight, n_cut_flow, 0, n_cut_flow);
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

  for (auto &jet : vec_sel_jet_match)
  {
    float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
    if (mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score)
      vec_btag_match.push_back(true);
    else
      vec_btag_match.push_back(false);
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

  if (!IsDATA)
  {
    // SF for b-tagging
    if (run_debug)
      weight_b_tag = mcCorr->GetBTaggingReweight_1a(vec_sel_jet, vec_jet_tagging_para.at(0));
    else
    {
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
    } //   if (run_debug)

    weight *= weight_b_tag;
  }

  if (!IsDATA)
  {
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
  }

  if (n_b_jet < 2)
    return;

  FillHist(param.Name + "/Cut_Flow", Cut_Flow::At_Least_Two_B_Tagged, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::At_Least_Two_B_Tagged, weight, n_cut_flow, 0, n_cut_flow);

  // cut on MET
  met_pt = met.Pt();
  met_phi = met.Phi();

  // if (met_pt < MET_PT)
  //   return;

  pt_ratio = lepton_pt / met_pt;

  FillHist(param.Name + "/Cut_Flow", Cut_Flow::MET, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::MET, weight, n_cut_flow, 0, n_cut_flow);

  FillHist(param.Name + "/N_Real_Jets_Before_PUJet_Veto", n_real_jet_before, weight, 20, 0, 20);
  FillHist(param.Name + "/N_PU_Jets_Before_PUJet_Veto", n_pu_jet_before, weight, 20, 0, 20);
  FillHist(param.Name + "/N_Real_Jets_After_PUJet_Veto", n_real_jet_after, weight, 20, 0, 20);
  FillHist(param.Name + "/N_PU_Jets_After_PUJet_Veto", n_pu_jet_after, weight, 20, 0, 20);

  // Gen matching
  vector<int> vec_hf_flavour;
  vector<int> vec_hf_origin;
  int index_gen[4];
  bool surely_matched[4];
  float matched_jet_dr[4];
  int switch_included = -1;
  if (!IsDATA)
  {
    // scan GenHFHadronMatcher results and construct JER for gen match
    vector<float> vec_jer_match;
    for (auto &jet : vec_sel_jet_match)
    {
      int jet_flavour = jet.GenHFHadronMatcherFlavour();
      int jet_origin = jet.GenHFHadronMatcherOrigin();

      vec_hf_flavour.push_back(jet_flavour);
      vec_hf_origin.push_back(jet_origin);

      float resolution_pt = jet_resolution.getResolution({{JME::Binning::JetPt, jet.Pt()}, {JME::Binning::JetEta, jet.Eta()}, {JME::Binning::Rho, Rho}});
      float resolution_pt_sf = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jet.Pt()}, {JME::Binning::JetEta, jet.Eta()}}, Variation::NOMINAL);

      vec_jer_match.push_back(resolution_pt * resolution_pt_sf);
    }

    // try to find jets from W using closed delta R matching
    vector<Gen> vec_gen = GetGens();
    // PrintGen(vec_gen);

    Gen_Match_Lepton(lepton, vec_gen, chk_gentau_conta);
    Gen_Match_TT(vec_sel_jet_match, vec_gen, vec_hf_flavour, vec_hf_origin, vec_jer_match, index_gen, index_matched_jet_match, surely_matched, matched_jet_dr);

    bool chk_gen_acceptance = true;
    for (int i = 0; i < 4; ++i)
    {
      int index = index_gen[i];

      float pt = vec_gen[index].Pt();
      float eta = vec_gen[index].Eta();

      if (pt < JET_PT_MATCH || JET_ETA_MATCH < TMath::Abs(eta))
      {
        chk_gen_acceptance = false;
        break;
      }
    }

    if (chk_gen_acceptance)
      FillHist(param.Name + "/Gen_Acceptance", 1, weight, 2, 0, 2);
    else
      FillHist(param.Name + "/Gen_Acceptance", 0, weight, 2, 0, 2);

    // cout << "test index_matched_jet_match[i]" << endl;
    // for(int i=0; i<4; i++){ cout << index_matched_jet_match[i] << " " << matched_jet_dr[i] << endl; }

    // if all four jets are successfully match
    if (-1 < index_matched_jet_match[0] && -1 < index_matched_jet_match[1] && -1 < index_matched_jet_match[2] && -1 < index_matched_jet_match[3])
    {
      FillHist(param.Name + "/PF_Gen_Matched", 1, weight, 2, 0, 2);

      // BJetRegression Test
      TLorentzVector w_gen_matched = vec_sel_jet_match.at(index_matched_jet_match[1]) + vec_sel_jet_match.at(index_matched_jet_match[2]);
      TLorentzVector t_gen_matched = vec_sel_jet_match.at(index_matched_jet_match[0]) + w_gen_matched;

      FillHist(param.Name + "/W_Gen_Matched_Mass", w_gen_matched.M(), weight, 100, 0, 400);
      FillHist(param.Name + "/T_Gen_Matched_Mass", t_gen_matched.M(), weight, 100, 0, 600);

      // correction
      Jet j0 = vec_sel_jet_match.at(index_matched_jet_match[0]);
      Jet j1 = vec_sel_jet_match.at(index_matched_jet_match[1]);
      Jet j2 = vec_sel_jet_match.at(index_matched_jet_match[2]);

      float j0_corr = j0.BJetNNCorrection();
      j0.SetPtEtaPhiM(j0_corr * j0.Pt(), j0.Eta(), j0.Phi(), j0.M());

      float j1_corr = 1;
      if (vec_hf_flavour[index_matched_jet_match[1]] == 4)
        j1_corr = j1.CJetNNCorrection();
      else if (vec_hf_flavour[index_matched_jet_match[1]] == 5)
        j1_corr = j1.BJetNNCorrection();

      j1.SetPtEtaPhiM(j1_corr * j1.Pt(), j1.Eta(), j1.Phi(), j1.M());

      float j2_corr = 1;
      if (vec_hf_flavour[index_matched_jet_match[2]] == 4)
        j2_corr = j2.CJetNNCorrection();
      else if (vec_hf_flavour[index_matched_jet_match[2]] == 5)
        j2_corr = j2.BJetNNCorrection();

      j2.SetPtEtaPhiM(j2_corr * j2.Pt(), j2.Eta(), j2.Phi(), j2.M());

      TLorentzVector w_gen_matched_corr = j1 + j2;
      TLorentzVector t_gen_matched_corr = j0 + w_gen_matched_corr;

      FillHist(param.Name + "/W_Gen_Matched_Mass_CorrNN", w_gen_matched_corr.M(), weight, 100, 0, 400);
      FillHist(param.Name + "/T_Gen_Matched_Mass_CorrNN", t_gen_matched_corr.M(), weight, 100, 0, 600);

      // Delta R
      for (int i = 0; i < 4; i++)
      {
        if (surely_matched[i] == true)
          FillHist(param.Name + Form("/DR_Surely_Matched_%d", i), matched_jet_dr[i], weight, 30, 0, 3);
      }

      if (matched_jet_dr[1] < matched_jet_dr[2])
      {
        FillHist(param.Name + "/W_DR_Small", matched_jet_dr[1], weight, 30, 0, 3);
        FillHist(param.Name + "/W_DR_Large", matched_jet_dr[2], weight, 30, 0, 3);
      }
      else
      {
        FillHist(param.Name + "/W_DR_Small", matched_jet_dr[2], weight, 30, 0, 3);
        FillHist(param.Name + "/W_DR_Large", matched_jet_dr[1], weight, 30, 0, 3);
      }

      // if all matched four jets from tt system are included in vec_jet_sel, chk_included=true
      switch_included = Chk_Included(index_matched_jet);
    }
    // at least of matching failed
    else
    {
      FillHist(param.Name + "/PF_Gen_Matched", 0, weight, 2, 0, 2);
    }

    //-1 gen matching failed
    // 0 gen matching succeeded && four jets are in selection
    // 0< gen matching succeeded &&
    // first bjet is out of selection -> 1 bit
    // first wjet is out of selection -> 2 bit
    // second wjet is out of selection -> 3 bit
    // second bjet is out of selection -> 4 bit
    FillHist(param.Name + "/Objects_In_Sel_Jet", switch_included, weight, 18, -2, 16);

    // To match index of jet_match (loose jet wo/ lepton veto) and jet. If
    Index_Converter(vec_sel_jet, vec_sel_jet_match, index_matched_jet_match, index_matched_jet);

    // cout << "test index_matched_jet[i]" << endl;
    // for(int i=0; i<4; i++){ cout << index_matched_jet[i] << endl; }

    // To estimate the fraction of the included events
    if (Chk_Included(index_matched_jet) == 0)
      FillHist(param.Name + "/Included", 1, weight, 2, 0, 2);
    else
      FillHist(param.Name + "/Included", 0, weight, 2, 0, 2);
  } // if(!IsDATA)

  // run_Template_truth
  if (!IsDATA && run_template_truth)
  {
    Make_Template_Truth_Tree();
    return;
  } // if(!IsDATA && run_template_truth)

  // kinematic fitter
  vector<float> vec_resolution_pt;
  for (auto &jet : vec_sel_jet)
  {
    // JEC
    float resolution_pt = jet_resolution.getResolution({{JME::Binning::JetPt, jet.Pt()}, {JME::Binning::JetEta, jet.Eta()}, {JME::Binning::Rho, Rho}});
    float resolution_pt_sf = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jet.Pt()}, {JME::Binning::JetEta, jet.Eta()}}, Variation::NOMINAL);

    vec_resolution_pt.push_back(resolution_pt * resolution_pt_sf);
  }

  fitter_driver->Set_Objects(vec_sel_jet, vec_resolution_pt, vec_btag, lepton, met, chk_matched_jets_only, index_matched_jet);
  fitter_driver->Scan();

  if (!run_permutation_tree && !run_chi)
    KF_Ambiguity_Remover(vec_sel_jet, index_matched_jet);

  // Permutation Tree for signal study
  if (!IsDATA && run_permutation_tree)
  {
    Make_Permutation_Tree();

    return;
  }

  // HF Contamination
  if (!IsDATA && run_hf_contamination_tree)
  {
    Make_HF_Contamination_Tree();

    return;
  } // if(!IsDATA && run_hf_contamination_tree)

  // run_template
  if (!IsDATA && run_template)
  {
    Make_Template_Tree();

    return;
  } // if(!IsDATA && run_template)

  if (run_result)
  {
    Make_Result_Tree(param);

    return;
  } // if(run_result)

  return;
} // Vcb::executeEventFromParameter(AnalyzerParameter param)

//////////

float Vcb::Calculate_HT(const vector<Jet> &vec_jet)
{
  float ht = 0;

  for (unsigned int i = 0; i < vec_jet.size(); i++)
    ht += vec_jet[i].Pt();

  return ht;
} // float Vcb::Calculate_HT(const vector<Jet>& vec_sel)

//////////

float Vcb::Calculate_Mt(const Particle &lepton, const float &neu_px, const float &neu_py)
{
  // float px = lepton.Px() + neu_px;
  // float py = lepton.Py() + neu_py;

  // float mt = TMath::Sqrt(W_MASS*W_MASS + px*px + py*py);

  float lepton_px = lepton.Px();
  float lepton_py = lepton.Py();
  float lepton_et = Sqrt(lepton_px * lepton_px + lepton_py * lepton_py);

  float neu_et = Sqrt(neu_px * neu_px + neu_py * neu_py);

  TVector3 lepton3(lepton_px, lepton_py, 0);
  TVector3 neu3(neu_px, neu_py, 0);

  float angle = lepton3.Angle(neu3);

  float mt = Sqrt(lepton_et * neu_et * (1 - Cos(angle)));

  return mt;
} // float Vcb::Calculate_Mt(const Particle& lepton, const Particle& met)

//////////

void Vcb::Clear()
{
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

  weight_hem_veto = 1;
  weight_lumi = 1;
  weight_mc = 1;

  weight_pdf_alternative = 1;
  memset(weight_pdf_error_set, 1, sizeof(weight_pdf_error_set));
  weight_pdf_as_up = 1;
  weight_pdf_as_down = 1;

  weight_pileup = 1;
  weight_pileup_down = 1;
  weight_pileup_up = 1;

  weight_pujet_veto = 1;
  weight_pujet_veto_down = 1;
  weight_pujet_veto_up = 1;

  weight_prefire = 1;
  weight_prefire_down = 1;
  weight_prefire_up = 1;

  weight_scale_variation_1 = 1;
  weight_scale_variation_2 = 1;
  weight_scale_variation_3 = 1;
  weight_scale_variation_4 = 1;
  weight_scale_variation_6 = 1;
  weight_scale_variation_8 = 1;

  weight_top_pt = 1;

  weight_mu_id = 1;
  weight_mu_id_down = 1;
  weight_mu_id_up = 1;

  weight_mu_iso = 1;
  weight_mu_iso_down = 1;
  weight_mu_iso_up = 1;

  weight_el_id = 1;
  weight_el_id_down = 1;
  weight_el_id_up = 1;

  weight_el_reco = 1;
  weight_el_reco_down = 1;
  weight_el_reco_up = 1;

  weight_sl_trig = 1;
  weight_sl_trig_down = 1;
  weight_sl_trig_up = 1;

  decay_mode = 999;

  chk_gentau_conta = false;

  vec_this_muon.clear();
  vec_this_electron.clear();
  vec_this_jet.clear();

  vec_sel_jet.clear();
  vec_sel_jet_match.clear();

  vec_btag.clear();
  vec_btag_match.clear();
  vec_ctag.clear();

  memset(index_matched_jet, -999, sizeof(index_matched_jet));
  memset(index_matched_jet_match, -999, sizeof(index_matched_jet_match));

  return;
} // void Vcb::Clear()

//////////

int Vcb::Chk_Included(const int index_matched_jet[4])
{
  int result = 0;
  unsigned int tmp = 1;
  for (int i = 0; i < 4; ++i)
  {
    if (index_matched_jet[i] < 0)
      result += tmp;

    tmp = tmp << 1;
  } // for(int i=0; i<4; i++)

  // 0: every gen jet is matched to loose jets wo/ lepton veto and the loose jets passes baseline selection
  return result;
} // int Vcb::Chk_Included(const int index_matched_jet[4])

//////////

bool Vcb::Compare_Jet(const Jet &jet0, const Jet &jet1)
{
  float pt0 = jet0.Pt();
  float pt1 = jet1.Pt();

  if (Abs(pt0 - pt1) < 1e-8)
    return true;
  else
    return false;
} // bool Vcb::Compare_Jet(const Jet& jet0, const Jet& jet1)

//////////

int Vcb::Compare_Jet_Pair(const Jet jet0[2], const Jet jet1[2])
{
  Double_t pt0[2] = {jet0[0].Pt(), jet0[1].Pt()};
  // float eta0[2] = {jet0[0].Eta(), jet0[1].Eta()};
  // float phi0[2] = {jet0[0].Phi(), jet0[1].Phi()};

  Double_t pt1[2] = {jet1[0].Pt(), jet1[1].Pt()};
  // float eta1[2] = {jet1[0].Eta(), jet1[1].Eta()};
  // float phi1[2] = {jet1[0].Phi(), jet1[1].Phi()};

  // Same jet pair
  if (Abs(pt0[0] - pt1[0]) < 1e-8 && Abs(pt0[1] - pt1[1]) < 1e-8)
    return 1;

  // Same jet pair but, swapped_truth
  else if (Abs(pt0[0] - pt1[1]) < 1e-8 && Abs(pt0[1] - pt1[0]) < 1e-8)
    return 2;

  // different jet pair
  return 0;
} // bool Vcb::Compare_Jet_Pair(const Jet& jet0, const Jet& jet1)

//////////

bool Vcb::Gen_Match_Lepton(const Lepton &lepton, const vector<Gen> &vec_gen, bool &chk_gentau_conta)
{
  // PrintGen(vec_gen);

  const int w_pdg_id = 24;
  const int aw_pdg_id = -24;

  // find leptonically decaying W
  int index_last_t = -999;
  int index_last_at = -999;
  int index_first_w = -999;
  int index_first_aw = -999;
  int index_last_w = -999;
  int index_last_aw = -999;
  int index_first_lepton = -999;
  int index_last_lepton = -999;
  int lepton_pid = -999;
  for (unsigned int i = 0; i < vec_gen.size(); i++)
  {
    Gen gen = vec_gen.at(i);

    int pid = gen.PID();
    int m_index = gen.MotherIndex();

    // find last index of t and tbar
    if (pid == 6)
      index_last_t = i;
    if (pid == -6)
      index_last_at = i;

    // find w(aw) from t decay
    if (m_index == index_last_t && pid == w_pdg_id)
    {
      index_first_w = i;
      index_last_w = index_first_w;
    }
    if (m_index == index_last_at && pid == aw_pdg_id)
    {
      index_first_aw = i;
      index_last_aw = index_first_aw;
    }

    // find last index of w(aw) from t(tbar) decay
    if (m_index == index_last_w && pid == w_pdg_id)
    {
      index_last_w = i;
    }
    if (m_index == index_last_aw && pid == aw_pdg_id)
    {
      index_last_aw = i;
    }

    // find first index of lepton from W decay
    if (m_index == index_last_w && (pid == -11 || pid == -13 || pid == -15))
    {
      index_first_lepton = i;
      index_last_lepton = index_first_lepton;
      lepton_pid = pid;
    }
    if (m_index == index_last_aw && (pid == 11 || pid == 13 || pid == 15))
    {
      index_first_lepton = i;
      index_last_lepton = index_first_lepton;
      lepton_pid = pid;
    }

    // find last index of lepton from W decay
    if (m_index == index_last_lepton && pid == lepton_pid)
      index_last_lepton = i;
  }

  // GenTau contamination
  if (abs(lepton_pid) == 15)
    chk_gentau_conta = true;

  // // dr matching
  // Gen gen_lepton = vec_gen.at(index_last_lepton);

  // float gen_lepton_pt = gen_lepton.Pt();
  // float gen_lepton_eta = gen_lepton.Eta();
  // float gen_lepton_phi = gen_lepton.Phi();

  // float lepton_pt = lepton.Pt();
  // float lepton_eta = lepton.Eta();
  // float lepton_phi = lepton.Phi();

  // float dr = Sqrt(Power(gen_lepton_eta - lepton_eta, 2.) + Power(gen_lepton_phi - lepton_phi, 2.));

  // cout << "index_first_w " << index_first_w << ", index_first_aw = " << index_first_aw << endl;
  // cout << "index_last_w " << index_last_w << ", index_last_aw = " << index_last_aw << endl;
  // cout << "index_first_lepton = " << index_first_lepton << ", index_last_lepton = " << index_last_lepton << ", lepton_pid = " << lepton_pid << endl;
  // cout << "gen_lepton_pt = " << gen_lepton_pt << ", gen_lepton_eta = " << gen_lepton_eta << ", gen_lepton_phi = " << gen_lepton_phi << endl;
  // cout << "lepton_pt = " << lepton_pt << ", lepton_eta = " << lepton_eta << ", lepton_phi = " << lepton_phi << endl;
  // cout << "dr = " << dr << endl;

  return false;
} // bool Vcb::Gen_Match_Lepton(const Lepton &lepton, const vector<Gen> &vec_gen)

//////////

void Vcb::Gen_Match_Residual(const vector<Jet> &vec_jet, const vector<Gen> &vec_gen, const vector<int> &vec_hf_flavour, const vector<int> &vec_hf_origin, const vector<float> &vec_jer, int index_gen[4], int index_matched_jet[4], bool surely_matched[4], float dr_return[4])
{
  // For gen particles which are not matched surely, this method will try to match those to jeco jet using DR matching

  const float rel_pt_diff_cut = 3;
  const float dr_cut = 0.4;

  int possible_matched_index[4] = {-999, -999, -999, -999};
  float dr_smallest[4] = {999, 999, 999, 999};
  for (int i = 0; i < 4; ++i)
  {
    // cout << "gen " << i << " " << index_gen[i] << " " << index_matched_jet[i] << endl;

    if (index_matched_jet[i] != -999)
      continue;

    Gen gen = vec_gen[index_gen[i]];

    float gen_pt = gen.Pt();
    float gen_eta = gen.Eta();
    float gen_phi = gen.Phi();

    // cout << "gen pt = " << gen_pt << " " << gen_eta << " " << gen_phi << endl;

    // acceptance test. If gen particle is out of acceptance, matching is not tried
    //  if(JET_ETA_MATCH<Abs(gen_eta))
    //   	{
    //  	  //cout << "Out of acceptance" << endl;
    //  	  index_matched_jet[i] = -1;//gen is out of acceptance
    //   	  continue;
    //   	}

    for (unsigned int j = 0; j < vec_jet.size(); ++j)
    {
      // if jet is aleady matched surely with other gen particles, continue
      bool chk_allocated = false;
      for (int k = 0; k < 4; ++k)
      {
        if ((signed)j == index_matched_jet[k])
        {
          chk_allocated = true;
          break;
        }
      }
      if (chk_allocated)
        continue;

      // if origin is tagged by GenHFHadronMatcher, continue
      // int jet_flavour = vec_hf_flavour.at(j);
      int jet_origin = vec_hf_origin.at(j);

      if (jet_origin != -999)
        continue;

      Jet jet = vec_jet.at(j);

      float jet_pt = jet.Pt();
      float jet_eta = jet.Eta();
      float jet_phi = jet.Phi();

      float del_pt = gen_pt - jet_pt;
      float del_eta = gen_eta - jet_eta;
      float del_phi = gen_phi - jet_phi;

      float rel_pt_diff = Abs(del_pt / vec_jer[i] / jet_pt);

      float dr = Sqrt(del_eta * del_eta + del_phi * del_phi);

      if (rel_pt_diff < rel_pt_diff_cut && dr < dr_smallest[i])
      {
        dr_smallest[i] = dr;
        possible_matched_index[i] = j;
      }

      // cout << j << " " << jet_flavour << " " << jet_origin << " " << rel_pt_diff << " " << dr << " " << possible_matched_index[i] << endl;

    } // loop over jet

    // if no jet is allocated to gen particle
    if (possible_matched_index[i] == -999 || dr_cut < dr_smallest[i])
    {
      // cout << "impossible to allocate jet to gen " << i << endl;
      index_matched_jet[i] = -2; // impossible to allocate any jet to gen particle
    }

    // cout << endl;
  } // loop over gen

  // if two or more gen particles are matched to same jet, gen particle with smallest dr is chosen
  // then recursively call this method to allcate jet to the second ranked gen
  for (int i = 0; i < 3; ++i)
  {
    if (index_matched_jet[i] != -999)
      continue;

    bool chk_overlap[4] = {false, false, false, false};
    chk_overlap[i] = true;

    for (int j = i + 1; j < 4; ++j)
    {
      if (index_matched_jet[j] != -999)
        continue;

      // same jet allocation found
      if (possible_matched_index[i] == possible_matched_index[j])
        chk_overlap[j] = true;
    }

    // find gen particle with smallest dr among overlaped gen
    int gen_index_best_among_overlap = -999;
    float dr_smallest_among_overlap = 999;
    for (int j = i; j < 4; ++j)
    {
      if (chk_overlap[j] != true)
        continue;

      if (dr_smallest[j] < dr_smallest_among_overlap)
      {
        gen_index_best_among_overlap = j;
        dr_smallest_among_overlap = dr_smallest[j];
      }
    }

    index_matched_jet[gen_index_best_among_overlap] = possible_matched_index[gen_index_best_among_overlap];
    dr_return[gen_index_best_among_overlap] = dr_smallest[gen_index_best_among_overlap];

    Gen_Match_Residual(vec_jet, vec_gen, vec_hf_flavour, vec_hf_origin, vec_jer, index_gen, index_matched_jet, surely_matched, dr_return);
  }

  return;
} // void Vcb::Gen_Match_Residual(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, int index_gen[4], int index_matched_jet[4], bool surely_matched[4])

//////////

void Vcb::Gen_Match_TT(const vector<Jet> &vec_jet, const vector<Gen> &vec_gen, const vector<int> &vec_hf_flavour, const vector<int> &vec_hf_origin, const vector<float> &vec_jer, int index_gen[4], int index_matched_jet[4], bool surely_matched[4], float dr_return[4])
{
  for (int i = 0; i < 4; i++)
  {
    index_gen[i] = -999;
    index_matched_jet[i] = -999;
    surely_matched[i] = false;
    dr_return[i] = -999;
  }

  int selected_w = Gen_Match_W(vec_jet, vec_gen, vec_hf_flavour, vec_hf_origin, vec_jer, &index_gen[1], &index_matched_jet[1], &surely_matched[1], &dr_return[1]);

  int index_last_t = -999;
  int index_last_at = -999;
  int index_first_b = -999;
  int index_first_ab = -999;
  for (unsigned int i = 0; i < vec_gen.size(); i++)
  {
    Gen gen = vec_gen.at(i);

    int pid = gen.PID();
    int m_index = gen.MotherIndex();

    // find last index of t and tbar
    if (pid == 6)
      index_last_t = i;
    if (pid == -6)
      index_last_at = i;

    // find b from t decay
    if (m_index == index_last_t && pid == 5)
      index_first_b = i;
    if (m_index == index_last_at && pid == -5)
      index_first_ab = i;
  }

  // W+ decays hadronically
  if (selected_w == 1)
  {
    index_gen[0] = index_first_b;
    index_gen[3] = index_first_ab;
  }
  // W- decays hadronically
  else if (selected_w == -1)
  {
    index_gen[0] = index_first_ab;
    index_gen[3] = index_first_b;
  }
  // no W decays hadronically
  else
  {
    if (run_debug)
      cout << "Can't find hadronically decayed W" << endl;

    return;
  }

  for (unsigned int i = 0; i < vec_jet.size(); i++)
  {
    int jet_flavour = vec_hf_flavour.at(i);
    int jet_origin = vec_hf_origin.at(i);

    float jet_eta = vec_jet.at(i).Eta();
    float jet_phi = vec_jet.at(i).Phi();

    // b from t
    if (jet_flavour == 5 && jet_origin == 6 && index_first_b != -999)
    {
      Gen gen_b = vec_gen.at(index_first_b);

      float del_eta = jet_eta - gen_b.Eta();
      float del_phi = jet_phi - gen_b.Phi();

      float dr = Sqrt(Power(del_eta, 2.0) + Power(del_phi, 2.0));

      if (selected_w == 1)
      {
        index_matched_jet[0] = i;
        surely_matched[0] = true;
        dr_return[0] = dr;
      }
      else
      {
        index_matched_jet[3] = i;
        surely_matched[3] = true;
        dr_return[3] = dr;
      }
    }

    // ab from at
    if (jet_flavour == 5 && jet_origin == -6 && index_first_ab != -999)
    {
      Gen gen_ab = vec_gen.at(index_first_ab);

      float del_eta = jet_eta - gen_ab.Eta();
      float del_phi = jet_phi - gen_ab.Phi();

      float dr = Sqrt(Power(del_eta, 2.0) + Power(del_phi, 2.0));

      if (selected_w == 1)
      {
        index_matched_jet[3] = i;
        surely_matched[3] = true;
        dr_return[3] = dr;
      }
      else
      {
        index_matched_jet[0] = i;
        surely_matched[0] = true;
        dr_return[0] = dr;
      }
    }
  } // for(unsigned int i=0; i<vec_jet.size(); i++)

  // if(run_debug)
  //   {
  //     cout << endl;
  //     cout << endl;
  //     if(selected_w == 1) cout << "W+ decayed hadronically" << endl;
  //     else cout << "W- decayed hadronically" << endl;

  //     for(unsigned int i=0; i<4; i++)
  //     	{
  //     	  if(index_gen[i]!=-999)
  //     	    {
  //     	      Gen gen = vec_gen.at(index_gen[i]);
  //     	      int pid = gen.PID();

  //     	      Gen gen_mother = vec_gen.at(gen.MotherIndex());
  //     	      int m_pid = gen_mother.PID();

  //     	      float gen_eta = gen.Eta();
  //     	      float  gen_phi = gen.Phi();
  //     	      float gen_pt = gen.Pt();

  //     	      cout << "(" << index_gen[i] << ", " << pid << ", " << m_pid << ", " << gen_pt << ", " << gen_eta << ", " << gen_phi << "), ";

  //     	    }
  //     	  else
  //     	    {
  //     	      cout << "(), ";
  //     	    }
  //     	}
  //     cout << endl;

  //     for(unsigned int i=0; i<vec_jet.size(); i++)
  //     	{
  //     	  int jet_flavour = vec_hf_flavour.at(i);
  //         int jet_origin = vec_hf_origin.at(i);

  //     	  Jet jet = vec_jet.at(i);

  //     	  float jet_eta = jet.Eta();
  //     	  float jet_phi = jet.Phi();
  //     	  float jet_pt =  jet.Pt();

  //     	  float tagging_score = jet.GetTaggerResult(JetTagging::DeepJet);
  //     	  bool b_tagged = false;
  //     	  if(mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < tagging_score) b_tagged = true;

  //     	  cout << "(" << jet_flavour << ", " << jet_origin << ", " << jet_pt << ", " << jet_eta << ", " << jet_phi << ", " << ", " << b_tagged << ")" << endl;
  //     	}
  //     cout << endl;

  //     for(int i=0; i<4; i++)
  // 	{
  // 	  cout << "Matched jet index = " << index_matched_jet[i] << ", surely matched = " << surely_matched[i] << ", d_r = " << dr_return[i] << endl;
  // 	}
  //     cout << endl;
  //   }

  Gen_Match_Residual(vec_jet, vec_gen, vec_hf_flavour, vec_hf_origin, vec_jer, index_gen, index_matched_jet, surely_matched, dr_return);

  return;
} // void Vcb::Gen_Match_TT(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour, const vector<int>& vec_hf_origin, int index_matched_jet[4], bool surely_matched[4], float dr_return[4])

//////////

int Vcb::Gen_Match_W(const vector<Jet> &vec_jet, const vector<Gen> &vec_gen, const vector<int> &vec_hf_flavour, const vector<int> &vec_hf_origin, const vector<float> &vec_jer, int index_gen[2], int index_matched_jet[2], bool surely_matched[2], float dr_return[2])
{
  // this method is designed to find gen particles from W and to try to match to reco jet via GenHFHadronMatcher results

  const int target_pdg_id_positive = 24;
  const int target_pdg_id_negative = -24;

  // For CHtoCB decay
  // const int target_pdg_id_positive = 37;
  // const int target_pdg_id_negative = -37;

  int index_last_w = -999;
  int index_last_aw = -999;
  int index_d0_w = -999;
  int index_d1_w = -999;
  int index_d0_aw = -999;
  int index_d1_aw = -999;

  int selected_w = -999;

  // scan gen to find W
  for (unsigned int i = 0; i < vec_gen.size(); i++)
  {
    Gen gen = vec_gen.at(i);

    int pid = gen.PID();
    int m_index = gen.MotherIndex();

    // find last index of W+ and W-
    if (pid == target_pdg_id_positive)
      index_last_w = i;
    if (pid == target_pdg_id_negative)
      index_last_aw = i;

    // find decay products of W+
    if (m_index == index_last_w && pid != target_pdg_id_positive)
    {
      if (index_d0_w == -999)
        index_d0_w = i;
      else
        index_d1_w = i;
    }

    // find decay products of W-
    if (m_index == index_last_aw && pid != target_pdg_id_negative)
    {
      if (index_d0_aw == -999)
        index_d0_aw = i;
      else
        index_d1_aw = i;
    }
  }

  // if input sample is not TT, so both of W couldn't be found
  if (index_last_w == -999 || index_last_aw == -999)
  {
    // if(run_debug) cout << "Can't find both of W" << endl;
    return -999;
  }

  // W+ decay hadronically
  if (Abs(vec_gen.at(index_d0_w).PID()) < 10)
  {
    if (vec_gen.at(index_d0_w).PID() % 2 == 0)
    {
      index_gen[0] = index_d0_w;
      index_gen[1] = index_d1_w;
    }
    else
    {
      index_gen[0] = index_d1_w;
      index_gen[1] = index_d0_w;
    }

    selected_w = 1; // mean W+ decayed hadronically
  }
  // W- decay hadronically
  else if (Abs(vec_gen.at(index_d0_aw).PID() < 10))
  {
    if (vec_gen.at(index_d0_aw).PID() % 2 == 0)
    {
      index_gen[0] = index_d0_aw;
      index_gen[1] = index_d1_aw;
    }
    else
    {
      index_gen[0] = index_d1_aw;
      index_gen[1] = index_d0_aw;
    }

    selected_w = -1; // mean W- decayed hadronically
  }
  // no W decays hadronically
  else
  {
    if (run_debug)
      cout << "No W decays hadronically" << endl;

    return -999;
  }

  // matching with GenHFHadronMatcher information
  for (unsigned int i = 0; i < vec_jet.size(); i++)
  {
    int hf_flavour = vec_hf_flavour.at(i);
    int hf_origin = vec_hf_origin.at(i);

    Jet jet = vec_jet.at(i);

    float jet_eta = jet.Eta();
    float jet_phi = jet.Phi();

    // w is tagged by GenHFHadronMatcher
    if (Abs(hf_origin) == target_pdg_id_positive)
    {
      for (int j = 0; j < 2; j++)
      {
        if (Abs(hf_flavour) == Abs(vec_gen.at(index_gen[j]).PID()))
        {
          index_matched_jet[j] = i;
          surely_matched[j] = true;

          Gen gen = vec_gen.at(index_gen[j]);

          float gen_eta = gen.Eta();
          float gen_phi = gen.Phi();

          float del_eta = jet_eta - gen_eta;
          float del_phi = jet_phi - gen_phi;

          float del_r = Sqrt(Power(del_eta, 2.) + Power(del_phi, 2.));

          dr_return[j] = del_r;
        }
      }
    }
  }

  // cout << endl;
  // cout << "Daughters from W";
  // if(selected_w>0) cout << "+" << endl;
  // else cout << "-" << endl;
  // cout << "(" << index_gen[0] << ", " << vec_gen.at(index_gen[0]).PID() <<  "), (" << index_gen[1] << ", " <<  vec_gen.at(index_gen[1]).PID() << ")" << endl;

  // for(unsigned int i=0; i<vec_jet.size(); i++)
  //   {
  //     int jet_flavour = vec_hf_flavour.at(i);
  //     int jet_origin = vec_hf_origin.at(i);

  //     cout << "(" << jet_flavour << "," <<  jet_origin << "), ";
  //   }
  // cout << endl;
  // for(int i=0; i<2; i++) cout << "Matched jet index = " << index_matched_jet[i] << ", surely matched = " << surely_matched[i] << ", d_r = " << dr_return[i] << endl;

  return selected_w;
} // int Vcb::Gen_Match_W(const vector<Jet>& vec_jet, const vector<Gen>& vec_gen, const vector<int>& vec_hf_flavour , const vector<int>& vec_hf_origin, int index_matched_jet[2])

//////////

void Vcb::Index_Converter(const vector<Jet> &vec_sel_jet, const vector<Jet> &vec_sel_jet_match, const int index_matched_jet_match[4], int index_matched_jet[4])
{
  for (int i = 0; i < 4; i++)
  {
    index_matched_jet[i] = -999;

    if (index_matched_jet_match[i] < 0)
      continue;

    Jet jet_match = vec_sel_jet_match[index_matched_jet_match[i]];

    float jet_match_pt = jet_match.Pt();
    float jet_match_eta = jet_match.Eta();
    float jet_match_phi = jet_match.Phi();

    for (unsigned int j = 0; j < vec_sel_jet.size(); j++)
    {
      Jet jet = vec_sel_jet.at(j);

      float jet_pt = jet.Pt();
      float jet_eta = jet.Eta();
      float jet_phi = jet.Phi();

      if (abs(jet_match_pt - jet_pt) < 1e-8 && abs(jet_match_eta - jet_eta) < 1e-8 && abs(jet_match_phi - jet_phi) < 1e-8)
      {
        index_matched_jet[i] = j;

        break;
      }
    } // for(unsigned int j=0; j<vec_sel_jet.size(); j++)
  }   // for(int i=0; i<4; i++)

  return;
} // void Vcb::Index_Converter(const vector<Jet>& vec_sel_jet, const vector<Jet>& vec_sel_jet_match, const int index_matched_jet_match[4], int index_matched_jet[4])

//////////

void Vcb::KF_Ambiguity_Remover(const vector<Jet> &vec_sel_jet, const int index_matched_jet[4])
{
  float mva_swapper = -999;
  float mva_swapper_swapped = -999;

  Results_Container results_container = fitter_driver->Get_Results();
  if (fitter_driver->Check_Status())
  {
    int best_index_had_t_b = results_container.best_index_had_t_b;
    int best_index_w_u = results_container.best_index_w_u;
    int best_index_w_d = results_container.best_index_w_d;
    int best_index_lep_t_b = results_container.best_index_lep_t_b;

    Jet jet_had_t_b = vec_sel_jet[best_index_had_t_b];
    Jet jet_w_u = vec_sel_jet[best_index_w_u];
    Jet jet_w_d = vec_sel_jet[best_index_w_d];
    Jet jet_lep_t_b = vec_sel_jet[best_index_lep_t_b];

    // tagging scores
    bvsc_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet);
    cvsb_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
    cvsl_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);

    bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
    cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
    cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);

    bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
    cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
    cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);

    bvsc_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet);
    cvsb_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
    cvsl_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);

    // BJetRegression should be considered so results from KF should be retrieved
    // pt
    pt_had_t_b = results_container.best_pt_had_t_b;
    pt_w_u = results_container.best_pt_w_u;
    pt_w_d = results_container.best_pt_w_d;
    pt_lep_t_b = results_container.best_pt_lep_t_b;

    // angles
    theta_w_u_w_d = results_container.best_theta_w_u_w_d;
    theta_had_w_had_t_b = results_container.best_theta_had_w_had_t_b;
    theta_lep_neu = results_container.best_theta_lep_neu;
    theta_lep_w_lep_t_b = results_container.best_theta_lep_w_lep_t_b;
    del_phi_had_t_lep_t = results_container.best_del_phi_had_t_lep_t;

    // masses
    had_t_mass = results_container.best_initial_had_t_mass;
    had_w_mass = results_container.best_initial_had_w_mass;
    lep_t_mass = results_container.best_initial_lep_t_mass;
    lep_t_partial_mass = results_container.best_initial_lep_t_partial_mass;

    // chi2s
    chi2_jet_had_t_b = results_container.best_chi2_jet_had_t_b;
    chi2_jet_w_u = results_container.best_chi2_jet_w_u;
    chi2_jet_w_d = results_container.best_chi2_jet_w_d;
    chi2_jet_lep_t_b = results_container.best_chi2_jet_lep_t_b;

    chi2_jet_extra = results_container.best_chi2_jet_extra;

    chi2_constraint_had_t = results_container.best_chi2_constraint_had_t;
    chi2_constraint_had_w = results_container.best_chi2_constraint_had_w;
    chi2_constraint_lep_t = results_container.best_chi2_constraint_lep_t;
    chi2_constraint_lep_w = results_container.best_chi2_constraint_lep_w;

    // swapper
    if (n_sel_jet == 4)
      mva_swapper = reader_swapper[0]->EvaluateMVA("BDTG_4Jets");
    else if (n_sel_jet == 5)
      mva_swapper = reader_swapper[1]->EvaluateMVA("BDTG_5Jets");
    else if (n_sel_jet == 6)
      mva_swapper = reader_swapper[1]->EvaluateMVA("BDTG_6Jets");
    else
      mva_swapper = reader_swapper[1]->EvaluateMVA("BDTG_7Jets");

    Results results;
    for (unsigned int i = 0; i < results_container.vec_results.size(); i++)
    {
      results = results_container.vec_results[i];

      int index_had_t_b = results.index_had_t_b;
      int index_w_u = results.index_w_u;
      int index_w_d = results.index_w_d;
      int index_lep_t_b = results.index_lep_t_b;

      if (best_index_had_t_b == index_had_t_b &&
          best_index_w_u == index_w_d &&
          best_index_w_d == index_w_u &&
          best_index_lep_t_b == index_lep_t_b)
      {
        // tagging scores are simply swapped
        float bvsc_temp = bvsc_w_u;
        float cvsb_temp = cvsb_w_u;
        float cvsl_temp = cvsl_w_u;

        bvsc_w_u = bvsc_w_d;
        cvsb_w_u = cvsb_w_d;
        cvsl_w_u = cvsl_w_d;

        bvsc_w_d = bvsc_temp;
        cvsb_w_d = cvsb_temp;
        cvsl_w_d = cvsl_temp;

        // BJetRegression should be considered so results from KF should be retrieved
        // pt
        pt_had_t_b = results.pt_had_t_b;
        pt_w_u = results.pt_w_u;
        pt_w_d = results.pt_w_d;
        pt_lep_t_b = results.pt_lep_t_b;

        // angles
        theta_w_u_w_d = results.theta_w_u_w_d;
        theta_had_w_had_t_b = results.theta_had_w_had_t_b;
        theta_lep_neu = results.theta_lep_neu;
        theta_lep_w_lep_t_b = results.theta_lep_w_lep_t_b;
        del_phi_had_t_lep_t = results.del_phi_had_t_lep_t;

        // masses
        had_t_mass = results.initial_had_t_mass;
        had_w_mass = results.initial_had_w_mass;
        lep_t_mass = results.initial_lep_t_mass;
        lep_t_partial_mass = results.initial_lep_t_partial_mass;

        // chi2s
        chi2_jet_had_t_b = results.chi2_jet_had_t_b;
        chi2_jet_w_u = results.chi2_jet_w_u;
        chi2_jet_w_d = results.chi2_jet_w_d;
        chi2_jet_lep_t_b = results.chi2_jet_lep_t_b;

        chi2_jet_extra = results.chi2_jet_extra;

        chi2_constraint_had_t = results.chi2_constraint_had_t;
        chi2_constraint_had_w = results.chi2_constraint_had_w;
        chi2_constraint_lep_t = results.chi2_constraint_lep_t;
        chi2_constraint_lep_w = results.chi2_constraint_lep_w;

        if (n_sel_jet == 4)
          mva_swapper_swapped = reader_swapper[0]->EvaluateMVA("BDTG_4Jets");
        else if (n_sel_jet == 5)
          mva_swapper_swapped = reader_swapper[1]->EvaluateMVA("BDTG_5Jets");
        else if (n_sel_jet == 6)
          mva_swapper_swapped = reader_swapper[1]->EvaluateMVA("BDTG_6Jets");
        else
          mva_swapper_swapped = reader_swapper[1]->EvaluateMVA("BDTG_7Jets");

        break;
      } // swapped case found
    }   // loop over vec_results

    if (mva_swapper < mva_swapper_swapped)
    {
      best_mva_score = mva_swapper_swapped;

      // update masses including BJetRegression
      results_container.best_initial_had_t_mass = results.initial_had_t_mass;
      results_container.best_initial_had_w_mass = results.initial_had_w_mass;
      results_container.best_initial_lep_t_mass = results.initial_lep_t_mass;
      results_container.best_initial_lep_t_partial_mass = results.initial_lep_t_partial_mass;
      results_container.best_initial_lep_w_mass = results.initial_lep_w_mass;

      // results_container.best_pt_w_u = ;
      // results_container.best_pt_w_d = ;

      swapped_mva = 1;
    }
    else
    {
      best_mva_score = mva_swapper;
      swapped_mva = 0;
    }
  } // if(fitter_driver->Check_Status())

  return;
} // void Vcb::KF_Ambiguity_Remover(const vector<Jet>& vec_sel_jet)

//////////

void Vcb::Make_HF_Contamination_Tree()
{
  // if(n_b_jet<3) return;

  Results_Container results_container = fitter_driver->Get_Results();
  if (!fitter_driver->Check_Status())
    return;

  chk_hf_contamination = false;

  int index[2];
  index[0] = results_container.best_index_w_u;
  index[1] = results_container.best_index_w_d;

  for (int i = 0; i < 2; i++)
  {
    Jet jet = vec_sel_jet[index[i]];

    // int hf_flavour = jet.GenHFHadronMatcherFlavour();
    int hf_origin = jet.GenHFHadronMatcherOrigin();

    if (hf_origin == 21)
      chk_hf_contamination = true;
  }

  // vec_bjet for easy handling
  vector<Jet> vec_bjet;
  for (int i = 0; i < n_sel_jet; i++)
  {
    if (vec_btag[i] == true)
      vec_bjet.push_back(vec_sel_jet[i]);
  }

  // vec_cjet for easy handling
  vector<Jet> vec_cjet;
  for (int i = 0; i < n_sel_jet; i++)
  {
    if (vec_ctag[i] == true)
      vec_cjet.push_back(vec_sel_jet[i]);
  }

  best_mva_score_pre = results_container.best_mva_score;
  del_phi_had_t_lep_t = results_container.best_del_phi_had_t_lep_t;

  theta_b_b = 999;
  for (int i = 0; i < n_b_jet; i++)
  {
    Jet b0 = vec_bjet[i];

    for (int j = i + 1; j < n_b_jet; j++)
    {
      Jet b1 = vec_bjet[j];

      float theta_b_b_new = b0.Angle(b1.Vect());

      if (theta_b_b_new < theta_b_b)
        theta_b_b = theta_b_b_new;
    }
  }

  // not useful...
  theta_c_c = 999;
  for (int i = 0; i < n_c_jet; i++)
  {
    Jet c0 = vec_cjet[i];

    for (int j = i + 1; j < n_c_jet; j++)
    {
      Jet c1 = vec_cjet[j];

      float theta_c_c_new = c0.Angle(c1.Vect());

      if (theta_c_c_new < theta_c_c)
        theta_c_c = theta_c_c_new;
    }
  }

  Jet jet_w_u = vec_sel_jet[index[0]];
  Jet jet_w_d = vec_sel_jet[index[1]];

  theta_p_had_w = results_container.best_theta_w_u_w_d * (jet_w_u + jet_w_d).P();

  for (int i = 0; i < 2; i++)
  {
    Jet jet_w_candi = vec_sel_jet[index[i]];

    int index_closest_b = -1;
    // int index_smallest_b = -1;
    float theta_b = 999;
    // float mass_b = 999;
    for (int j = 0; j < n_b_jet; j++)
    {
      Jet bjet = vec_bjet[j];

      // if two jets are same
      if (Compare_Jet(jet_w_candi, bjet))
        continue;

      float theta_b_new = jet_w_candi.Angle(bjet.Vect());
      // float mass_b_new = (jet_w_candi+bjet).M();

      if (theta_b_new < theta_b)
      {
        index_closest_b = j;
        theta_b = theta_b_new;
      }

      // if(mass_b_new<mass_b) mass_b = mass_b_new;
    } // loop over bjets

    if (i == 0)
    {
      theta_w_u_b = theta_b;

      w_u_b_bscore = vec_bjet[index_closest_b].GetTaggerResult(JetTagging::DeepJet);
      m_w_u_b = (jet_w_candi + vec_bjet[index_closest_b]).M();
    }
    else if (i == 1)
    {
      theta_w_d_b = theta_b;

      w_d_b_bscore = vec_bjet[index_closest_b].GetTaggerResult(JetTagging::DeepJet);
      m_w_d_b = (jet_w_candi + vec_bjet[index_closest_b]).M();
    }
  } // loop over w candidate

  if (chk_hf_contamination)
    hf_contamination_tree_correct->Fill();
  else
    hf_contamination_tree_wrong->Fill();

  return;
} // void Vcb::Make_HF_Contamination_Tree()

//////////

void Vcb::Make_Permutation_Tree()
{
  vector<Gen> vec_gen = GetGens();
  Gen gen_neutrino = Neutrino(vec_gen);

  gen_neutrino_px = gen_neutrino.Px();
  gen_neutrino_py = gen_neutrino.Py();
  gen_neutrino_pz = gen_neutrino.Pz();

  met_px = met.Px();
  met_py = met.Py();

  Results_Container results_container = fitter_driver->Get_Results();

  // search correct index first
  unsigned int index_correct = 999999; // should be very large number
  float diff_pz = 99999.;
  for (unsigned int i = 0; i < results_container.vec_results.size(); ++i)
  {
    Results results = results_container.vec_results[i];

    int index_had_t_b = results.index_had_t_b;
    int index_w_u = results.index_w_u;
    int index_w_d = results.index_w_d;
    int index_lep_t_b = results.index_lep_t_b;

    if (index_matched_jet[0] == index_had_t_b && index_matched_jet[1] == index_w_u &&
        index_matched_jet[2] == index_w_d && index_matched_jet[3] == index_lep_t_b)
    {
      float diff_current = TMath::Abs(results.neutrino_pz_sol - gen_neutrino_pz);

      if (diff_current < diff_pz)
      {
        index_correct = i;
        diff_pz = diff_current;
      }
    }
  } // search done

  // calculate neutrino pz solution with unrebalanced met and find best one
  float diff_pz_unrebal[2];
  float neu_pz_sol_unrebal[2];

  Sol_Neutrino_Pz(lepton, met, neu_pz_sol_unrebal);

  diff_pz_unrebal[0] = TMath::Abs(gen_neutrino_pz - neu_pz_sol_unrebal[0]);
  diff_pz_unrebal[1] = TMath::Abs(gen_neutrino_pz - neu_pz_sol_unrebal[1]);

  if (diff_pz_unrebal[0] < diff_pz_unrebal[1])
    neutrino_pz_sol_unrebal = neu_pz_sol_unrebal[0];
  else
    neutrino_pz_sol_unrebal = neu_pz_sol_unrebal[1];

  for (unsigned int i = 0; i < results_container.vec_results.size(); ++i)
  {
    Results results = results_container.vec_results[i];

    met_rebalance_px = results.met_rebalance_px;
    met_rebalance_py = results.met_rebalance_py;
    neutrino_pz_sol = results.neutrino_pz_sol;
    chk_real_neu_pz = results.chk_real_neu_pz;

    met_pt = sqrt(met_rebalance_px * met_rebalance_px + met_rebalance_py * met_rebalance_py);
    neutrino_p = sqrt(met_rebalance_px * met_rebalance_px + met_rebalance_py * met_rebalance_py + neutrino_pz_sol * neutrino_pz_sol);
    pt_ratio = lepton_pt / met_pt;

    mt_gen = Calculate_Mt(lepton, gen_neutrino_px, gen_neutrino_py);
    mt_met = Calculate_Mt(lepton, met_px, met_py);
    mt_met_rebalance = Calculate_Mt(lepton, met_rebalance_px, met_rebalance_py);

    int index_had_t_b = results.index_had_t_b;
    int index_w_u = results.index_w_u;
    int index_w_d = results.index_w_d;
    int index_lep_t_b = results.index_lep_t_b;

    n_matched_jets = 0;
    if (index_matched_jet[0] == index_had_t_b)
      n_matched_jets++;
    if (index_matched_jet[1] == index_w_u)
      n_matched_jets++;
    if (index_matched_jet[2] == index_w_d)
      n_matched_jets++;
    if (index_matched_jet[3] == index_lep_t_b)
      n_matched_jets++;

    Jet jet_had_t_b = vec_sel_jet[index_had_t_b];
    Jet jet_w_u = vec_sel_jet[index_w_u];
    Jet jet_w_d = vec_sel_jet[index_w_d];
    Jet jet_lep_t_b = vec_sel_jet[index_lep_t_b];

    // tagging scores, no jet energy regnession applied
    bvsc_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet);
    cvsb_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
    cvsl_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);

    bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
    cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
    cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);

    bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
    cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
    cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);

    bvsc_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet);
    cvsb_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
    cvsl_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);

    // pt
    pt_had_t_b = results.pt_had_t_b;
    pt_w_u = results.pt_w_u;
    pt_w_d = results.pt_w_d;
    pt_lep_t_b = results.pt_lep_t_b;

    // eta
    eta_had_t_b = results.eta_had_t_b;
    eta_w_u = results.eta_w_u;
    eta_w_d = results.eta_w_d;
    eta_lep_t_b = results.eta_lep_t_b;

    // angles
    del_phi_w_u_w_d = results.del_phi_w_u_w_d;
    del_phi_had_w_had_t_b = results.del_phi_had_w_had_t_b;
    del_phi_lep_neu = results.del_phi_lep_neu;
    del_phi_lep_w_lep_t_b = results.del_phi_lep_w_lep_t_b;
    del_phi_had_t_lep_t = results.del_phi_had_t_lep_t;

    del_eta_w_u_w_d = results.del_eta_w_u_w_d;
    del_eta_had_w_had_t_b = results.del_eta_had_w_had_t_b;
    del_eta_lep_neu = results.del_eta_lep_neu;
    del_eta_lep_w_lep_t_b = results.del_eta_lep_w_lep_t_b;
    del_eta_had_t_lep_t = results.del_eta_had_t_lep_t;

    del_r_w_u_w_d = results.del_r_w_u_w_d;
    del_r_had_w_had_t_b = results.del_r_had_w_had_t_b;
    del_r_lep_neu = results.del_r_lep_neu;
    del_r_lep_w_lep_t_b = results.del_r_lep_w_lep_t_b;
    del_r_had_t_lep_t = results.del_r_had_t_lep_t;

    theta_w_u_w_d = results.theta_w_u_w_d;
    theta_had_w_had_t_b = results.theta_had_w_had_t_b;
    theta_lep_neu = results.theta_lep_neu;
    theta_lep_w_lep_t_b = results.theta_lep_w_lep_t_b;
    theta_had_t_lep_t = results.theta_had_t_lep_t;

    // masses
    had_t_mass = results.initial_had_t_mass;
    had_w_mass = results.initial_had_w_mass;
    lep_t_mass = results.initial_lep_t_mass;

    lep_t_partial_mass = results.initial_lep_t_partial_mass;

    // chi2
    chi2_jet_had_t_b = results.chi2_jet_had_t_b;
    chi2_jet_w_u = results.chi2_jet_w_u;
    chi2_jet_w_d = results.chi2_jet_w_d;
    chi2_jet_lep_t_b = results.chi2_jet_lep_t_b;
    chi2_jet_extra = results.chi2_jet_extra;

    chi2_constraint_had_t = results.chi2_constraint_had_t;
    chi2_constraint_had_w = results.chi2_constraint_had_w;
    chi2_constraint_lep_t = results.chi2_constraint_lep_t;
    chi2_constraint_lep_w = results.chi2_constraint_lep_w;

    chi2 = results.chi2;

    if (i == index_correct && chk_gentau_conta == false)
      permutation_tree_correct->Fill();
    else
      permutation_tree_wrong->Fill();
  } // for(unsigned int i=0; i<result_container.vec_results.size(); ++i)

  return;
} // void Vcb::Make_Permutation_Tree()

//////////

void Vcb::Make_Result_Tree(const AnalyzerParameter &param)
{
  bool chk_fitter_status = fitter_driver->Check_Status();
  if (!chk_fitter_status)
    return;

  FillHist(param.Name + "/Cut_Flow", Cut_Flow::KF_Pass, weight, n_cut_flow, 0, n_cut_flow);
  FillHist(param.Name + Form("/Cut_Flow_%d", decay_mode), Cut_Flow::KF_Pass, weight, n_cut_flow, 0, n_cut_flow);

  Jet leading_jet = vec_sel_jet.at(0);
  Jet subleading_jet = vec_sel_jet.at(1);

  pt_leading_jet = leading_jet.Pt();
  pt_subleading_jet = subleading_jet.Pt();

  eta_leading_jet = leading_jet.Eta();
  eta_subleading_jet = subleading_jet.Eta();

  bvsc_leading_jet = leading_jet.GetTaggerResult(JetTagging::DeepJet);
  cvsb_leading_jet = leading_jet.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_leading_jet = leading_jet.GetTaggerResult(JetTagging::DeepJet_CvsL);

  bvsc_subleading_jet = subleading_jet.GetTaggerResult(JetTagging::DeepJet);
  cvsb_subleading_jet = subleading_jet.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_subleading_jet = subleading_jet.GetTaggerResult(JetTagging::DeepJet_CvsL);

  ht = Calculate_HT(vec_sel_jet);

  Results_Container results_container = fitter_driver->Get_Results();
  best_mva_score_pre = results_container.best_mva_score;
  best_chi2 = results_container.best_chi2;

  met_px = met.Px();
  met_py = met.Py();
  mt = Calculate_Mt(lepton, met_px, met_py);

  int index_had_t_b = results_container.best_index_had_t_b;
  int index_w_u = results_container.best_index_w_u;
  int index_w_d = results_container.best_index_w_d;
  int index_lep_t_b = results_container.best_index_lep_t_b;

  Jet jet_had_t_b = vec_sel_jet.at(index_had_t_b);
  Jet jet_w_u = vec_sel_jet.at(index_w_u);
  Jet jet_w_d = vec_sel_jet.at(index_w_d);
  Jet jet_lep_t_b = vec_sel_jet.at(index_lep_t_b);

  bvsc_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet);
  cvsb_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_had_t_b = jet_had_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);

  bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
  cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);

  bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
  cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);

  bvsc_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet);
  cvsb_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_lep_t_b = jet_lep_t_b.GetTaggerResult(JetTagging::DeepJet_CvsL);

  // should be BJetRegression corrected pt
  pt_had_t_b = results_container.best_pt_had_t_b;
  pt_w_u = results_container.best_pt_w_u;
  pt_w_d = results_container.best_pt_w_d;
  pt_lep_t_b = results_container.best_pt_lep_t_b;

  eta_had_t_b = jet_had_t_b.Eta();
  eta_w_u = jet_w_u.Eta();
  eta_w_d = jet_w_d.Eta();
  eta_lep_t_b = jet_lep_t_b.Eta();

  // should be BJetRegression corrected pt
  m_had_t = results_container.best_initial_had_t_mass;
  m_had_w = results_container.best_initial_had_w_mass;
  m_lep_t = results_container.best_initial_lep_t_mass;
  m_lep_w = results_container.best_initial_lep_w_mass;

  m_w_u = jet_w_u.GetM();
  m_w_d = jet_w_d.GetM();

  // For MC
  chk_reco_correct = false;
  chk_included = false;
  chk_hf_contamination = false;
  pu_conta_had_t_b = false;
  pu_conta_w_u = false;
  pu_conta_w_d = false;
  pu_conta_lep_t_b = false;
  swapped_truth = -1;

  if (!IsDATA)
  {
    // chk_reco_correct & swapped_truth
    if (index_matched_jet[1] == index_w_u && index_matched_jet[2] == index_w_d)
    {
      chk_reco_correct = true;
      swapped_truth = 0;
    }
    else if (index_matched_jet[1] == index_w_d && index_matched_jet[2] == index_w_u)
    {
      chk_reco_correct = true;
      swapped_truth = 1;
    }
    else
    {
      chk_reco_correct = false;
      swapped_truth = -1;
    }

    // chk_included
    if (Chk_Included(index_matched_jet) == 0)
      chk_included = true;
    else
      chk_included = false;

    // chk_hf_contamination
    chk_hf_contamination = false;

    int hf_origin = jet_w_u.GenHFHadronMatcherOrigin();
    if (hf_origin == 21)
      chk_hf_contamination = true;

    hf_origin = jet_w_d.GenHFHadronMatcherOrigin();
    if (hf_origin == 21)
      chk_hf_contamination = true;

    if (jet_had_t_b.GenHFHadronMatcherFlavour() == -999 && jet_had_t_b.GenHFHadronMatcherOrigin() == -999)
      pu_conta_had_t_b = true;
    else
      pu_conta_had_t_b = false;
    if (jet_w_u.GenHFHadronMatcherFlavour() == -999 && jet_w_u.GenHFHadronMatcherOrigin() == -999)
      pu_conta_w_u = true;
    else
      pu_conta_w_u = false;
    if (jet_w_d.GenHFHadronMatcherFlavour() == -999 && jet_w_d.GenHFHadronMatcherOrigin() == -999)
      pu_conta_w_d = true;
    else
      pu_conta_w_d = false;
    if (jet_lep_t_b.GenHFHadronMatcherFlavour() == -999 && jet_lep_t_b.GenHFHadronMatcherOrigin() == -999)
      pu_conta_lep_t_b = true;
    else
      pu_conta_lep_t_b = false;
  } // For MC

  map_result_tree[param.syst_]->Fill();

  Set_Region();

  FillHist(param.Name + "/" + region + "/N_Vertex", nPV, weight, 100, 0, 100);
  FillHist(param.Name + "/" + region + "/Lepton_Pt", lepton_pt, weight, 50, 0, 300);
  FillHist(param.Name + "/" + region + "/Lepton_Eta", lepton_eta, weight, 50, -3, 3);
  FillHist(param.Name + "/" + region + "/N_Jets", n_sel_jet, weight, 30, 0, 30);
  FillHist(param.Name + "/" + region + "/N_BJets", n_b_jet, weight, 30, 0, 30);
  FillHist(param.Name + "/" + region + "/N_CJets", n_c_jet, weight, 30, 0, 30);
  FillHist(param.Name + "/" + region + "/Pt_Leading_Jet", pt_leading_jet, weight, 50, 0, 300);
  FillHist(param.Name + "/" + region + "/Pt_Subleading_Jet", pt_subleading_jet, weight, 50, 0, 300);
  FillHist(param.Name + "/" + region + "/Eta_Leading_Jet", eta_leading_jet, weight, 50, -3, 3);
  FillHist(param.Name + "/" + region + "/Eta_Subleading_Jet", eta_subleading_jet, weight, 50, -3, 3);
  FillHist(param.Name + "/" + region + "/Met_Pt", met_pt, weight, 50, 0, 300);
  FillHist(param.Name + "/" + region + "/Met_Phi", met_phi, weight, 80, -4, 4);
  FillHist(param.Name + "/" + region + "/Best_MVA_Score_Pre", best_mva_score_pre, weight, 110, -1.1, 1.1);
  FillHist(param.Name + "/" + region + "/Best_MVA_Score", best_mva_score, weight, 110, -1.1, 1.1);

  return;
} // void Vcb::Make_Result_Tree()

//////////

void Vcb::Make_Template_Tree()
{
  Results_Container results_container = fitter_driver->Get_Results();
  if (!fitter_driver->Check_Status())
    return;

  best_mva_score_pre = results_container.best_mva_score;

  int index_had_t_b = results_container.best_index_had_t_b;
  int index_w_u = results_container.best_index_w_u;
  int index_w_d = results_container.best_index_w_d;
  int index_lep_t_b = results_container.best_index_lep_t_b;

  if (index_matched_jet[0] != index_had_t_b || index_matched_jet[1] != index_w_u ||
      index_matched_jet[2] != index_w_d || index_matched_jet[3] != index_lep_t_b)
    return;

  Jet jet_w_u = vec_sel_jet[index_w_u];
  Jet jet_w_d = vec_sel_jet[index_w_d];

  bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
  cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);

  bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
  cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);

  int index = -1;
  if (decay_mode == 21)
    index = 0;
  else if (decay_mode == 23)
    index = 1;
  else if (decay_mode == 41)
    index = 2;
  else if (decay_mode == 43)
    index = 3;
  else if (decay_mode == 45)
    index = 4;
  if (index != -1)
    template_tree[index]->Fill();

  return;
} // void Vcb::Make_Template_Tree()

//////////

void Vcb::Make_Template_Truth_Tree()
{
  if (index_matched_jet[1] < 0 || index_matched_jet[2] < 0)
    return;

  Jet jet_w_u = vec_sel_jet[index_matched_jet[1]];
  Jet jet_w_d = vec_sel_jet[index_matched_jet[2]];

  best_mva_score_pre = 1;

  bvsc_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet);
  cvsb_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_w_u = jet_w_u.GetTaggerResult(JetTagging::DeepJet_CvsL);
  m_w_u = jet_w_u.GetM();

  bvsc_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet);
  cvsb_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsB);
  cvsl_w_d = jet_w_d.GetTaggerResult(JetTagging::DeepJet_CvsL);
  m_w_d = jet_w_d.GetM();

  int index = -1;
  if (decay_mode == 21)
    index = 0;
  else if (decay_mode == 23)
    index = 1;
  else if (decay_mode == 41)
    index = 2;
  else if (decay_mode == 43)
    index = 3;
  else if (decay_mode == 45)
    index = 4;
  if (index != -1)
    template_truth_tree[index]->Fill();

  return;
} // void Vcb::Make_Template_Truth_Tree()

//////////

Gen Vcb::Neutrino(const vector<Gen> &vec_gen)
{
  // find neutrino from W
  int index_last_w = -999;
  int index_last_aw = -999;
  int index_d_w[4] = {-999, -999, -999, -999};

  for (unsigned int i = 0; i < vec_gen.size(); ++i)
  {
    Gen gen = vec_gen[i];

    int pid = gen.PID();
    int m_index = gen.MotherIndex();

    // find last index of W+ and W-
    if (pid == 24)
      index_last_w = i;
    if (pid == -24)
      index_last_aw = i;

    // find decay products of W+
    if (m_index == index_last_w && pid != 24)
    {
      if (index_d_w[0] == -999)
        index_d_w[0] = i;
      else
        index_d_w[1] = i;
    }

    // find decay products of W-
    if (m_index == index_last_aw && pid != -24)
    {
      if (index_d_w[2] == -999)
        index_d_w[2] = i;
      else
        index_d_w[3] = i;
    }
  } // for(unsigned int i=0; i<vec_gen.size(); ++i)

  // Find neutrino within daughters of leptonic decayed W
  int index_neutrino = 0;
  for (int i = 0; i < 4; ++i)
  {
    if (Abs(vec_gen[index_d_w[i]].PID()) == 12 || Abs(vec_gen[index_d_w[i]].PID()) == 14)
    {
      index_neutrino = index_d_w[i];
      break;
    }
  }

  Gen gen = vec_gen[index_neutrino];

  return gen;
} // Gen Neutrino(const vector<Gen>& vec_gen)

//////////

Particle Vcb::Rebalance_Met()
{
  // cout << "before " << met.Pt() << " " << met.Phi() << endl;

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

  // cout << "after " << met_rebal.Pt() << " " << met_rebal.Phi() << endl;

  return met_rebal;
} // Particle Vcb::Rebalance_Met()

//////////

void Vcb::Set_HF_Contamination_Tree()
{
  hf_contamination_tree_correct = new TTree("HF_Contamination_Tree_Correct", "HF_Contamination_Tree_Correct");
  hf_contamination_tree_correct->Branch("n_jets", &n_sel_jet);
  hf_contamination_tree_correct->Branch("n_bjets", &n_b_jet);
  hf_contamination_tree_correct->Branch("n_cjets", &n_c_jet);
  hf_contamination_tree_correct->Branch("best_mva_score_pre", &best_mva_score_pre);
  hf_contamination_tree_correct->Branch("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
  hf_contamination_tree_correct->Branch("theta_b_b", &theta_b_b);
  hf_contamination_tree_correct->Branch("theta_c_c", &theta_c_c);
  hf_contamination_tree_correct->Branch("theta_p_had_w", &theta_p_had_w);
  hf_contamination_tree_correct->Branch("theta_w_u_b", &theta_w_u_b);
  hf_contamination_tree_correct->Branch("theta_w_d_b", &theta_w_d_b);
  hf_contamination_tree_correct->Branch("w_u_b_bscore", &w_u_b_bscore);
  hf_contamination_tree_correct->Branch("w_d_b_bscore", &w_d_b_bscore);
  hf_contamination_tree_correct->Branch("m_w_u_b", &m_w_u_b);
  hf_contamination_tree_correct->Branch("m_w_d_b", &m_w_d_b);

  hf_contamination_tree_wrong = new TTree("HF_Contamination_Tree_Wrong", "HF_Contamination_Tree_Wrong");
  hf_contamination_tree_wrong->Branch("n_jets", &n_sel_jet);
  hf_contamination_tree_wrong->Branch("n_bjets", &n_b_jet);
  hf_contamination_tree_wrong->Branch("n_cjets", &n_c_jet);
  hf_contamination_tree_wrong->Branch("best_mva_score_pre", &best_mva_score_pre);
  hf_contamination_tree_wrong->Branch("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
  hf_contamination_tree_wrong->Branch("theta_p_had_w", &theta_p_had_w);
  hf_contamination_tree_wrong->Branch("theta_b_b", &theta_b_b);
  hf_contamination_tree_wrong->Branch("theta_c_c", &theta_c_c);
  hf_contamination_tree_wrong->Branch("theta_w_u_b", &theta_w_u_b);
  hf_contamination_tree_wrong->Branch("theta_w_d_b", &theta_w_d_b);
  hf_contamination_tree_wrong->Branch("theta_b_b", &theta_b_b);
  hf_contamination_tree_wrong->Branch("w_u_b_bscore", &w_u_b_bscore);
  hf_contamination_tree_wrong->Branch("w_d_b_bscore", &w_d_b_bscore);
  hf_contamination_tree_wrong->Branch("m_w_u_b", &m_w_u_b);
  hf_contamination_tree_wrong->Branch("m_w_d_b", &m_w_d_b);

  return;
} // void Vcb::Set_HF_Contamination_Tree()

//////////

void Vcb::Set_Permutation_Tree()
{
  permutation_tree_correct = new TTree("Permutation_Correct", "Permutation_Correct");
  permutation_tree_correct->Branch("weight", &weight);
  permutation_tree_correct->Branch("n_jets", &n_sel_jet);
  permutation_tree_correct->Branch("n_bjets", &n_b_jet);
  permutation_tree_correct->Branch("n_cjets", &n_c_jet);
  permutation_tree_correct->Branch("lepton_pt", &lepton_pt);
  permutation_tree_correct->Branch("pt_ratio", &pt_ratio);
  permutation_tree_correct->Branch("chk_gentau_conta", &chk_gentau_conta);
  permutation_tree_correct->Branch("n_matched_jets", &n_matched_jets);
  permutation_tree_correct->Branch("gen_neutrino_px", &gen_neutrino_px);
  permutation_tree_correct->Branch("gen_neutrino_py", &gen_neutrino_py);
  permutation_tree_correct->Branch("gen_neutrino_pz", &gen_neutrino_pz);
  permutation_tree_correct->Branch("met_pt", &met_pt);
  permutation_tree_correct->Branch("met_px", &met_px);
  permutation_tree_correct->Branch("met_py", &met_py);
  permutation_tree_correct->Branch("met_rebalance_px", &met_rebalance_px);
  permutation_tree_correct->Branch("met_rebalance_py", &met_rebalance_py);
  permutation_tree_correct->Branch("neutrino_pz_sol", &neutrino_pz_sol);
  permutation_tree_correct->Branch("neutrino_pz_sol_unrebal", &neutrino_pz_sol_unrebal);
  permutation_tree_correct->Branch("neutrino_p", &neutrino_p);
  permutation_tree_correct->Branch("chk_real_neu_pz", &chk_real_neu_pz);
  permutation_tree_correct->Branch("mt_gen", &mt_gen);
  permutation_tree_correct->Branch("mt_met", &mt_met);
  permutation_tree_correct->Branch("mt_met_rebalance", &mt_met_rebalance);
  permutation_tree_correct->Branch("pt_had_t_b", &pt_had_t_b);
  permutation_tree_correct->Branch("pt_w_u", &pt_w_u);
  permutation_tree_correct->Branch("pt_w_d", &pt_w_d);
  permutation_tree_correct->Branch("pt_lep_t_b", &pt_lep_t_b);
  permutation_tree_correct->Branch("eta_had_t_b", &eta_had_t_b);
  permutation_tree_correct->Branch("eta_w_u", &eta_w_u);
  permutation_tree_correct->Branch("eta_w_d", &eta_w_d);
  permutation_tree_correct->Branch("eta_lep_t_b", &eta_lep_t_b);
  permutation_tree_correct->Branch("bvsc_had_t_b", &bvsc_had_t_b);
  permutation_tree_correct->Branch("cvsb_had_t_b", &cvsb_had_t_b);
  permutation_tree_correct->Branch("cvsl_had_t_b", &cvsl_had_t_b);
  permutation_tree_correct->Branch("bvsc_w_u", &bvsc_w_u);
  permutation_tree_correct->Branch("cvsb_w_u", &cvsb_w_u);
  permutation_tree_correct->Branch("cvsl_w_u", &cvsl_w_u);
  permutation_tree_correct->Branch("bvsc_w_d", &bvsc_w_d);
  permutation_tree_correct->Branch("cvsb_w_d", &cvsb_w_d);
  permutation_tree_correct->Branch("cvsl_w_d", &cvsl_w_d);
  permutation_tree_correct->Branch("bvsc_lep_t_b", &bvsc_lep_t_b);
  permutation_tree_correct->Branch("cvsb_lep_t_b", &cvsb_lep_t_b);
  permutation_tree_correct->Branch("cvsl_lep_t_b", &cvsl_lep_t_b);
  permutation_tree_correct->Branch("met_rebalance_px", &met_rebalance_px);
  permutation_tree_correct->Branch("met_rebalance_py", &met_rebalance_py);
  permutation_tree_correct->Branch("del_phi_w_u_w_d", &del_phi_w_u_w_d);
  permutation_tree_correct->Branch("del_phi_had_w_had_t_b", &del_phi_had_w_had_t_b);
  permutation_tree_correct->Branch("del_phi_lep_neu", &del_phi_lep_neu);
  permutation_tree_correct->Branch("del_phi_lep_w_lep_t_b", &del_phi_lep_w_lep_t_b);
  permutation_tree_correct->Branch("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
  permutation_tree_correct->Branch("del_eta_w_u_w_d", &del_eta_w_u_w_d);
  permutation_tree_correct->Branch("del_eta_had_w_had_t_b", &del_eta_had_w_had_t_b);
  permutation_tree_correct->Branch("del_eta_lep_neu", &del_eta_lep_neu);
  permutation_tree_correct->Branch("del_eta_lep_w_lep_t_b", &del_eta_lep_w_lep_t_b);
  permutation_tree_correct->Branch("del_eta_had_t_lep_t", &del_eta_had_t_lep_t);
  permutation_tree_correct->Branch("del_r_w_u_w_d", &del_r_w_u_w_d);
  permutation_tree_correct->Branch("del_r_had_w_had_t_b", &del_r_had_w_had_t_b);
  permutation_tree_correct->Branch("del_r_lep_neu", &del_r_lep_neu);
  permutation_tree_correct->Branch("del_r_lep_w_lep_t_b", &del_r_lep_w_lep_t_b);
  permutation_tree_correct->Branch("del_r_had_t_lep_t", &del_r_had_t_lep_t);
  permutation_tree_correct->Branch("theta_w_u_w_d", &theta_w_u_w_d);
  permutation_tree_correct->Branch("theta_had_w_had_t_b", &theta_had_w_had_t_b);
  permutation_tree_correct->Branch("theta_lep_neu", &theta_lep_neu);
  permutation_tree_correct->Branch("theta_lep_w_lep_t_b", &theta_lep_w_lep_t_b);
  permutation_tree_correct->Branch("theta_had_t_lep_t", &theta_had_t_lep_t);
  permutation_tree_correct->Branch("had_t_mass", &had_t_mass);
  permutation_tree_correct->Branch("had_w_mass", &had_w_mass);
  permutation_tree_correct->Branch("lep_t_mass", &lep_t_mass);
  permutation_tree_correct->Branch("lep_t_partial_mass", &lep_t_partial_mass);
  permutation_tree_correct->Branch("chi2_jet_had_t_b", &chi2_jet_had_t_b);
  permutation_tree_correct->Branch("chi2_jet_w_u", &chi2_jet_w_u);
  permutation_tree_correct->Branch("chi2_jet_w_d", &chi2_jet_w_d);
  permutation_tree_correct->Branch("chi2_jet_lep_t_b", &chi2_jet_lep_t_b);
  permutation_tree_correct->Branch("chi2_jet_extra", &chi2_jet_extra);
  permutation_tree_correct->Branch("chi2_constraint_had_t", &chi2_constraint_had_t);
  permutation_tree_correct->Branch("chi2_constraint_had_w", &chi2_constraint_had_w);
  permutation_tree_correct->Branch("chi2_constraint_lep_t", &chi2_constraint_lep_t);
  permutation_tree_correct->Branch("chi2_constraint_lep_w", &chi2_constraint_lep_w);
  permutation_tree_correct->Branch("chi2", &chi2);

  permutation_tree_wrong = new TTree("Permutation_Wrong", "Permutation_Wrong");
  permutation_tree_wrong->Branch("weight", &weight);
  permutation_tree_wrong->Branch("n_jets", &n_sel_jet);
  permutation_tree_wrong->Branch("n_bjets", &n_b_jet);
  permutation_tree_wrong->Branch("n_cjets", &n_c_jet);
  permutation_tree_wrong->Branch("lepton_pt", &lepton_pt);
  permutation_tree_wrong->Branch("pt_ratio", &pt_ratio);
  permutation_tree_wrong->Branch("chk_gentau_conta", &chk_gentau_conta);
  permutation_tree_wrong->Branch("n_matched_jets", &n_matched_jets);
  permutation_tree_wrong->Branch("gen_neutrino_px", &gen_neutrino_px);
  permutation_tree_wrong->Branch("gen_neutrino_py", &gen_neutrino_py);
  permutation_tree_wrong->Branch("gen_neutrino_pz", &gen_neutrino_pz);
  permutation_tree_wrong->Branch("met_pt", &met_pt);
  permutation_tree_wrong->Branch("met_px", &met_px);
  permutation_tree_wrong->Branch("met_py", &met_py);
  permutation_tree_wrong->Branch("met_rebalance_px", &met_rebalance_px);
  permutation_tree_wrong->Branch("met_rebalance_py", &met_rebalance_py);
  permutation_tree_wrong->Branch("neutrino_pz_sol", &neutrino_pz_sol);
  permutation_tree_wrong->Branch("neutrino_pz_sol_unrebal", &neutrino_pz_sol_unrebal);
  permutation_tree_wrong->Branch("neutrino_p", &neutrino_p);
  permutation_tree_wrong->Branch("chk_real_neu_pz", &chk_real_neu_pz);
  permutation_tree_wrong->Branch("mt_gen", &mt_gen);
  permutation_tree_wrong->Branch("mt_met", &mt_met);
  permutation_tree_wrong->Branch("mt_met_rebalance", &mt_met_rebalance);
  permutation_tree_wrong->Branch("pt_had_t_b", &pt_had_t_b);
  permutation_tree_wrong->Branch("pt_w_u", &pt_w_u);
  permutation_tree_wrong->Branch("pt_w_d", &pt_w_d);
  permutation_tree_wrong->Branch("pt_lep_t_b", &pt_lep_t_b);
  permutation_tree_wrong->Branch("eta_had_t_b", &eta_had_t_b);
  permutation_tree_wrong->Branch("eta_w_u", &eta_w_u);
  permutation_tree_wrong->Branch("eta_w_d", &eta_w_d);
  permutation_tree_wrong->Branch("eta_lep_t_b", &eta_lep_t_b);
  permutation_tree_wrong->Branch("bvsc_had_t_b", &bvsc_had_t_b);
  permutation_tree_wrong->Branch("cvsb_had_t_b", &cvsb_had_t_b);
  permutation_tree_wrong->Branch("cvsl_had_t_b", &cvsl_had_t_b);
  permutation_tree_wrong->Branch("bvsc_w_u", &bvsc_w_u);
  permutation_tree_wrong->Branch("cvsb_w_u", &cvsb_w_u);
  permutation_tree_wrong->Branch("cvsl_w_u", &cvsl_w_u);
  permutation_tree_wrong->Branch("bvsc_w_d", &bvsc_w_d);
  permutation_tree_wrong->Branch("cvsb_w_d", &cvsb_w_d);
  permutation_tree_wrong->Branch("cvsl_w_d", &cvsl_w_d);
  permutation_tree_wrong->Branch("bvsc_lep_t_b", &bvsc_lep_t_b);
  permutation_tree_wrong->Branch("cvsb_lep_t_b", &cvsb_lep_t_b);
  permutation_tree_wrong->Branch("cvsl_lep_t_b", &cvsl_lep_t_b);
  permutation_tree_wrong->Branch("met_rebalance_px", &met_rebalance_px);
  permutation_tree_wrong->Branch("met_rebalance_py", &met_rebalance_py);
  permutation_tree_wrong->Branch("del_phi_w_u_w_d", &del_phi_w_u_w_d);
  permutation_tree_wrong->Branch("del_phi_had_w_had_t_b", &del_phi_had_w_had_t_b);
  permutation_tree_wrong->Branch("del_phi_lep_neu", &del_phi_lep_neu);
  permutation_tree_wrong->Branch("del_phi_lep_w_lep_t_b", &del_phi_lep_w_lep_t_b);
  permutation_tree_wrong->Branch("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
  permutation_tree_wrong->Branch("del_eta_w_u_w_d", &del_eta_w_u_w_d);
  permutation_tree_wrong->Branch("del_eta_had_w_had_t_b", &del_eta_had_w_had_t_b);
  permutation_tree_wrong->Branch("del_eta_lep_neu", &del_eta_lep_neu);
  permutation_tree_wrong->Branch("del_eta_lep_w_lep_t_b", &del_eta_lep_w_lep_t_b);
  permutation_tree_wrong->Branch("del_eta_had_t_lep_t", &del_eta_had_t_lep_t);
  permutation_tree_wrong->Branch("del_r_w_u_w_d", &del_r_w_u_w_d);
  permutation_tree_wrong->Branch("del_r_had_w_had_t_b", &del_r_had_w_had_t_b);
  permutation_tree_wrong->Branch("del_r_lep_neu", &del_r_lep_neu);
  permutation_tree_wrong->Branch("del_r_lep_w_lep_t_b", &del_r_lep_w_lep_t_b);
  permutation_tree_wrong->Branch("del_r_had_t_lep_t", &del_r_had_t_lep_t);
  permutation_tree_wrong->Branch("theta_w_u_w_d", &theta_w_u_w_d);
  permutation_tree_wrong->Branch("theta_had_w_had_t_b", &theta_had_w_had_t_b);
  permutation_tree_wrong->Branch("theta_lep_neu", &theta_lep_neu);
  permutation_tree_wrong->Branch("theta_lep_w_lep_t_b", &theta_lep_w_lep_t_b);
  permutation_tree_wrong->Branch("theta_had_t_lep_t", &theta_had_t_lep_t);
  permutation_tree_wrong->Branch("had_t_mass", &had_t_mass);
  permutation_tree_wrong->Branch("had_w_mass", &had_w_mass);
  permutation_tree_wrong->Branch("lep_t_mass", &lep_t_mass);
  permutation_tree_wrong->Branch("lep_t_partial_mass", &lep_t_partial_mass);
  permutation_tree_wrong->Branch("chi2_jet_had_t_b", &chi2_jet_had_t_b);
  permutation_tree_wrong->Branch("chi2_jet_w_u", &chi2_jet_w_u);
  permutation_tree_wrong->Branch("chi2_jet_w_d", &chi2_jet_w_d);
  permutation_tree_wrong->Branch("chi2_jet_lep_t_b", &chi2_jet_lep_t_b);
  permutation_tree_wrong->Branch("chi2_jet_extra", &chi2_jet_extra);
  permutation_tree_wrong->Branch("chi2_constraint_had_t", &chi2_constraint_had_t);
  permutation_tree_wrong->Branch("chi2_constraint_had_w", &chi2_constraint_had_w);
  permutation_tree_wrong->Branch("chi2_constraint_lep_t", &chi2_constraint_lep_t);
  permutation_tree_wrong->Branch("chi2_constraint_lep_w", &chi2_constraint_lep_w);
  permutation_tree_wrong->Branch("chi2", &chi2);

  return;
} // void Vcb::Set_Permutation_Tree()

//////////

void Vcb::Set_Reader_HF_Contamination()
{
  reader_hf_contamination_lessthantwo = new TMVA::Reader("!Color:!Silent");
  reader_hf_contamination_lessthantwo->AddSpectator("n_jets", &n_sel_jet);
  reader_hf_contamination_lessthantwo->AddVariable("n_bjets", &n_b_jet_f); // strange...
  reader_hf_contamination_lessthantwo->AddVariable("n_cjets", &n_c_jet_f);
  reader_hf_contamination_lessthantwo->AddVariable("best_mva_score", &best_mva_score_pre);
  reader_hf_contamination_lessthantwo->AddVariable("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
  reader_hf_contamination_lessthantwo->AddVariable("theta_b_b", &theta_b_b);
  reader_hf_contamination_lessthantwo->AddVariable("theta_p_had_w", &theta_p_had_w);
  reader_hf_contamination_lessthantwo->AddVariable("theta_w_u_b", &theta_w_u_b);
  reader_hf_contamination_lessthantwo->AddVariable("theta_w_d_b", &theta_w_d_b);
  reader_hf_contamination_lessthantwo->AddVariable("w_u_b_bscore", &w_u_b_bscore);
  reader_hf_contamination_lessthantwo->AddVariable("w_d_b_bscore", &w_d_b_bscore);
  reader_hf_contamination_lessthantwo->AddVariable("m_w_u_b", &m_w_u_b);
  reader_hf_contamination_lessthantwo->AddVariable("m_w_d_b", &m_w_d_b);

  TString weight_file_base = getenv("SKFlat_WD");
  weight_file_base += "/data/Run2UltraLegacy_v3/";
  weight_file_base += DataEra;
  weight_file_base += "/HF_Contamination/";

  TString weight_file = weight_file_base + "LessThanTwo/weights/TMVAClassification_BDTG.weights.xml";

  cout << weight_file << endl;
  reader_hf_contamination_lessthantwo->BookMVA("LessThanTwo", weight_file);

  reader_hf_contamination_morethantwo = new TMVA::Reader("!Color:!Silent");
  reader_hf_contamination_morethantwo->AddSpectator("n_jets", &n_sel_jet);
  reader_hf_contamination_morethantwo->AddVariable("n_bjets", &n_b_jet_f); // strange...
  reader_hf_contamination_morethantwo->AddVariable("n_cjets", &n_c_jet_f);
  reader_hf_contamination_morethantwo->AddVariable("best_mva_score", &best_mva_score_pre);
  reader_hf_contamination_morethantwo->AddVariable("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);
  reader_hf_contamination_morethantwo->AddVariable("theta_b_b", &theta_b_b);
  reader_hf_contamination_morethantwo->AddVariable("theta_c_c", &theta_c_c);
  reader_hf_contamination_morethantwo->AddVariable("theta_p_had_w", &theta_p_had_w);
  reader_hf_contamination_morethantwo->AddVariable("theta_w_u_b", &theta_w_u_b);
  reader_hf_contamination_morethantwo->AddVariable("theta_w_d_b", &theta_w_d_b);
  reader_hf_contamination_morethantwo->AddVariable("w_u_b_bscore", &w_u_b_bscore);
  reader_hf_contamination_morethantwo->AddVariable("w_d_b_bscore", &w_d_b_bscore);
  reader_hf_contamination_morethantwo->AddVariable("m_w_u_b", &m_w_u_b);
  reader_hf_contamination_morethantwo->AddVariable("m_w_d_b", &m_w_d_b);

  weight_file = weight_file_base + "MoreThanTwo/weights/TMVAClassification_BDTG.weights.xml";

  cout << weight_file << endl;
  reader_hf_contamination_morethantwo->BookMVA("MoreThanTwo", weight_file);

  return;
} // void Vcb::Set_Reader_HF_Contamination()

//////////

void Vcb::Set_Reader_Swapper()
{
  for (int i = 0; i < 2; i++)
  {
    reader_swapper[i] = new TMVA::Reader("!Color:!Silent");

    reader_swapper[i]->AddSpectator("n_jets", &n_sel_jet);

    reader_swapper[i]->AddVariable("met_pt", &met_pt);
    reader_swapper[i]->AddVariable("neutrino_p", &neutrino_p);
    reader_swapper[i]->AddVariable("lepton_pt", &lepton_pt);
    reader_swapper[i]->AddVariable("pt_ratio", &pt_ratio);

    reader_swapper[i]->AddVariable("pt_had_t_b", &pt_had_t_b);
    reader_swapper[i]->AddVariable("pt_w_u", &pt_w_u);
    reader_swapper[i]->AddVariable("pt_w_d", &pt_w_d);
    reader_swapper[i]->AddVariable("pt_lep_t_b", &pt_lep_t_b);

    reader_swapper[i]->AddVariable("bvsc_had_t_b", &bvsc_had_t_b);
    reader_swapper[i]->AddVariable("cvsb_had_t_b", &cvsb_had_t_b);
    reader_swapper[i]->AddVariable("cvsl_had_t_b", &cvsl_had_t_b);

    reader_swapper[i]->AddVariable("bvsc_w_u", &bvsc_w_u);
    reader_swapper[i]->AddVariable("cvsb_w_u", &cvsb_w_u);
    reader_swapper[i]->AddVariable("cvsl_w_u", &cvsl_w_u);

    reader_swapper[i]->AddVariable("bvsc_w_d", &bvsc_w_d);
    reader_swapper[i]->AddVariable("cvsb_w_d", &cvsb_w_d);
    reader_swapper[i]->AddVariable("cvsl_w_d", &cvsl_w_d);

    reader_swapper[i]->AddVariable("bvsc_lep_t_b", &bvsc_lep_t_b);
    reader_swapper[i]->AddVariable("cvsb_lep_t_b", &cvsb_lep_t_b);
    reader_swapper[i]->AddVariable("cvsl_lep_t_b", &cvsl_lep_t_b);

    reader_swapper[i]->AddVariable("theta_w_u_w_d", &theta_w_u_w_d);
    reader_swapper[i]->AddVariable("theta_had_w_had_t_b", &theta_had_w_had_t_b);
    reader_swapper[i]->AddVariable("theta_lep_neu", &theta_lep_neu);
    reader_swapper[i]->AddVariable("theta_lep_w_lep_t_b", &theta_lep_w_lep_t_b);
    reader_swapper[i]->AddVariable("del_phi_had_t_lep_t", &del_phi_had_t_lep_t);

    reader_swapper[i]->AddVariable("had_t_mass", &had_t_mass);
    reader_swapper[i]->AddVariable("had_w_mass", &had_w_mass);
    reader_swapper[i]->AddVariable("lep_t_mass", &lep_t_mass);
    reader_swapper[i]->AddVariable("lep_t_partial_mass", &lep_t_partial_mass);

    reader_swapper[i]->AddVariable("chi2_jet_had_t_b", &chi2_jet_had_t_b);
    reader_swapper[i]->AddVariable("chi2_jet_w_u", &chi2_jet_w_u);
    reader_swapper[i]->AddVariable("chi2_jet_w_d", &chi2_jet_w_d);
    reader_swapper[i]->AddVariable("chi2_jet_lep_t_b", &chi2_jet_lep_t_b);

    if (i == 0)
      reader_swapper[i]->AddSpectator("chi2_jet_extra", &chi2_jet_extra);
    else
      reader_swapper[i]->AddVariable("chi2_jet_extra", &chi2_jet_extra);

    reader_swapper[i]->AddVariable("chi2_constraint_had_t", &chi2_constraint_had_t);
    reader_swapper[i]->AddVariable("chi2_constraint_had_w", &chi2_constraint_had_w);
    reader_swapper[i]->AddVariable("chi2_constraint_lep_t", &chi2_constraint_lep_t);
    reader_swapper[i]->AddVariable("chi2_constraint_lep_w", &chi2_constraint_lep_w);
  }

  for (int i = 0; i < 4; i++)
  {
    TString weight_file = getenv("SKFlat_WD");
    weight_file += "/data/Run2UltraLegacy_v3/";
    weight_file += DataEra;
    weight_file += "/Swapper/";
    if (run_mu_ch)
      weight_file += "Mu/";
    else if (run_el_ch)
      weight_file += "El/";
    weight_file += "Permutation_";
    weight_file += to_string(i + 4);
    weight_file += "Jets/weights/TMVAClassification_BDTG.weights.xml";

    TString mva_name = "BDTG_" + to_string(i + 4) + "Jets";
    if (i == 0)
      reader_swapper[0]->BookMVA(mva_name, weight_file);
    else
      reader_swapper[1]->BookMVA(mva_name, weight_file);
  }

  return;
} // void Vcb::Set_Reader_Swapper()

//////////

void Vcb::Set_Region()
{
  bool chk_b_tagged_had_t_b = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < bvsc_had_t_b ? true : false;
  bool chk_b_tagged_lep_t_b = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < bvsc_lep_t_b ? true : false;
  bool chk_best_mva_score = 0.9 < best_mva_score ? true : false;
  bool chk_n_b_tagged = 3 <= n_b_jet ? true : false;
  bool chk_n_c_tagged = n_c_jet == 1 ? true : false;
  bool chk_c_tagged_w_u = (mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet_CvsB, JetTagging::Medium) < cvsb_w_u && mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet_CvsL, JetTagging::Medium) < cvsl_w_u) ? true : false;
  bool chk_b_tagged_w_d = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium) < bvsc_w_d ? true : false;

  // Control 0 : b-inversion
  if (chk_b_tagged_had_t_b && chk_b_tagged_lep_t_b &&
      chk_best_mva_score &&
      chk_n_c_tagged &&
      chk_c_tagged_w_u &&
      !chk_n_b_tagged &&
      !chk_b_tagged_w_d)
    region = "Control0";

  // Control 1 : c-inversion
  else if (chk_b_tagged_had_t_b && chk_b_tagged_lep_t_b &&
           chk_best_mva_score &&
           !chk_n_c_tagged &&
           !chk_c_tagged_w_u &&
           chk_n_b_tagged &&
           chk_b_tagged_w_d)
    region = "Control1";

  // Control 2 : score inversion
  else if (chk_b_tagged_had_t_b && chk_b_tagged_lep_t_b &&
           !chk_best_mva_score &&
           chk_n_c_tagged &&
           chk_c_tagged_w_u &&
           chk_n_b_tagged &&
           chk_b_tagged_w_d)
    region = "Control2";

  else
    region = "Other";

  return;
} // void Vcb::Set_Region()

//////////

void Vcb::Set_Result_Tree()
{
  for (unsigned int i = 0; i < vec_syst_type.size(); i++)
  {
    AnalyzerParameter::Syst syst_type = vec_syst_type.at(i);

    cout << "Set_Result_Tree: vec_syst_type.size() = " << vec_syst_type.size() << endl;

    TTree *result_tree = new TTree("Result_Tree", "Result_Tree");

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
      result_tree->Branch("weight_sl_trig_down", &weight_sl_trig_down);
      result_tree->Branch("weight_sl_trig_up", &weight_sl_trig_up);
    }

    result_tree->Branch("weight_top_pt", &weight_top_pt);

    result_tree->Branch("n_vertex", &nPV);

    result_tree->Branch("lepton_pt", &lepton_pt);
    result_tree->Branch("lepton_eta", &lepton_eta);

    result_tree->Branch("n_jets", &n_sel_jet);
    result_tree->Branch("n_bjets", &n_b_jet);
    result_tree->Branch("n_cjets", &n_c_jet);

    result_tree->Branch("pt_leading_jet", &pt_leading_jet);
    result_tree->Branch("pt_subleading_jet", &pt_subleading_jet);

    result_tree->Branch("eta_leading_jet", &eta_leading_jet);
    result_tree->Branch("eta_subleading_jet", &eta_subleading_jet);

    result_tree->Branch("bvsc_leading_jet", &bvsc_leading_jet);
    result_tree->Branch("cvsb_leading_jet", &cvsb_leading_jet);
    result_tree->Branch("cvsl_leading_jet", &cvsl_leading_jet);

    result_tree->Branch("bvsc_subleading_jet", &bvsc_subleading_jet);
    result_tree->Branch("cvsb_subleading_jet", &cvsb_subleading_jet);
    result_tree->Branch("cvsl_subleading_jet", &cvsl_subleading_jet);

    result_tree->Branch("met_pt", &met_pt);
    result_tree->Branch("met_phi", &met_phi);

    result_tree->Branch("best_mva_score_pre", &best_mva_score_pre);
    result_tree->Branch("best_mva_score", &best_mva_score);
    result_tree->Branch("best_chi2", &best_chi2);
    result_tree->Branch("ht", &ht);
    result_tree->Branch("mt", &mt);
    result_tree->Branch("mva_hf_score", &mva_hf_score);

    result_tree->Branch("bvsc_had_t_b", &bvsc_had_t_b);
    result_tree->Branch("cvsb_had_t_b", &cvsb_had_t_b);
    result_tree->Branch("cvsl_had_t_b", &cvsl_had_t_b);
    result_tree->Branch("bvsc_w_u", &bvsc_w_u);
    result_tree->Branch("cvsb_w_u", &cvsb_w_u);
    result_tree->Branch("cvsl_w_u", &cvsl_w_u);
    result_tree->Branch("bvsc_w_d", &bvsc_w_d);
    result_tree->Branch("cvsb_w_d", &cvsb_w_d);
    result_tree->Branch("cvsl_w_d", &cvsl_w_d);
    result_tree->Branch("bvsc_lep_t_b", &bvsc_lep_t_b);
    result_tree->Branch("cvsb_lep_t_b", &cvsb_lep_t_b);
    result_tree->Branch("cvsl_lep_t_b", &cvsl_lep_t_b);

    result_tree->Branch("pt_had_t_b", &pt_had_t_b);
    result_tree->Branch("pt_w_u", &pt_w_u);
    result_tree->Branch("pt_w_d", &pt_w_d);
    result_tree->Branch("pt_lep_t_b", &pt_lep_t_b);

    result_tree->Branch("eta_had_t_b", &eta_had_t_b);
    result_tree->Branch("eta_w_u", &eta_w_u);
    result_tree->Branch("eta_w_d", &eta_w_d);
    result_tree->Branch("eta_lep_t_b", &eta_lep_t_b);

    result_tree->Branch("m_had_t", &m_had_t);
    result_tree->Branch("m_had_w", &m_had_w);
    result_tree->Branch("m_lep_t", &m_lep_t);
    result_tree->Branch("m_lep_w", &m_lep_w);

    result_tree->Branch("m_w_u", &m_w_u);
    result_tree->Branch("m_w_d", &m_w_d);

    result_tree->Branch("swapped_mva", &swapped_mva);

    // For MC
    result_tree->Branch("decay_mode", &decay_mode);
    result_tree->Branch("chk_reco_correct", &chk_reco_correct);
    result_tree->Branch("chk_included", &chk_included);
    result_tree->Branch("chk_gentau_conta", &chk_gentau_conta);
    result_tree->Branch("chk_hf_contamination", &chk_hf_contamination);

    result_tree->Branch("pu_conta_had_t_b", &pu_conta_had_t_b);
    result_tree->Branch("pu_conta_w_u", &pu_conta_w_u);
    result_tree->Branch("pu_conta_w_d", &pu_conta_w_d);
    result_tree->Branch("pu_conta_lep_t_b", &pu_conta_lep_t_b);

    result_tree->Branch("swapped_truth", &swapped_truth);

    result_tree->Branch("Gen_HF_Flavour", &vec_gen_hf_flavour);
    result_tree->Branch("Gen_HF_Origin", &vec_gen_hf_origin);

    result_tree->Branch("Sel_Gen_HF_Flavour", &vec_sel_gen_hf_flavour);
    result_tree->Branch("Sel_Gen_HF_Origin", &vec_sel_gen_hf_origin);

    map_result_tree.insert({syst_type, result_tree});
  }

  return;
} // void Vcb::Set_Result_Tree()

//////////

void Vcb::Set_Template_Truth_Tree()
{
  template_truth_tree[0] = new TTree("Template_Truth_21", "Template_Truth_21");
  template_truth_tree[1] = new TTree("Template_Truth_23", "Template_Truth_23");
  template_truth_tree[2] = new TTree("Template_Truth_41", "Template_Truth_41");
  template_truth_tree[3] = new TTree("Template_Truth_43", "Template_Truth_43");
  template_truth_tree[4] = new TTree("Template_Truth_45", "Template_Truth_45");

  for (int i = 0; i < 5; i++)
  {
    template_truth_tree[i]->Branch("best_mva_score_pre", &best_mva_score_pre); // dummy
    template_truth_tree[i]->Branch("bvsc_w_u", &bvsc_w_u);
    template_truth_tree[i]->Branch("cvsb_w_u", &cvsb_w_u);
    template_truth_tree[i]->Branch("cvsl_w_u", &cvsl_w_u);
    template_truth_tree[i]->Branch("m_w_u", &m_w_u);
    template_truth_tree[i]->Branch("bvsc_w_d", &bvsc_w_d);
    template_truth_tree[i]->Branch("cvsb_w_d", &cvsb_w_d);
    template_truth_tree[i]->Branch("cvsl_w_d", &cvsl_w_d);
    template_truth_tree[i]->Branch("m_w_d", &m_w_d);
  }

  return;
} // void Vcb::Set_Template_Truth_Tree()

//////////

void Vcb::Sol_Neutrino_Pz(const Particle &lepton, const Particle &met, float neutrino_pz_sol[2])
{
  float lepton_mass = lepton.M();

  float met_px = met.Px();
  float met_py = met.Py();
  float met_pt = met.Pt();

  double k = TMath::Power(W_MASS, 2.) / 2.0 - lepton_mass * lepton_mass / 2.0 + lepton.Px() * met_px + lepton.Py() * met_py;
  double a = TMath::Power(lepton.Px(), 2.0) + TMath::Power(lepton.Py(), 2.0);
  double b = -2 * k * lepton.Pz();
  double c = TMath::Power(lepton.Pt(), 2.0) * TMath::Power(met_pt, 2.0) - TMath::Power(k, 2.0);

  double determinant = TMath::Power(b, 2.0) - 4 * a * c;

  // real solution
  if (determinant >= 0)
  {
    neutrino_pz_sol[0] = (-b + TMath::Sqrt(determinant)) / (2 * a);
    neutrino_pz_sol[1] = (-b - TMath::Sqrt(determinant)) / (2 * a);
  }
  // complex solution
  else
  {
    neutrino_pz_sol[0] = -b / (2 * a);
    // Resol_Neutrino_Pt();
  }

  return;
} // void Vcb::Sol_Neutrino_Pz(const Particle& lepton, const Particle& met)

//////////
