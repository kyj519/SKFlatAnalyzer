#ifndef Vcb_Def_h
#define Vcb_Def_h

const float T_MASS = 172.5;
const float T_WIDTH = 1.5;
const float W_MASS = 80.379;
const float W_WIDTH = 2.085;

const float Z_MASS = 91.1876;

const float MUON_PT = 20.;
const float MUON_PT_VETO = 15.;
const float MUON_ETA = 2.4;

const float ELECTRON_PT = 15.;
const float ELECTRON_PT_VETO = 15;
const float ELECTRON_ETA = 2.5;

const float JET_PT = 20.;
const float JET_ETA = 2.4;
const float JET_ETA_2016 = 2.4;
const float DR_LEPTON_VETO = 0.3;

const float JET_PT_MATCH = 10;
const float JET_ETA_MATCH = 10;

const float MET_PT = 20.;
const float MET_PT_DL = 40.;

//selection variables for making SkimTree.
//to veto events with at least two tight leptons
const float MUON_PT_SKIM = 15;
const float MUON_ETA_SKIM = 2.4;

const float ELECTRON_PT_SKIM = 15.;
const float ELECTRON_ETA_SKIM = 2.5;

const float JET_PT_SKIM = 10.;
const float JET_ETA_SKIM = 10;
const float DR_LEPTON_VETO_SKIM = 0.3;

const float MET_PT_SKIM = 20.;

const float cvsl_2016a_m = 0.098;
const float cvsb_2016a_m = 0.370; 

const float cvsl_2016b_m = 0.099;
const float cvsb_2016b_m = 0.353;

const float cvsl_2017_m = 0.085;
const float cvsb_2017_m = 0.34;

const float cvsl_2018_m = 0.099;
const float cvsb_2018_m = 0.325;

const float REL_ISO_MUON = 0.15;

//medium
const float REL_ISO_ELECTRON_BARREL_A = 0.0478;
const float REL_ISO_ELECTRON_BARREL_B = 0.506;
const float REL_ISO_ELECTRON_ENDCAP_A = 0.0658;
const float REL_ISO_ELECTRON_ENDCAP_B = 0.963;

#endif /* Vcb_Def_h */
