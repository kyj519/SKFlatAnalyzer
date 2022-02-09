#ifndef Vcb_Def_h
#define Vcb_Def_h

const float MUON_PT = 20.;
const float MUON_PT_VETO = 15.;
const float MUON_ETA = 2.4;

const float ELECTRON_PT = 15.;
const float ELECTRON_PT_VETO = 15;
const float ELECTRON_ETA = 2.5;

const float JET_PT = 30.;
const float JET_ETA = 2.4;
const float DR_LEPTON_VETO = 0.3;

const float MET_PT = 20.;

//selection variables for making SkimTree.
//to veto events with at least two tight leptons
const float MUON_PT_SKIM = 15;
const float MUON_ETA_SKIM = 2.4;

const float ELECTRON_PT_SKIM = 15.;
const float ELECTRON_ETA_SKIM = 2.5;

const float JET_PT_SKIM = 30.;
const float JET_ETA_SKIM = 2.4;
const float DR_LEPTON_VETO_SKIM = 0.3;

const float MET_PT_SKIM = 20.;

#endif /* Vcb_Def_h */
