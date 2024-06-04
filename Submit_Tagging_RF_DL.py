#!/usr/bin/env python3

import argparse
import os

##########

parser = argparse.ArgumentParser(description='SKFlat -a Vcb Command')
parser.add_argument('-e', dest='era', default="2018")
#parser.add_argument('-mc', action='store_true', default="")
parser.add_argument('-nmax', dest='nmax', default="300")
args = parser.parse_args()

if args.era=="2016a": args.era="2016preVFP"
if args.era=="2016b": args.era="2016postVFP"

# key:[njob for SkimTree, njob Vcb Analyzer]
#data_list = {"SingleMuon":[300,100], "EGamma":[300,50]}

mc_list = {
    "TTLJ_WtoCB_powheg":[40,20],
    "TTLJ_powheg":[300,200], 
    "TTLL_powheg":[300,100], 
    #"TTBB":[200,100], 
    "SingleTop_sch_Lep":[20,40],
    "SingleTop_tch_antitop_Incl":[20,40],
    "SingleTop_tch_top_Incl":[20,40],
    "SingleTop_tW_antitop_NoFullyHad":[30,60], 
    "SingleTop_tW_top_NoFullyHad":[30,60], 
    "DYJets_MG":[50,40], 
    "WJets_MG":[50,40],
    "DYJets":[50,80],
    "WJets_Sherpa":[50,80],
    "QCD_bEnriched_HT100to200": [10,20],
    "QCD_bEnriched_HT200to300": [10,20],
    "QCD_bEnriched_HT300to500": [10,15],
    "QCD_bEnriched_HT500to700": [10,10],
    "QCD_bEnriched_HT700to1000": [10,10],
    "QCD_bEnriched_HT1000to1500": [10,10],
    "QCD_bEnriched_HT1500to2000": [10,10],
    "QCD_bEnriched_HT2000toInf": [10,10],
    "TTLJ_powheg_CP5Down":[300,200],
    "TTLJ_powheg_CP5Up":[300,200],
    "TTLJ_powheg_hdampDown":[300,200],
    "TTLJ_powheg_hdampUp":[300,200],
    "TTLJ_powheg_mtop171p5":[300,200],
    "TTLJ_powheg_mtop173p5":[300,200],
    "TTLL_powheg_CP5Down":[300,100],
    "TTLL_powheg_CP5Up":[300,100],
    "TTLL_powheg_hdampDown":[300,100],
    "TTLL_powheg_hdampUp":[300,100],
    "TTLL_powheg_mtop171p5":[300,100],
    "TTLL_powheg_mtop173p5":[300,100],
    "ttWToLNu": [20, 40],
    "ttWToQQ": [20, 10],
    "ttZToLLNuNu": [20, 60],
    "ttZToQQ": [20, 60],
    "ttHToNonbb": [20, 40],
    "ttHTobb": [20, 60],
    "WW_pythia": [20, 10],
    "WZ_pythia": [20, 10],
    "ZZ_pythia": [20, 10],
}


## Vcb_Tagging_RF
# data
#if args.data == True:
#    data  = "SingleMuon"
#        
#    operation = f"nohup SKFlat.py -a Vcb_Tagging_RF_DL -i {data} -n {data_list[data][1]} -e {args.era} --userflag RunMu &"
#      
#    print(operation)
#    os.system(operation)
    
# mc
for mc in mc_list:
    operation =  f"nohup SKFlat.py -a Vcb_Tagging_RF_DL -i {mc} -n {mc_list[mc][1]} -e {args.era} --nmax {args.nmax} &"
        
    print(operation)
    os.system(operation)
    
