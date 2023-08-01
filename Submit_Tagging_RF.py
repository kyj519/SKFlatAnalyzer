#!/usr/bin/env python3

import argparse
import os

##########

parser = argparse.ArgumentParser(description='SKFlat -a Vcb Command')
parser.add_argument('-e', dest='era', default="2018")
parser.add_argument('-data', action='store_true', default="")
parser.add_argument('-mc', action='store_true', default="")
parser.add_argument('-ch', dest='channel', default="")
args = parser.parse_args()

if args.era=="2016a": args.era="2016preVFP"
if args.era=="2016b": args.era="2016postVFP"

# key:[njob for SkimTree, njob Vcb Analyzer]
data_list = {"SingleMuon":[300,100], "EGamma":[300,50]}
mc_list = {
    "TTLJ_WtoCB_powheg":[40,10],
    "TTLJ_powheg":[300,200], 
    "TTLL_powheg":[300,100], 
    #"TTBB":[200,100], 
    "SingleTop_sch_Lep":[20,20],
    "SingleTop_tch_antitop_Incl":[20,20],
    "SingleTop_tch_top_Incl":[20,20],
    "SingleTop_tW_antitop_NoFullyHad":[100,20], 
    "SingleTop_tW_top_NoFullyHad":[100,20], 
    "DYJets_MG":[50,50], 
    "WJets_MG":[50,50],
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
    "ttWToLNu": [20, 20],
    "ttWToQQ": [20, 5],
    "ttZToLLNuNu": [20, 30],
    "ttZToQQ": [20, 30],
    "ttHToNonbb": [20, 20],
    "ttHTobb": [20, 30],
    "WW_pythia": [20, 5],
    "WZ_pythia": [20, 5],
    "ZZ_pythia": [20, 5],
}


## Vcb_Tagging_RF
# data
if args.data == True:
    data  = "SingleMuon"
        
    operation = f"nohup SKFlat.py -a Vcb_Tagging_RF -i {data} -n {data_list[data][1]} -e {args.era} --userflag RunMu --reduction 100 &"
      
    print(operation)
    os.system(operation)
    
# mc
if args.mc == True:
    for mc in mc_list:
        operation =  f"nohup SKFlat.py -a Vcb_Tagging_RF -i {mc} -n {mc_list[mc][1]} -e {args.era} --userflag Run{args.channel} &"
        
        print(operation)
        os.system(operation)
    
