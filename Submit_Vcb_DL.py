#!/usr/bin/env python3

import argparse
import os

##########

parser = argparse.ArgumentParser(description='SKFlat -a Vcb Command')
parser.add_argument('-e', dest='era', default="2018")
parser.add_argument('-data', action='store_true', default="")
parser.add_argument('-mc', action='store_true', default="")
parser.add_argument('-flag', dest='flag', default="All")
parser.add_argument('-nmax', dest='nmax', default="400")
args = parser.parse_args()

if args.era=="2016a": args.era="2016preVFP"
if args.era=="2016b": args.era="2016postVFP"

# key:[njob for SkimTree, njob Vcb Analyzer]

# data
if args.data == True:
    data_list = {"SingleMuon":[40,40], "EGamma":[40,40]}
    if args.era != "2018":
        del data_list["EGamma"]
        data_list["SingleElectron"] = [40, 40]
    
    for data in data_list:
        operation = f"nohup SKFlat.py -a Vcb_DL -i {data} -n {data_list[data][1]} -e {args.era} --nmax {args.nmax} &"
            
        print(operation)
        os.system(operation)
    
# mc
if args.mc == True:
    if args.flag == "Syst" or args.flag == "All":
        mc_list = {"TTLJ_WtoCB_powheg":[40,30],
                   "TTLJ_powheg":[300,300],
                   "TTLL_powheg":[300,100],
                   #"TTBB":[200,100],
                   "SingleTop_sch_Lep":[100,30],
                   "SingleTop_tW_antitop_NoFullyHad":[100,30],
                   "SingleTop_tW_top_NoFullyHad":[100,30],
                   "SingleTop_tch_antitop_Incl":[100,50],
                   "SingleTop_tch_top_Incl":[100,50],
                   "DYJets_MG":[20,40],
                   "WJets_MG":[20,40],
                   "DYJets":[20,50],
                   "WJets_Sherpa":[20,50],
                   "QCD_bEnriched_HT100to200": [3,10],
                   "QCD_bEnriched_HT200to300": [2,10],
                   "QCD_bEnriched_HT300to500": [2,5],
                   "QCD_bEnriched_HT500to700": [2,2],
                   "QCD_bEnriched_HT700to1000": [2,2],
                   "QCD_bEnriched_HT1000to1500": [2,2],
                   "QCD_bEnriched_HT1500to2000": [2,2],
                   "QCD_bEnriched_HT2000toInf": [2,2],
                   "ttWToLNu": [20, 20],
                   "ttWToQQ": [20, 10],
                   "ttZToLLNuNu": [20, 30],
                   "ttZToQQ": [20, 30],
                   "ttHToNonbb": [20, 20],
                   "ttHTobb": [20, 30],
                   "WW_pythia": [20, 20],
                   "WZ_pythia": [20, 20],
                   "ZZ_pythia": [20, 20],
        }

        for mc in mc_list:
            operation =  f"nohup SKFlat.py -a Vcb_DL -i {mc} -n {mc_list[mc][1]} -e {args.era} --userflag RunSyst --nmax {args.nmax} &"
            print(operation)
            os.system(operation)

    if args.flag == "Syst_Top" or args.flag == "All":
        mc_list = {"TTLJ_powheg_CP5Down": [300, 200],
                   "TTLJ_powheg_CP5Up": [300, 200],
                   "TTLJ_powheg_hdampDown": [300, 200],
                   "TTLJ_powheg_hdampUp": [300, 200],
                   "TTLJ_powheg_mtop171p5": [300, 200],
                   "TTLJ_powheg_mtop173p5": [300, 200],
                   "TTLL_powheg_CP5Down": [300, 100],
                   "TTLL_powheg_CP5Up": [300, 100],
                   "TTLL_powheg_hdampDown": [300, 100],
                   "TTLL_powheg_hdampUp": [300, 100],
                   "TTLL_powheg_mtop171p5": [300, 100],
                   "TTLL_powheg_mtop173p5": [300, 100],
               }
        
        for mc in mc_list:
            operation = f"nohup SKFlat.py -a Vcb_DL -i {mc} -n {mc_list[mc][1]} -e {args.era} --nmax {args.nmax} &"
            print(operation)
            os.system(operation)
    
