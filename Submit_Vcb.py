#!/usr/bin/env python3

import argparse
import os

##########

parser = argparse.ArgumentParser(description='SKFlat -a Vcb Command')
parser.add_argument('-e', dest='era', default="2018")
parser.add_argument('-ch', dest='channel', default="Mu")
parser.add_argument('-data', action='store_true', default="")
parser.add_argument('-mc', action='store_true', default="")
parser.add_argument('-flag', dest='flag', default="Syst")
args = parser.parse_args()

if args.era=="2016a": args.era="2016preVFP"
if args.era=="2016b": args.era="2016postVFP"

if args.channel == "Mu":
    if args.flag == "EecDown" or args.flag == "EecUp" or args.flag == "EerDown" or args.flag == "EerUp":
        print("Systematic variations for electron are not compatible with Muon channel")
        exit()

# key:[njob for SkimTree, njob Vcb Analyzer]
data_list = {"SingleMuon":[300,50], "EGamma":[300,50]}
if args.era != "2018":
    del data_list["EGamma"]
    data_list["SingleElectron"] = [300, 50]

mc_list = {
    "TTLJ_WtoCB_powheg":[40,20],
    "TTLJ_powheg":[300,300], 
    "TTLL_powheg":[300,100], 
    #"TTBB":[200,100], 
    "SingleTop_sch_Lep":[100,10],
    "SingleTop_tW_antitop_NoFullyHad":[100,10], 
    "SingleTop_tW_top_NoFullyHad":[100,10],
    "SingleTop_tch_antitop_Incl":[100,20],
    "SingleTop_tch_top_Incl":[100,20],
    "DYJets_MG":[20,10], 
    "WJets_MG":[20,5],
    "QCD_bEnriched_HT100to200": [3,3],
    "QCD_bEnriched_HT200to300": [2,2],
    "QCD_bEnriched_HT300to500": [2,1],
    "QCD_bEnriched_HT500to700": [2,1],
    "QCD_bEnriched_HT700to1000": [2,1],
    "QCD_bEnriched_HT1000to1500": [2,1],
    "QCD_bEnriched_HT1500to2000": [2,1],
    "QCD_bEnriched_HT2000toInf": [2,1],
    "ttWToLNu": [20, 20],
    "ttWToQQ": [20, 5],
    "ttZToLLNuNu": [20, 30],
    "ttZToQQ": [20, 30],
    "ttHToNonbb": [20, 20],
    "ttHTobb": [30, 30],
    "WW_pythia": [20, 5],
    "WZ_pythia": [20, 5],
    "ZZ_pythia": [20, 5],
}

operation = ""

## SkimTree
# obsolete
if args.channel == "Skim":
    # data
    for data in data_list:
        operation = f"nohup SKFlat.py -a SkimTree_Vcb -i {data} -n {data_list[data][0]} -e {args.era} &"
                
        if args.data == True:
            print(operation)
            os.system(operation)
    
    # mc
    for mc in mc_list:
        operation = f"nohup SKFlat.py -a SkimTree_Vcb -i {mc} -n {mc_list[mc][0]} -e {args.era} &"
        
        if args.mc == True: 
            print(operation)
            os.system(operation)

## Vcb
if args.channel == "El" or args.channel == "Mu":
    # data
    if args.data == True:
        data  = ""
        
        if args.channel == "El":
            if args.era == "2018": data = "EGamma"
            else: data = "SingleElectron"
        elif args.channel == "Mu": data = "SingleMuon"
        
        operation = f"nohup SKFlat.py -a Vcb -i {data} -n {data_list[data][1]} -e {args.era} --userflag Run{args.channel},RunResult &"
      
        print(operation)
        os.system(operation)
    
    # mc
    if args.mc == True:
        if args.flag == "Syst_Top":
            #del mc_list["TTLJ_powheg"] 
            mc_list = dict()
            mc_list["TTLJ_powheg_CP5Down"] = [300, 100]
            mc_list["TTLJ_powheg_CP5Up"] = [300, 100]
            mc_list["TTLJ_powheg_hdampDown"] = [300, 100]
            mc_list["TTLJ_powheg_hdampUp"] = [300, 100]
            mc_list["TTLJ_powheg_mtop171p5"] = [300, 100]
            mc_list["TTLJ_powheg_mtop173p5"] = [300, 100]
            mc_list["TTLL_powheg_CP5Down"] = [300, 30]
            mc_list["TTLL_powheg_CP5Up"] = [300, 30]
            mc_list["TTLL_powheg_hdampDown"] = [300, 30]
            mc_list["TTLL_powheg_hdampUp"] = [300, 30]
            mc_list["TTLL_powheg_mtop171p5"] = [300, 30]
            mc_list["TTLL_powheg_mtop173p5"] = [300, 30]

        for mc in mc_list:
            if args.flag == "Result":
                operation =  f"nohup SKFlat.py -a Vcb -i {mc} -n {mc_list[mc][1]} -e {args.era} --userflag Run{args.channel},RunResult &"
            elif args.flag == "Syst":
                operation =  f"nohup SKFlat.py -a Vcb -i {mc} -n {mc_list[mc][1]} -e {args.era} --userflag Run{args.channel},RunResult,RunSyst &"
            elif args.flag == "Syst_Top":
                operation = f"nohup SKFlat.py -a Vcb -i {mc} -n {mc_list[mc][1]} -e {args.era} --userflag Run{args.channel},RunResult &"
            
            print(operation)
            os.system(operation)
    
