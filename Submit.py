#!/usr/bin/env python3

import argparse
import os

##########

parser = argparse.ArgumentParser(description='SKFlat -a Vcb Command')
parser.add_argument('-e', dest='Era', default="")
parser.add_argument('-ch', dest='Channel', default="")
parser.add_argument('-data', action='store_true', default="")
parser.add_argument('-mc', action='store_true', default="")
args = parser.parse_args()

if args.Era=="2016a": args.Era="2016preVFP"
if args.Era=="2016b": args.Era="2016postVFP"

analyzer = ""
if args.Channel == "Mu":
    analyzer = "Vcb_Mu"
elif args.Channel == "El":
    analyzer = "Vcb_El"
elif args.Channel == "Skim":
    analyzer = "SkimTree_Vcb"
else:
    print("Wrong argument -ch")

print(f"Selected analyzer is {analyzer}")

# key:[njob for SkimTree, njob Vcb Analyzer]
data_list = {"SingleMuon":[300,20], "EGamma":[300,20]}
mc_list = {"TTLJ_powheg":[300,200], "TTLL_powheg":[300,100], "TTBB":[200,100], "SingleTop_tW_antitop_NoFullyHad":[100,20], "SingleTop_tW_top_NoFullyHad":[100,20], "DYJets_MG":[20,5], "WJets_MG":[20,1], "TTLJ_WtoCB_powheg":[40:20] }

operation = ""

## SkimTree
if args.Channel == "Skim":
    # data
    for data in data_list:
        operation = f"nohup SKFlat.py -a SkimTree_Vcb -i {data} -n {data_list[data][0]} -e {args.Era} &"
                
        if args.data == True:
            print(operation)
            os.system(operation)
    
    # mc
    for mc in mc_list:
        operation = f"nohup SKFlat.py -a SkimTree_Vcb -i {mc} -n {mc_list[mc][0]} -e {args.Era} &"
        
        if args.mc == True: 
            print(operation)
            os.system(operation)

## Vcb_El or Vcb_Mu
if args.Channel == "El" or args.Channel == "Mu":
    # data
    data  = ""
    if args.Channel == "El": data = "EGamma"
    elif args.Channel == "Mu": data = "SingleMuon"
        
    operation = f"nohup SKFlat.py -a Vcb_{args.Channel} -i SkimTree_Vcb_{data} -n {data_list[data][1]} -e {args.Era} &"
      
    if args.data == True:
        print(operation)
        os.system(operation)
    
    # mc
    for mc in mc_list:
        operation =  f"nohup SKFlat.py -a Vcb_{args.Channel} -i {mc}_Skim -n {mc_list[mc][1]} -e {args.Era} &"
        
        if args.mc == True: 
            print(operation)
            os.system(operation)
    
