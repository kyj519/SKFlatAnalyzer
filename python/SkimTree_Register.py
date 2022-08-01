#!/usr/bin/env python3

import argparse
import os
import shutil
from pathlib import Path

##########

parser = argparse.ArgumentParser(description='SkimTree Register Command')
parser.add_argument('-e', dest='Era', default="")
parser.add_argument('-s', dest='Sample', default="")
args = parser.parse_args()

if args.Era=="2016a": args.Era="2016preVFP"
if args.Era=="2016b": args.Era="2016postVFP"

AvailableDataPeriods = []
if args.Era == "2016preVFP":
    AvailableDataPeriods = ["B_ver2","C","D","E","F"]
elif args.Era == "2016postVFP":
    AvailableDataPeriods = ["F","G","H"]
elif args.Era == "2017":
    AvailableDataPeriods = ["B","C","D","E","F"]
elif args.Era == "2018":
    AvailableDataPeriods = ["A", "B","C","D"]
else:
    print("Wrong Era : "+args.Era)
    exit()

chk_data = True
data_sets = []
if (args.Sample.find("EGamma")!=-1) or (args.Sample.find("SingleMuon")!=-1):
    chk_data = True
    print("Data")
    
    if(args.Sample.find(":")!=-1):
        data_sets = [args.Sample]
        data = args.Sample.split(":")[0]
    else:
        print("All period")
        data = args.Sample
        for period in AvailableDataPeriods:
            data_sets.append(f"{args.Sample}:{period}")
else:
    chk_data = False
    print("MC")

SKFlat_WD = os.environ['SKFlat_WD']
txt_path = f"{SKFlat_WD}/data/Run2UltraLegacy_v2/{args.Era}/Sample"

##########

## edit CommonSampleInfo.txt

if chk_data:
    pass
else:
    path_common_f = f"{txt_path}/CommonSampleInfo/{args.Sample}.txt"
    path_common_copy = f"{txt_path}/CommonSampleInfo/{args.Sample}_Skim.txt"
    
    common_f = open(path_common_f, "rt")
    common_copy_f = open(path_common_copy, "wt")
    
    for line in common_f:
        common_copy_f.write(line.replace(args.Sample, f"{args.Sample}_Skim"))
        
        if '#' in line: 
            continue
        else:
            common_val = line.split()
        
    full_name = common_val[1]
    
    
## root file list

path_list = []
if chk_data:
    for data_set in data_sets:
        data_period = data_set.split(":")
        path_list.append(f"/gv0/Users/isyoon/Run2UltraLegacy_v2/{args.Era}/DATA_SkimTree_Vcb/{data_period[0]}/period{data_period[1]}")
else:
    path_list.append(f"/gv0/Users/isyoon/Run2UltraLegacy_v2/{args.Era}/MC_SkimTree_Vcb/{full_name}")

for path in path_list:
    #choose newest folder
    folder_list = sorted(Path(path).iterdir(), key=os.path.getmtime)
        
    root_list = os.listdir(folder_list[-1])
    root_list = [root_f for root_f in root_list if root_f.endswith(".root")]

    if chk_data:
        period = path[-1]
        path_for_snu = f"{txt_path}/ForSNU/SkimTree_Vcb_{data}_{period}.txt"
    else:
        path_for_snu = f"{txt_path}/ForSNU/{args.Sample}_Skim.txt"
        
    root_list_f = open(path_for_snu, "w")
    
    for root_f in root_list:
        root_list_f.write(f"{folder_list[-1]}/{root_f}\n")

## update SampleSummary_*.txt

if chk_data:
    path_summary_f = f"{txt_path}/SampleSummary_DATA.txt"

    summary_fin = open(path_summary_f, "rt")
    lines = summary_fin.readlines()
    summary_fin.close()

    vals = []
    for line in lines:
        vals.append(line.split())
        
    for data_set in data_sets:
                
        chk_update = True
        for val in vals:
            period = data_set.split(":")[1]
            if f"SkimTree_Vcb_{data}" == val[0] and period == val[1]:
                chk_update = False
            
        if chk_update:
            
            for val in vals:
                period = data_set.split(":")[1]
                if data == val[0] and period == val[1]: 
                    n_event = val[2]

            path_summary_fout = open(path_summary_f, "a")
            path_summary_fout.write(f"SkimTree_Vcb_{data}\t{period}\t{n_event}\n")
            path_summary_fout.close()
        
else:
    path_summary_f = f"{txt_path}/SampleSummary_MC.txt"

    summary_fin = open(path_summary_f, "rt")
    lines = summary_fin.readlines() 
    summary_fin.close()
    
    chk = False
    for line in lines:
        val = line.split()
        
        if(val[0] == f"{args.Sample}_Skim"):
            chk = True
            
    #update
    if(chk==False):
        path_summary_fout = open(path_summary_f, "a")
        path_summary_fout.write(f"{args.Sample}_Skim\t{common_val[1]}\t{common_val[2]}\t{common_val[3]}\t{common_val[4]}\n")
