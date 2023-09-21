#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-e', dest='Era', default='2018')
parser.add_argument('-p', dest='Chk_Private', type=bool, default=False)
parser.add_argument('-s', dest='Sample', default='')
args = parser.parse_args()

if args.Era=="2016a": args.Era="2016preVFP"
if args.Era=="2016b": args.Era="2016postVFP"

mc_list = {'TTLJ_powheg':['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'],
           'TTLL_powheg':['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8'],
           'tw_antitop':['ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8'],
           'tW_top':['ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8'],
           'DYJets_MG':['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'],
           'WJet_MG':['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8'],
           'QCD_bEn_HT100to200':['QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8'],
           'QCD_bEn_HT200to300':['QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8'],
           'QCD_bEn_HT300to500':['QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8'],
           'QCD_bEn_HT500to700':['QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8'],
           'QCD_bEn_HT700to1000':['QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8'],
           'QCD_bEn_HT1000to1500':['QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8'],
           'QCD_bEn_HT1500to2000':['QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8'],
           'QCD_bEn_HT2000toInf':['QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8'],
           #'TTLJ_WtoCB_powheg':['ttj_Vcb_NLO_FXFX_Summer20UL18_LHEGEN', 'ttj_Vcb_NLO_FXFX_Summer20UL18_GEN_extra'],
           'TTLJ_WtoCB_powheg':['TTToSemiLeptonic_Vcb_TuneCP5_13TeV-powheg-pythia8'],
           'TTLJ_powheg_hdampDown':['TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8'],
           'TTLJ_powheg_hdampUp':['TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8'],
           'TTLJ_powheg_mtop171p5':['TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8'],
           'TTLJ_powheg_mtop173p5':['TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8'],
           'TTLJ_powheg_CP5Down':['TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8'],
           'TTLJ_powheg_CP5Up':['TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8'],
           'TTLL_powheg_hdampDown':['TTTo2L2Nu_hdampDOWN_TuneCP5_13TeV-powheg-pythia8'],
           'TTLL_powheg_hdampUp':['TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8'],
           'TTLL_powheg_mtop171p5':['TTTo2L2Nu_mtop171p5_TuneCP5_13TeV-powheg-pythia8'],
           'TTLL_powheg_mtop173p5':['TTTo2L2Nu_mtop173p5_TuneCP5_13TeV-powheg-pythia8'],
           'TTLL_powheg_CP5Down':['TTTo2L2Nu_TuneCP5down_13TeV-powheg-pythia8'],
           'TTLL_powheg_CP5Up':['TTTo2L2Nu_TuneCP5up_13TeV-powheg-pythia8'],
           'TTbb':['TTbb_4f_TTToSemiLeptonic_TuneCP5-Powheg-Openloops-Pythia8'],
           'ttHTobb':['ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8'],
}

#for private only sample
if args.Chk_Private == False:
    if args.Sample == 'TTLJ_WtoCB_powheg' or args.Sample == 'TTLJ_powheg_hdampDown' or args.Sample == 'TTLJ_powheg_hdampUp' or args.Sample == 'TTLJ_powheg_m171p5' or args.Sample == 'TTLJ_powheg_m173p5' or args.Sample == 'TTLJ_powheg_CP5Down' or args.Sample == 'TTLJ_powheg_CP5Up' or args.Sample == 'TTbb' or args.Sample == 'ttHtobb':
        print(f"{args.Sample} is private only sample. Check argument first. Add -p 1")
        exit(1)

#if args.Sample == 'TTLJ_WtoCB_powheg':
#    if args.Era == '2017':
#        mc_list[args.Sample] = ['ttj_Vcb_NLO_FXFX_Summer20UL17_9M_Events'];
#    elif args.Era == '2016postVFP':
#        mc_list[args.Sample] = ['ttj_Vcb_NLO_FXFX_Summer20UL16_9M_Events'];
#    elif args.Era == '2016preVFP':
#        mc_list[args.Sample] = ['ttj_Vcb_NLO_FXFX_Summer20UL16APV_9M_Events'];

target_mc_list = list()
for mc in mc_list:
    if args.Sample == mc:
        target_mc_list = mc_list[mc]
        break

if len(target_mc_list) == 0:
    print(f"Can't find {args.Sample} in mc_list")
    exit(1)

import os 
import subprocess 

if os.path.exists(f"{args.Sample}.txt"):
    os.remove(f"{args.Sample}.txt")

fout = open(f"{args.Sample}.txt", "w")

if args.Chk_Private:
    for mc in target_mc_list:
        path=f"/gv0/Users/isyoon/Run2UltraLegacy_v3/{args.Era}/{mc}/SKFlat_Run2UltraLegacy_v3/"
        
        #find lastest crab output
        result = os.listdir(path)

        #ingore folder end with Skip. It's only for archiving
        result = [folder for folder in result if not folder.endswith("Skip")]
        
        path += str(result[-1])
        
        find = subprocess.Popen(['find', path, '-type', 'f'], stdout=subprocess.PIPE)
        subprocess.run(['grep', 'root'], stdin=find.stdout, stdout=fout)
else:
    path=f"/gv0/DATA/SKFlat/Run2UltraLegacy_v3/{args.Era}/MC/{target_mc_list[0]}/"
    
    #find lastest crab output
    result = os.listdir(path)
    path += str(result[-1])

    find = subprocess.Popen(['find', path, '-type', 'f'], stdout=subprocess.PIPE)
    subprocess.run(['grep', 'root'], stdin=find.stdout, stdout=fout)

fout.close()

import shutil
dir_path = os.environ["DATA_DIR"]
dir_path += f"/{args.Era}/Sample/ForSNU/"
print(dir_path)
shutil.move(f"./{args.Sample}.txt", dir_path+f"{args.Sample}.txt")
