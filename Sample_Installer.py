#!/usr/bin/env python3

import argparse
import os
import subprocess
import ROOT

parser = argparse.ArgumentParser(description='Sample_Installer')
parser.add_argument('-e', dest='era', default="2018")
parser.add_argument('-s', dest='sample', default="")
parser.add_argument('-p', dest='Chk_Private', type=bool, default=False)
args = parser.parse_args()

if args.era=="2016a": args.era="2016preVFP"
if args.era=="2016b": args.era="2016postVFP"

sample_dict = {'TTLJ_WtoCB_powheg':'ttj_Vcb_NLO_FXFX_Summer20UL18',
               'TTLJ_powheg_hdampDown':'TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8',
               'TTLJ_powheg_hdampUp':'TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8',
               'TTLJ_powheg_mtop171p5':'TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8',
               'TTLJ_powheg_mtop173p5':'TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8',
               'TTLJ_powheg_CP5Down':'TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8',
               'TTLJ_powheg_CP5Up':'TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8',
               'TTLL_powheg_hdampDown':'TTTo2L2Nu_hdampDOWN_TuneCP5_13TeV-powheg-pythia8',
               'TTLL_powheg_hdampUp':'TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8',
               'TTLL_powheg_mtop171p5':'TTTo2L2Nu_mtop171p5_TuneCP5_13TeV-powheg-pythia8',
               'TTLL_powheg_mtop173p5':'TTTo2L2Nu_mtop173p5_TuneCP5_13TeV-powheg-pythia8',
               'TTLL_powheg_CP5Down':'TTTo2L2Nu_TuneCP5down_13TeV-powheg-pythia8',
               'TTLL_powheg_CP5Up':'TTTo2L2Nu_TuneCP5up_13TeV-powheg-pythia8',
               'TTbb':'TTbb_4f_TTToSemiLeptonic_TuneCP5-Powheg-Openloops-Pythia8',
               'ttHTobb':'ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8',}

cross_dict = {'TTLJ_WtoCB_powheg':1.281e-01,
              'TTLJ_powheg_hdampDown':365.34,
              'TTLJ_powheg_hdampUp':365.34,
              'TTLJ_powheg_mtop171p5':365.34,
              'TTLJ_powheg_mtop173p5':365.34,
              'TTLJ_powheg_CP5Down':365.34,
              'TTLJ_powheg_CP5Up':365.34,
              'TTLL_powheg_hdampDown':88.29,
              'TTLL_powheg_hdampUp':88.29,
              'TTLL_powheg_mtop171p5':88.29,
              'TTLL_powheg_mtop173p5':88.29,
              'TTLL_powheg_CP5Down':88.29,
              'TTLL_powheg_CP5Up':88.29,
              'TTbb':1.7466,
              'ttHTobb':0.2951,}

if args.sample == 'TTLJ_WtoCB_powheg':
    if args.era == '2017':
        sample_dict['TTLJ_WtoCB_powheg'] = 'ttj_Vcb_NLO_FXFX_Summer20UL17'
    elif args.era == '2016preVFP':
        sample_dict['TTLJ_WtoCB_powheg'] = 'ttj_Vcb_NLO_FXFX_Summer20UL16APV'
    elif args.era == '2016postVFP':
        sample_dict['TTLJ_WtoCB_powheg'] = 'ttj_Vcb_NLO_FXFX_Summer20UL16'

path=os.environ['SKFlat_WD']
path_sk_out=os.environ['SKFlatOutputDir']

f_get_eff_lumi=f"{path_sk_out}/Run2UltraLegacy_v3/GetEffLumi/{args.era}/GetEffLumi_{args.sample}.root"
#print(f_get_eff_lumi)

if os.path.isfile(f_get_eff_lumi):
    fin = ROOT.TFile(f_get_eff_lumi)
    histo_sumW = fin.Get("sumW")
    histo_sumSign = fin.Get("sumSign")
    
    nmc = histo_sumW.GetEntries()
    sumsign = histo_sumSign.GetBinContent(1)
    sumw = histo_sumW.GetBinContent(1)
    
    #print(nmc, sumsign, sumw)
    
    f_common_path=f"{path}/data/Run2UltraLegacy_v3/{args.era}/Sample/CommonSampleInfo/{args.sample}.txt"

    f_common = open(f_common_path, 'w')
    f_common.write('# alias PD xsec nmc sumsign sumw\n')
    f_common.write(f"{args.sample}\t{sample_dict[args.sample]}\t{cross_dict[args.sample]}\t{nmc}\t{sumsign}\t{sumw}\n")

    f_common.close()

else:
    f_common_path=f"{path}/data/Run2UltraLegacy_v3/{args.era}/Sample/CommonSampleInfo/{args.sample}.txt"

    f_common = open(f_common_path, 'w')
    f_common.write('# alias PD xsec nmc sumsign sumw\n')
    f_common.write(f"{args.sample}\t{sample_dict[args.sample]}\t{cross_dict[args.sample]}\t1\t1\t1\n")

    f_common.close()
    
    #make file list
    subprocess.check_output([f"{path}/List_SKFlat.py", '-e', f"{args.era}", '-s', f"{args.sample}", '-p', f"{args.Chk_Private}"])

    #run get eff lumi
    operation = f"nohup SKFlat.py -a GetEffLumi -e {args.era} -i {args.sample} -n 20 &"

    #print(operation)
    os.system(operation)
