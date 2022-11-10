#!/usr/bin/env python3

import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='')
parser.add_argument('-e', dest='Era', default='2018')
args = parser.parse_args()

mc_list = {#'TTLJ_powheg':'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
           #'TTLL_powheg':'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
           #'tw_antitop':'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8',
           #'tW_top':'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8',
           #'DYJets_MG':'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',
           #'WJet_MG':'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8',
           #'QCD_bEn_HT100to200':'QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8',
           #'QCD_bEn_HT200to300':'QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8',
           #'QCD_bEn_HT300to500':'QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8',
           #'QCD_bEn_HT500to700':'QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8',
           #'QCD_bEn_HT700to1000':'QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8',
           #'QCD_bEn_HT1000to1500':'QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8',
           #'QCD_bEn_HT1500to2000':'QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8',
           #'QCD_bEn_HT2000toInf':'QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8',
           #'TTLJ_WtoCB_powheg0':'ttj_Vcb_NLO_FXFX_Summer20UL18_LHEGEN',
           #'TTLJ_WtoCB_powheg1':'ttj_Vcb_NLO_FXFX_Summer20UL18_GEN_extra',
           #'TTLJ_powheg_hdampDown':'TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8',
           #'TTLJ_powheg_hdampUp':'TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8',
           #'TTLJ_powheg_m171p5':'TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8',
           #'TTLJ_powheg_m173p5':'TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8',
           'TTLJ_powheg_CP5Down':'TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8',
           'TTLJ_powheg_CP5Up':'TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8',
           #'TTbb':'TTbb_4f_TTToSemiLeptonic_TuneCP5-Powheg-Openloops-Pythia8',
}

ssh_key = "ssh isyoon@cms01.knu.ac.kr"
crab_path = "/pnfs/knu.ac.kr/data/cms/store/user/iyoon/SKFlat/"

for mc in mc_list:
    print(mc)
    
    target = f"{crab_path}/{args.Era}/{mc_list[mc]}/SKFlat_Run2UltraLegacy_v3"

    #find lastest crab output
    result = subprocess.Popen(f"{ssh_key} ls {target}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    result_list = result[0].decode().split('\n')

    print(result_list[-2])

    source = f"{target}/{result_list[-2]}"
    target = f"/gv0/Users/isyoon/Run2UltraLegacy_v3/{args.Era}/{mc_list[mc]}/SKFlat_Run2UltraLegacy_v3"
    print(source)
    print(target)

    os.makedirs(target, exist_ok=True)
    
    #make new screen
    #clean session first
    os.system(f"screen -X -S {mc} kill")
    
    print(f"screen -dmS {mc} rsync -rv --progress isyoon@cms01.knu.ac.kr:{source} {target}")
    
    proc = subprocess.Popen(f"screen -dmS {mc} rsync -rv --progress isyoon@kcms-t2.knu.ac.kr:{source} {target}", shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    
    
