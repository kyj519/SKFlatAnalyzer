#!/usr/bin/env python3

import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='')
parser.add_argument('-e', dest='Era', default='2018')
args = parser.parse_args()

mc_list = {#'DYJets_MG':'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',
           #'tw_antitop':'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8',
           #'tW_top':'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8',
           #'TTbb':'TTbb_4f_TTToSemiLeptonic_TuneCP5-Powheg-Openloops-Pythia8',
           #'TTLL_powheg':'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
           #'TTLJ_powheg':'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
           'WJet_MG':'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8',
}

ssh_key = "ssh isyoon@kcms-t2.knu.ac.kr"
crab_path = "/pnfs/knu.ac.kr/data/cms/store/user/iyoon/SKFlat/"

for mc in mc_list:
    target = f"{crab_path}/{args.Era}/{mc_list[mc]}/SKFlat_Run2UltraLegacy_v2"

    #find lastest crab output
    result = subprocess.Popen(f"{ssh_key} ls {target}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    result_list = result[0].decode().split('\n')

    print(result_list[-2])

    source = f"{target}/{result_list[-2]}"
    target = f"/gv0/Users/isyoon/Run2UltraLegacy_v2/{args.Era}/{mc_list[mc]}/SKFlat_Run2UltraLegacy_v2"
    print(source)
    print(target)
    
    #make new screen
    #clean session first
    os.system(f"screen -X -S {mc} kill")
    
    proc = subprocess.Popen(f"screen -dmS {mc} rsync -rv --progress isyoon@kcms-t2.knu.ac.kr:{source} {target}", shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    
    
