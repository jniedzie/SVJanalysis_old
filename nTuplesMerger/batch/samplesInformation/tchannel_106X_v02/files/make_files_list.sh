#!/bin/bash

source ../../make_files_list.sh

pf_nano_skim_file_list_directory=/t3home/fleble/SVJ/SVJanalysis/branchesProducer/batch/samplesInformation/tchannel_106X_v02/files/

models=(
    t-channel_mMed-500_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
    t-channel_mMed-1000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
    t-channel_mMed-1500_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
    t-channel_mMed-2000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
    t-channel_mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
)

for model in ${models[@]}; do
    file_name=${model}_files.txt
    pf_nano_skim_file_full_path=${pf_nano_skim_file_list_directory}/${file_name}

    $(make_files_list ${pf_nano_skim_file_full_path} ${file_name})
done
