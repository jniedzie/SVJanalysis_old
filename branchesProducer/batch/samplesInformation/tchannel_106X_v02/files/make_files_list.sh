#!/bin/bash

path=root://t3se01.psi.ch:1094//store/user/fleble/SVJ/t_channel/samples/2018/PFNANOAODSKIM/

models=(
    t-channel_mMed-500_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
    t-channel_mMed-1000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
    t-channel_mMed-1500_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
    t-channel_mMed-2000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
    t-channel_mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8
)

for model in ${models[@]}; do
    genRange 100 ${path}/PFNANOAODSKIM_${model}_part- .root > ${model}_files.txt
done
