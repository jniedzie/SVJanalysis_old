#!/bin/bash

outputPath=/work/fleble/t_channel_samples/102X/

#python makePreSelection.py -i /pnfs/psi.ch/cms/trivcat/store/user/fleble/SVJ/t_channel/samples/102X/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR/merged.root -o ${outputPath}/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim1/merged.root -t PFnano102X -c 10000
python makePreSelection.py -i /pnfs/psi.ch/cms/trivcat/store/user/fleble/SVJ/t_channel/samples/102X/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR/1.root -o ${outputPath}/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim1/tmp.root -t PFnano102X -c 10000
