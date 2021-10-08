#!/bin/bash

input_file=/pnfs/psi.ch/cms/trivcat/store/user/fleble/SVJ/t_channel/samples/2018/PFNANOAOD/step5_PFNANOAOD_t-channel_mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000_part-1.root

python makePreSelection.py -i ${input_file} -o test.root -t PFNanoAOD_106X_v02 -c 100 -m 4 -p Preselection_tchannel -xsec 0.01891
