#!/bin/bash

binning=binning.json
samples_description=samples_skim1.json
processor=HistogramPFNanoAODBase


## For tchannel

output_directory=/work/fleble/t_channel_histograms/skim1/
samples=(
    tchannel_mMed-1000_mDark-20_rinv-0.3_alpha-peak
    tchannel_mMed-3000_mDark-20_rinv-0.3_alpha-peak
    tchannel_mMed-4000_mDark-20_rinv-0.3_alpha-peak
    tchannel_mMed-6000_mDark-20_rinv-0.3_alpha-peak
)


for sample in ${samples[@]}; do
  python makeHistograms.py --binning ${binning} --samplesDescription ${samples_description} --samples ${sample} --processor ${processor} --outputDirectory ${output_directory} -e Efficiencies/totalEfficiency
done


## For QCD

output_directory=/work/fleble/QCD_histograms/skim1/
samples=(
    QCD_170_300
    QCD_300_470
    QCD_470_600
    QCD_600_800
    QCD_800_1000
    QCD_1000_1400
    QCD_1400_1800
    QCD_1800_2400
    QCD_2400_3200
    QCD_3200_Inf
)


for sample in ${samples[@]}; do
  python makeHistograms.py --binning ${binning} --samplesDescription ${samples_description} --samples ${sample} --processor ${processor} --outputDirectory ${output_directory} -e Cuts/Efficiency
done
