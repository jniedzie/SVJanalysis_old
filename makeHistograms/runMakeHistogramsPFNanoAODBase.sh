#!/bin/bash


make_file_arg() {

    str=""
    for file in ${@:1:($#-1)}; do
        str=${str}${file},
    done
    str=${str}${@:$#}
    echo ${str}
}


binning=binning.json
processor=HistogramPFNanoAODBase

model=t-channel_mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8

input_file=/pnfs/psi.ch/cms/trivcat//store/user/fleble/SVJ/t_channel/samples/2018/PFNANOAODSKIM/step6_PFNANOAODSKIM_${model}_n-1000_part-1to100.root
# Example using file name expansion
# Not the best t-channel because files are small, better using a merged file because individual files are so small
#input_file=$(make_file_arg /pnfs/psi.ch/cms/trivcat//store/user/fleble/SVJ/t_channel/samples/2018/PFNANOAODSKIM/step6_PFNANOAODSKIM_${model}_n-1000_part-[0-9]*.root)

output_file=/pnfs/psi.ch/cms/trivcat//store/user/fleble/SVJ/t_channel/histograms/2018/PFNANOAODSKIM/HISTOGRAMS_${model}.root

python makeHistograms.py -i ${input_file} -o ${output_file} --binning ${binning} --processor ${processor}
