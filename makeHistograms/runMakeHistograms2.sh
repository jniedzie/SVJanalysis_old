binning=binning.json
samples_description=samples_skim1.json
processor=Histogram2
output_directory=./


## For tchannel
variables=(
    1000
    3000
    4000
    6000
)

for var in ${variables[@]}; do
  sample=tchannel_mMed-${var}_mDark-20_rinv-0.3_alpha-peak_extended
  str=`rdump -f /work/fleble/t_channel_samples/102X/mMed-${var}_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim1/merged.root -b Efficiencies/totalEfficiency`
  efficiency=`echo ${str} | cut -d[ -f2 | cut -d] -f1`
  python makeHistograms.py --binning ${binning} --samplesDescription ${samples_description} --samples ${sample} --processor ${processor} --outputDirectory ${output_directory} -e ${efficiency} -t PFnano102X
done


## For QCD
variables=(
    170_300
    300_470
    470_600
    600_800
    800_1000
    1000_1400
    1400_1800
    1800_2400
    2400_3200
    3200_Inf
)

for var in ${variables[@]}; do
  sample=QCD_${var}_extended
  var=`python -c 'print("'${var}'".replace("_", "to"))'`
  str=`rdump -f /work/fleble/QCD_samples/106X/QCD_Pt_${var}_TuneCP5_13TeV_pythia8/NANOAODJMAR_skim1/merged_lite1.root -b Cuts/Efficiency`
  efficiency=`echo ${str} | cut -d[ -f2 | cut -d] -f1`
  python makeHistograms.py --binning ${binning} --samplesDescription ${samples_description} --samples ${sample} --processor ${processor} --outputDirectory ${output_directory} -e ${efficiency} -t PFnano106X
done
