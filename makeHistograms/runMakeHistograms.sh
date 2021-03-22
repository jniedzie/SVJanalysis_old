binning=binning.json
samplesDescription=samples.json
outputDir=./

processor=Histogram1

samples=(tchannel_mMed-3000_mDark-20_rinv-0p3)


for sample in ${samples[@]}; do
  python makeHistograms.py --mode recreate --binning ${binning} --samplesDescription ${samplesDescription} --samples ${sample} --processor ${processor} --outputDirectory ${outputDir}
done
