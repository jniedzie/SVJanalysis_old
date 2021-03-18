variables=variables.json
dataframes=filtersDefines.json
binning=binning.json
samplesDescription=samples.json
outputDir=/eos/user/f/fleble/SVJ/data/histograms/

samples=(QCD_Pt_170to300)

for sample in ${samples[@]}; do
  echo $sample
  python makeHistograms.py --mode RECREATE --variables ${variables} --dataframes ${dataframes} --binning ${binning} --samplesDescription ${samplesDescription} --samples ${sample} --outputDirectory ${outputDir}
done
