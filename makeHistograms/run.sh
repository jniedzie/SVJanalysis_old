variables=variables.json
dataframes=filtersDefines.json
binning=binning.json
samplesDescription=samples.json
outputDir=/eos/user/f/fleble/SVJ/data/histograms_test/

python makeHistograms.py --mode RECREATE --variables ${variables} --dataframes ${dataframes} --binning ${binning} --samplesDescription ${samplesDescription} --samples tchannel_mMed-3000_mDark-20_rinv-0p3 --outputDirectory ${outputDir}
#--lumi 59725
