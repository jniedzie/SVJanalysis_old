file=/work/fleble/t_channel_samples/102X/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim1/merged.root

python produceBranches.py --inputfile ${file} --chunksize 100000 -o test.root
