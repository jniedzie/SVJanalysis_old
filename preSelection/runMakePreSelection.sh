#!/bin/bash

sepath=/pnfs/psi.ch/cms/trivcat/store/user/fleble

chunksize=10000

mMed=3000
mDark=20
rinv=0.3
alpha=peak

var=mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1_13TeV-madgraphMLM-pythia8
ifile=${sepath}/SVJ/t_channel/samples/102X/${var}/NANOAODJMAR/merged.root
ofile=./test.root

python makePreSelection.py -i ${ifile} -o ${ofile} -t PFnano102X -c ${chunksize} -p Preselection_tchannel
