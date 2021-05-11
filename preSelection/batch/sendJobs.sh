#!/bin/bash

# Sample to run
sample=QCD_schannelCuts


make_dir() {
  ### Make directory given as argument if it does not exist ###

  local directory=$1
  if [ ! -z ${directory} ]; then mkdir -p ${directory}; fi
}

make_slurmjob_file() {
  ###
  # Receive
  #   * a sample name
  #   * a set of arguments
  #   * an output filename
  # and make the corresponding slurmjob file
  ###

  echo "Creating slurm job file ${@: -1}..."
  python make_slurmjob_file.py $@
  if [ ! "$?" == "0"  ]; then exit 1; fi
}


submit_job() {
  ###
  # Submit job
  # Arguments: * slurmjob files directory
  #            * slurmjob file name in that directory
  ###

  local slurmjob_dir=$1
  local slurmjob_file=$2

  echo "Submitting job from file ${slurmjob_dir}/${slurmjob_file}..."
  cd ${slurmjob_dir}
  sbatch ${slurmjob_file}
  cd -

}


## Make slurmjob files directory
slurmjob_dir=slurmjob_files
make_dir ${slurmjob_dir}


#############
###  QCD  ###
#############

if [ "${sample}" == "QCD" ]; then

  pt1List=(3200)
  pt2List=(Inf)
  #pt1List=(170 300 470 600 800 1000 1400 1800 2400 3200)
  #pt2List=(300 470 600 800 1000 1400 1800 2400 3200 Inf)
  dataset=_lite1
  
  for ((i=0; i<1; i++)); do
    echo ""
    pt1=${pt1List[${i}]}
    pt2=${pt2List[${i}]}
  
    # Make output directory if it does not exist
    outDir=/work/fleble/QCD_samples/106X/QCD_Pt_${pt1}to${pt2}_TuneCP5_13TeV_pythia8/NANOAODJMAR_skim1
    make_dir ${outDir}
  
    # Make the slurm job file
    slurmjob_file=${sample}_${pt1}_${pt2}${dataset}.sh
    make_slurmjob_file ${sample} ${pt1} ${pt2} ${dataset} ${slurmjob_dir}/${slurmjob_file}
  
    # Submit job
    submit_job ${slurmjob_dir} ${slurmjob_file}
  done

fi


#######################
###  QCD  schannel  ###
#######################

if [ "${sample}" == "QCD_schannelCuts" ]; then

  pt1List=(300 470 600 800)
  pt2List=(470 600 800 1000)
  f1List=(10 20 30 40)
  f2List=(10 10 10 10)
  
  
  for ((i=0; i<4; i++)); do
    echo ""
    pt1=${pt1List[${i}]}
    pt2=${pt2List[${i}]}
    for ((j=0; j<4; j++)); do
      echo ""
      f1=${f1List[${j}]}
      f2=${f2List[${j}]}
  
      # Make output directory if it does not exist
      outDir=/work/fleble/QCD_samples/106X/QCD_Pt_${pt1}to${pt2}_TuneCP5_13TeV_pythia8/NANOAODJMAR_skimschannel
      make_dir ${outDir}
  
      # Make the slurm job file
      slurmjob_file=${sample}_${pt1}_${pt2}_${f1}_${f2}.sh
      make_slurmjob_file ${sample} ${pt1} ${pt2} ${f1} ${f2} ${slurmjob_dir}/${slurmjob_file}
  
      # Submit job
      submit_job ${slurmjob_dir} ${slurmjob_file}
    done
  done

fi


###################
###  t-channel  ###
###################

if [ "${sample}" == "tchannel" ]; then
  mMedList=(3000)
  mDarkList=(20)
  rinvList=(0.3)
  alphaList=(peak)
  
  
  for ((i=0; i<1; i++)); do
    echo ""
    mMed=${mMedList[${i}]}
    mDark=${mDarkList[${i}]}
    rinv=${rinvList[${i}]}
    alpha=${alphaList[${i}]}
  
    # Make output directory if it does not exist
    outDir=/work/fleble/t_channel_samples/102X/mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim1
    make_dir ${outDir}
  
    # Make the slurm job file
    slurmjob_file=${sample}_mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1.sh
    make_slurmjob_file ${sample} ${mMed} ${mDark} ${rinv} ${alpha} ${slurmjob_dir}/${slurmjob_file}
  
    # Submit job
    submit_job ${slurmjob_dir} ${slurmjob_file}
  done

fi


###################
###  s-channel  ###
###################

if [ "${sample}" == "schannel" ]; then

  mMedList=(3500)
  mDarkList=(40)
  rinvList=(0.3)
  alphaList=(peak)
  
  
  for ((i=0; i<1; i++)); do
    echo ""
    mMed=${mMedList[${i}]}
    mDark=${mDarkList[${i}]}
    rinv=${rinvList[${i}]}
    alpha=${alphaList[${i}]}
  
    # Make output directory if it does not exist
    outDir=/work/fleble/s_channel_samples/102X/mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim
    make_dir ${outDir}
  
    # Make the slurm job file
    slurmjob_file=${sample}_mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1.sh
    make_slurmjob_file ${sample} ${mMed} ${mDark} ${rinv} ${alpha} ${slurmjob_dir}/${slurmjob_file}
  
    # Submit job
    submit_job ${slurmjob_dir} ${slurmjob_file}
  done

fi
