
## Make the slurm job file
#python make_slurmjob_file.py 170 300 _lite10 slurmjob_file.sh
#
## Submit job
#sbatch slurmjob_file.sh

#############
###  QCD  ###
#############

#pt1List=(3200)
#pt2List=(Inf)
##pt1List=(170 300 470 600 800 1000 1400 1800 2400 3200)
##pt2List=(300 470 600 800 1000 1400 1800 2400 3200 Inf)
#dataset=_lite1
#
#for ((i=0; i<1; i++)); do
#    echo ""
#    pt1=${pt1List[${i}]}
#    pt2=${pt2List[${i}]}
#
#    # Make output directory if it does not exist
#    out_dir=/work/fleble/QCD_samples/106X/QCD_Pt_${pt1}to${pt2}_TuneCP5_13TeV_pythia8/NANOAODJMAR_skim1
#    if [ ! -z ${out_dir} ]; then mkdir -p ${out_dir}; fi
#
#    # Make the slurm job file
#    slurmjob_dir=slurmjob_files
#    slurmjob_file=QCD_${pt1}_${pt2}${dataset}.sh
#    echo "Creating slurm job file ${slurmjob_file}..."
#    python make_slurmjob_file.py QCD ${pt1} ${pt2} ${dataset} ${slurmjob_dir}/${slurmjob_file}
#    if [ ! "$?" == "0"  ]; then exit 1; fi
#
#    # Submit job
#    echo "Submitting job from file ${slurmjob_file}..."
#    cd ${slurmjob_dir}
#    sbatch ${slurmjob_file}
#    cd ..
#done


#######################
###  QCD  schannel  ###
#######################

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
    out_dir=/work/fleble/QCD_samples/106X/QCD_Pt_${pt1}to${pt2}_TuneCP5_13TeV_pythia8/NANOAODJMAR_skimschannel
    if [ ! -z ${out_dir} ]; then mkdir -p ${out_dir}; fi

    # Make the slurm job file
    slurmjob_dir=slurmjob_files
    sample=QCD_schannelCuts
    slurmjob_file=${sample}_${pt1}_${pt2}_${f1}_${f2}.sh
    echo "Creating slurm job file ${slurmjob_dir}/${slurmjob_file}..."
    python make_slurmjob_file.py ${sample} ${pt1} ${pt2} ${f1} ${f2} ${slurmjob_dir}/${slurmjob_file}
    if [ ! "$?" == "0"  ]; then exit 1; fi

    # Submit job
    echo "Submitting job from file ${slurmjob_dir}/${slurmjob_file}..."
    cd ${slurmjob_dir}
    sbatch ${slurmjob_file}
    cd ..
  done
done



###################
###  t-channel  ###
###################

#mMedList=(3000)
#mDarkList=(20)
#rinvList=(0.3)
#alphaList=(peak)
#
#
#for ((i=0; i<1; i++)); do
#    echo ""
#    mMed=${mMedList[${i}]}
#    mDark=${mDarkList[${i}]}
#    rinv=${rinvList[${i}]}
#    alpha=${alphaList[${i}]}
#
#    # Make output directory if it does not exist
#    out_dir=/work/fleble/t_channel_samples/102X/mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim1
#    if [ ! -z ${out_dir} ]; then mkdir -p ${out_dir}; fi
#
#    # Make the slurm job file
#    slurmjob_dir=slurmjob_files
#    slurmjob_file=tchannel_mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1.sh
#    echo "Creating slurm job file ${slurmjob_file}..."
#    python make_slurmjob_file.py tchannel ${mMed} ${mDark} ${rinv} ${alpha} ${slurmjob_dir}/${slurmjob_file}
#
#    # Submit job
#    echo "Submitting job from file ${slurmjob_file}..."
#    cd ${slurmjob_dir}
#    sbatch ${slurmjob_file}
#    cd ..
#done

###################
###  s-channel  ###
###################

#mMedList=(3500)
#mDarkList=(40)
#rinvList=(0.3)
#alphaList=(peak)
#
#
#for ((i=0; i<1; i++)); do
#    echo ""
#    mMed=${mMedList[${i}]}
#    mDark=${mDarkList[${i}]}
#    rinv=${rinvList[${i}]}
#    alpha=${alphaList[${i}]}
#
#    # Make output directory if it does not exist
#    out_dir=/work/fleble/s_channel_samples/102X/mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim
#    if [ ! -z ${out_dir} ]; then mkdir -p ${out_dir}; fi
#
#    # Make the slurm job file
#    slurmjob_dir=slurmjob_files
#    slurmjob_file=schannel_mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1.sh
#    echo "Creating slurm job file ${slurmjob_file}..."
#    python make_slurmjob_file.py schannel ${mMed} ${mDark} ${rinv} ${alpha} ${slurmjob_dir}/${slurmjob_file}
#
#    # Submit job
#    echo "Submitting job from file ${slurmjob_file}..."
#    cd ${slurmjob_dir}
#    sbatch ${slurmjob_file}
#    cd ..
#done
