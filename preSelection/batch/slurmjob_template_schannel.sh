#!/bin/bash
# 
#SBATCH -p wn
#SBATCH --account=t3
#SBATCH --job-name=schannel_mMed-{mMed}_mDark-{mDark}_rinv-{rinv}_alpha-{alpha}
#SBATCH --mem=3000M
#SBATCH --time 02:00:00
#SBATCH -o /work/fleble/s_channel_samples/102X/mMed-{mMed}_mDark-{mDark}_rinv-{rinv}_alpha-{alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim/%x-%j.out    # replace default slurm-SLURM_JOB_ID.out; %x is a job-name (or script name when there is no job-name)
#SBATCH -e /work/fleble/s_channel_samples/102X/mMed-{mMed}_mDark-{mDark}_rinv-{rinv}_alpha-{alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim/%x-%j.err    # replace default slurm-SLURM_JOB_ID.err 


echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME
echo ""

## Each worker node has local /scratch space to be used during job run
tmpdir=/scratch/fleble/schannel_61-100/
if [ -d ${tmpdir} ]; then rm -rf ${tmpdir}; fi
mkdir -p ${tmpdir}

cd ${tmpdir}

files=`cat /t3home/fleble/SVJ/SVJanalysis/preSelection/batch/samplesLocation/schannel_61-100.txt`
for file in ${files[@]}; do
    xrdcp ${file} ${tmpdir}
done

echo `ls ${tmpdir}` | tr " " "\n" > samples0.txt
awk -v var="${tmpdir}" '//f{$0 = var $0}{print}' samples0.txt > samples.txt

## Make computation
out_dir=/work/fleble/s_channel_samples/102X/mMed-{mMed}_mDark-{mDark}_rinv-{rinv}_alpha-{alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim

python /t3home/fleble/SVJ/SVJanalysis/preSelection/makePreSelection.py -i samples.txt -o ${out_dir}/61-100.root -t PFnano102X -c 200000 -p Preselection_schannel


## Cleaning of temporal working dir when job was completed:
rm -rf ${tmpdir}
