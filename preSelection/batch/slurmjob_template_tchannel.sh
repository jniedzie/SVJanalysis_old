#!/bin/bash
# 
#SBATCH -p wn
#SBATCH --account=t3
#SBATCH --job-name=tchannel_mMed-{mMed}_mDark-{mDark}_rinv-{rinv}_alpha-{alpha}_yukawa-1
#SBATCH --mem=3000M
#SBATCH --time 1:00:00
#SBATCH -o /work/fleble/t_channel_samples/102X/mMed-{mMed}_mDark-{mDark}_rinv-{rinv}_alpha-{alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim1/%x-%j.out    # replace default slurm-SLURM_JOB_ID.out; %x is a job-name (or script name when there is no job-name)
#SBATCH -e /work/fleble/t_channel_samples/102X/mMed-{mMed}_mDark-{mDark}_rinv-{rinv}_alpha-{alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim1/%x-%j.err    # replace default slurm-SLURM_JOB_ID.err 


echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME
echo ""

## Each worker node has local /scratch space to be used during job run
#mkdir -p /scratch/$USER/${SLURM_JOB_ID}
#export TMPDIR=/scratch/$USER/${SLURM_JOB_ID}


## Make computation
out_dir=/work/fleble/t_channel_samples/102X/mMed-{mMed}_mDark-{mDark}_rinv-{rinv}_alpha-{alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR_skim1

input_file=/pnfs/psi.ch/cms/trivcat/store/user/fleble/SVJ/t_channel/samples/102X/mMed-{mMed}_mDark-{mDark}_rinv-{rinv}_alpha-{alpha}_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR/merged.root

python ../../makePreSelection.py -i ${input_file} -o ${out_dir}/merged.root -t PFnano102X -c 10000 -p Preselection_tchannel


## Cleaning of temporal working dir when job was completed:
#rmdir  -rf /scratch/$USER/${SLURM_JOB_ID}
