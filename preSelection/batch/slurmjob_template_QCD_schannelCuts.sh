#!/bin/bash
# 
#SBATCH -p wn
#SBATCH --account=t3
#SBATCH --job-name=QCD_{pt1}_{pt2}
#SBATCH --mem=4000M
#SBATCH --time 04:00:00
#SBATCH -o /work/fleble/QCD_samples/106X/QCD_Pt_{pt1}to{pt2}_TuneCP5_13TeV_pythia8/NANOAODJMAR_skimschannel/%x-%j.out    # replace default slurm-SLURM_JOB_ID.out; %x is a job-name (or script name when there is no job-name)
#SBATCH -e /work/fleble/QCD_samples/106X/QCD_Pt_{pt1}to{pt2}_TuneCP5_13TeV_pythia8/NANOAODJMAR_skimschannel/%x-%j.err    # replace default slurm-SLURM_JOB_ID.err 


echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME
echo ""

## Each worker node has local /scratch space to be used during job run
tmpdir=/scratch/fleble/QCD_${SLURM_JOB_ID}/
if [ -d ${tmpdir} ]; then rm -rf ${tmpdir}; fi
mkdir -p ${tmpdir}

cd ${tmpdir}

## Copy input files into the tmp directory
files=`cat /t3home/fleble/SVJ/SVJanalysis/preSelection/batch/samplesLocation/QCD_{pt1}_{pt2}_tmp.txt | head -n {f1} | tail -n {f2}`
for file in ${files[@]}; do
    xrdcp ${file} ${tmpdir}
done

## Make list of samples pointing to the files that were copied
echo `ls ${tmpdir}` | tr " " "\n" > samples0.txt
awk -v var="${tmpdir}" '//f{$0 = var $0}{print}' samples0.txt > samples.txt

## Make computation
outDir=/work/fleble/QCD_samples/106X/QCD_Pt_{pt1}to{pt2}_TuneCP5_13TeV_pythia8/NANOAODJMAR_skimschannel

python /t3home/fleble/SVJ/SVJanalysis/preSelection/makePreSelection.py -i samples.txt -o ${outDir}/merged_{f1}_{f2}.root -t PFnano106X -c 200000 -p Preselection_schannel -m 15


## Clean temporary working dir when job is completed
rm -rf ${tmpdir}
