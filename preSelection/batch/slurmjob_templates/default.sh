#!/bin/bash
# 
#SBATCH --partition {partition}
#SBATCH --account=t3
#SBATCH --job-name={jobName}
#SBATCH --mem={memory}
#SBATCH --time {timeLimit}
#SBATCH -o {logDir}/%x-%j.out    # replace default slurm-SLURM_JOB_ID.out; %x is a job-name (or script name when there is no job-name)
#SBATCH -e {logDir}/%x-%j.err    # replace default slurm-SLURM_JOB_ID.err 


get_last_characters() {

    local str=$1
    local n=$2
    local ichar=(${#str}-${n})
    local last_chars=${str:$ichar:$n}
    echo ${last_chars}
}


copy_file() {

    local input_file=$1
    local destination=$2  # filename or directory

    echo "Copying ${input_file} into ${destination} ..."
    xrdcp -f ${input_file} ${destination}
}


copy_files_in_directory() {

    local file_with_files_to_copy=$1
    local destination=$2  # a directory
    files_to_copy=$(cat ${file_with_files_to_copy})
    for file in ${files_to_copy[@]}; do
        copy_file ${file} ${destination}
    done
}


list_files_in_directory() {

   local directory=$1
   local output_file=$2

   echo `ls ${directory}/` | tr " " "\n" > ${directory}/tmp.txt
   awk -v var="${directory}/" '//f{$0 = var $0}{print}' ${directory}/tmp.txt > ${output_file}
}


echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME
echo ""


# Define temporary working directory
tmpdir=/scratch/$USER/${SLURM_JOB_ID}
#tmpdir=/scratch/$USER/tmp/
echo -e "\nMaking temporary working directory ${tmpdir}"
mkdir -p ${tmpdir}


# Copy input files into working directory...
if {xrdcpInputFile}; then
    input_file0={inputFile}
    if [ "$(get_last_characters ${input_file0} 5)" == ".root" ]; then
        input_file=${tmpdir}/input.root
        copy_file ${input_file0} ${input_file}
    elif [ "$(get_last_characters ${input_file0} 4)" == ".txt" ]; then
        input_file=${tmpdir}/input_file.txt
        copy_files_in_directory ${input_file0} ${tmpdir}/
        list_files_in_directory ${tmpdir} ${input_file}
    else
        echo "ERROR: ${input_file0} must have extension .root or .txt!"
        exit 1
    fi

# ... or read directly input files from where they are
else
    input_file={inputFile}
fi


# Make pre-selections
echo -e "\nMaking pre-selections..."
cd ../../
output_file=${tmpdir}/output.root
python makePreSelection.py -i ${input_file} -xsec {xsec} -o ${output_file} -t PFNanoAOD_106X_v02 -c 300 -m 20 -p Preselection_tchannel


# Copy produced file into storage element
echo -e "\nCopying output file into storage element..."
copy_file ${output_file} {outputFile}


# Cleaning of temporary working directory
echo -e "\nDeleting temporary working directory ${tmpdir}"
rm -rf ${tmpdir}
