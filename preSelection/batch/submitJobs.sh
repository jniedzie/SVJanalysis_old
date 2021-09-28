#!/bin/bash

source ../../utilities/bashUtilities.sh
source ../../utilities/slurmJobSubmission/slurmJobSubmissionUtilities.sh


# Sample to run
sample=tchannel


get_output_directory_for_sample() {

    local sample=$1
    
    if [ "${sample}" == "tchannel" ]; then
        output_directory=root://t3se01.psi.ch:1094//store/user/fleble/SVJ/t_channel/samples/2018
    elif [ "${sample}" == "QCD" ]; then
        output_directory=root://t3se01.psi.ch:1094//store/user/fleble/SVJ/QCD/samples/UL2018
    elif [ "${sample}" == "ttbar" ]; then
        output_directory=root://t3se01.psi.ch:1094//store/user/fleble/SVJ/ttbar/samples/UL2018
    else
        exit 1
    fi

    echo ${output_directory}
}


# Make slurm job files directory
slurmjob_dir=slurmjob_files
make_dir ${slurmjob_dir}

# Make log files directory
log_dir=${PWD}/logs
make_dir ${log_dir}

datatier=PFNANOAODSKIM
file_prefix=${datatier}_

output_directory=$(get_output_directory_for_sample ${sample})


# Change these numbers to only run for some files
first_part_number=-1
last_part_number=-1

xsec_file=cross_sections.txt


###################
###  t-channel  ###
###################

if [ "${sample}" == "tchannel" ]; then

    slurm_job_template_name=default
    partition=short
    time_limit=01:00:00
    memory=3000M

    xrdcp_files=false
    force_recreate=false
  
    # List of models to process
    mMedList=(  500 1000 1500 2000 3000)
    mDarkList=( 20  20   20   20   20  )
    rinvList=(  0.3 0.3  0.3  0.3  0.3 )
    alphaList=(peak peak peak peak peak)
    
    # Loop over all models
    for ((imodel=0; imodel<5; imodel++)); do
        echo ""
  
        # Get model parameters
        mMed=${mMedList[${imodel}]}
        mDark=${mDarkList[${imodel}]}
        rinv=${rinvList[${imodel}]}
        alpha=${alphaList[${imodel}]}

        # Define model name
        model=t-channel_mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1_13TeV-madgraphMLM-pythia8
    
        relative_path_to_files_list=inputFilesLists/tchannel_106X_v02/${model}_files.txt
        files_list_file=${PWD}/${relative_path_to_files_list}
        xsec=$(cat ${xsec_file} | grep ${relative_path_to_files_list} | cut -d' ' -f 2-)

        specific_flags=-xsec:${xsec}

        # Submit jobs
        submit_one_job_per_file ${slurm_job_template_name} ${partition} ${time_limit} ${memory} ${model} ${files_list_file} ${output_directory} ${datatier} ${file_prefix} ${xrdcp_files} ${force_recreate} ${first_part_number} ${last_part_number} ${specific_flags}
    done
fi


#############
###  QCD  ###
#############

if [ "${sample}" == "QCD" ]; then

    # Overiding default last part number to loop over 5% of the dataset
    last_part_number=.05

    slurm_job_template_name=default
    partition=standard
    time_limit=12:00:00
    memory=4000M

    xrdcp_files=true
    force_recreate=false

    # List of models to process
    pt1List=(170 300 470 600 800  1000 1400 1800 2400 3200)
    pt2List=(300 470 600 800 1000 1400 1800 2400 3200 Inf )
    
    # Loop over all models
    for ((imodel=0; imodel<10; imodel++)); do
        echo ""

        # Get model parameters
        pt1=${pt1List[${imodel}]}
        pt2=${pt2List[${imodel}]}

        # Define model name
        model=QCD_pt_${pt1}to${pt2}

        relative_path_to_files_list=inputFilesLists/QCD_106X_v02/${model}_files.txt
        files_list_file=${PWD}/${relative_path_to_files_list}
        xsec=$(cat ${xsec_file} | grep ${relative_path_to_files_list} | cut -d' ' -f 2-)

        specific_flags=-xsec:${xsec}

        # Submit jobs
        submit_one_job_per_file ${slurm_job_template_name} ${partition} ${time_limit} ${memory} ${model} ${files_list_file} ${output_directory} ${datatier} ${file_prefix} ${xrdcp_files} ${force_recreate} ${first_part_number} ${last_part_number} ${specific_flags}
    done
fi


###############
###  ttbar  ###
###############

if [ "${sample}" == "ttbar" ]; then

    # Overiding default last part number to loop over 1% of the dataset
    last_part_number=.01

    slurm_job_template_name=default
    partition=standard
    time_limit=12:00:00
    memory=4000M

    xrdcp_files=true
    force_recreate=false

    # Define model name
    model=ttbar

    relative_path_to_files_list=inputFilesLists/ttbar_106X_v02/${model}_files.txt
    files_list_file=${PWD}/${relative_path_to_files_list}
    xsec=$(cat ${xsec_file} | grep ${relative_path_to_files_list} | cut -d' ' -f 2-)

    specific_flags=-xsec:${xsec}

    # Submit jobs
    submit_one_job_per_file ${slurm_job_template_name} ${partition} ${time_limit} ${memory} ${model} ${files_list_file} ${output_directory} ${datatier} ${file_prefix} ${xrdcp_files} ${force_recreate} ${first_part_number} ${last_part_number} ${specific_flags}
fi

