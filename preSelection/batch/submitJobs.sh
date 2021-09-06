#!/bin/bash

# Sample to run
sample=ttbar


get_output_directory_for_sample() {

    sample=$1
    
    if [ "${sample}" == "tchannel" ]; then
        output_directory=root://t3se01.psi.ch:1094//store/user/fleble/SVJ/t_channel/samples/2018
    elif [ "${sample}" == "QCD" ]; then
        output_directory=root://t3se01.psi.ch:1094//store/user/fleble/SVJ/QCD/samples/2017
    elif [ "${sample}" == "ttbar" ]; then
        output_directory=root://t3se01.psi.ch:1094//store/user/fleble/SVJ/ttbar/samples/2017
    else
        exit 1
    fi

    echo ${output_directory}
}


make_dir() {
  ### Make directory given as argument if it does not exist ###

  local directory=$1
  if [ ! -z ${directory} ]; then mkdir -p ${directory}; fi
}


get_file_list() {

  local file_with_file_list=$1
  echo $(cat ${file_with_file_list})
}


get_number_of_files() {

  echo $#
}


get_file_number_from_file_list() {
  
  local file_number=$(($1+2))
  local file_list=$2
  echo ${!file_number}
}


get_redirector() {
  echo ${output_file} | cut -d/ -f-4
}


get_path() {
  echo ${output_file} | cut -d/ -f4-
}


file_exist() {
  # Return 0 if file exist, else 54

  redirector=$1
  file_path=$2
  
  xrdfs ${redirector} ls ${file_path} > /dev/null 2>&1
  echo $?
}


make_slurmjob_file() {
  ###
  # Receive
  #   * a sample name
  #   * a set of arguments
  #   * an output filename
  # and make the corresponding slurmjob file
  ###

  echo "Creating slurm job file ${@: -1} ..."
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

  echo "Submitting job from file ${slurmjob_dir}/${slurmjob_file} ..."
  cd ${slurmjob_dir}
  sbatch ${slurmjob_file}
  cd - > /dev/null

}


submit_one_job_per_file() {

    # Read arguments
    sample=$1
    model=$2
    file_list_directory=$3
    xsec_directory=$4
    xrdcp_files=$5
    use_custom_slurmjob_template=$6

    if [ ${use_custom_slurmjob_template} == "true" ]; then
        sample_flag="-s ${sample}"
    else
        sample_flag=""
    fi

    # Get sample cross-section, file list and number of files
    xsec=$(cat ${xsec_directory}/${model}_xsec.txt)
    file_with_file_list=${file_list_directory}/${model}_files.txt
    file_list=$(get_file_list ${file_with_file_list})
    number_of_files=$(get_number_of_files ${file_list})

    # Submit one job per file
    #for ((part=1; part<(${number_of_files}+1); part++)); do
    for ((part=1; part<2; part++)); do  # for tests
        echo ""

        # Make the slurm job file
        ifile=$((${part}-1))
        id=${model}_part-${part}
        input_file=$(get_file_number_from_file_list ${ifile} ${file_list[@]})
        output_directory=$(get_output_directory_for_sample ${sample})
        output_file=${output_directory}/PFNANOAODSKIM/PFNANOAODSKIM_${id}.root
        slurmjob_file=${id}

        # Check if output file already exists
        #redirector=$(echo ${output_file} | cut -d/ -f-4)
        redirector=$(get_redirector ${output_file})
        file_path=$(get_path ${output_file})
        
        if [ $(file_exist ${redirector} ${file_path}) -ne 0 ]; then

            make_slurmjob_file ${sample_flag} -j ${id} -l ${log_dir} -i ${input_file} -xsec ${xsec} -o ${output_file} -xrdcp ${xrdcp_files} -f ${slurmjob_dir}/${slurmjob_file}

            # Submit job
            submit_job ${slurmjob_dir} ${slurmjob_file}
            #cat ${slurmjob_dir}/${slurmjob_file}
        fi
    done
}


# Make slurm job files directory
slurmjob_dir=slurmjob_files
make_dir ${slurmjob_dir}

# Make log files directory
log_dir=${PWD}/logs
make_dir ${log_dir}


###################
###  t-channel  ###
###################

if [ "${sample}" == "tchannel" ]; then
    file_list_directory=${PWD}/samplesInformation/tchannel_106X_v02
    xsec_directory=${PWD}/samplesInformation/tchannel_106X_v02
    xrdcp_files=false
    use_custom_slurmjob_template=false
  
    # List of models to process
    mMedList=(  500 1000 1500 2000 3000)
    mDarkList=( 20  20   20   20   20  )
    rinvList=(  0.3 0.3  0.3  0.3  0.3 )
    alphaList=(peak peak peak peak peak)
    
    # Loop over all models
    for ((imodel=1; imodel<4; imodel++)); do
        echo ""
  
        # Get model parameters
        mMed=${mMedList[${imodel}]}
        mDark=${mDarkList[${imodel}]}
        rinv=${rinvList[${imodel}]}
        alpha=${alphaList[${imodel}]}

        # Define model name
        model=t-channel_mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-${alpha}_yukawa-1_13TeV-madgraphMLM-pythia8
    
        # Submit jobs
        submit_one_job_per_file ${sample} ${model} ${file_list_directory} ${xsec_directory} ${xrdcp_files} ${use_custom_slurmjob_template}
    done
fi


#############
###  QCD  ###
#############

if [ "${sample}" == "QCD" ]; then
    file_list_directory=${PWD}/samplesInformation/QCD_106X_v02/lite
    xsec_directory=${PWD}/samplesInformation/QCD_106X_v02
    xrdcp_files=true
    use_custom_slurmjob_template=true

    # List of models to process
    pt1List=(170 300 470 600 800  1000 1400 1800 2400 3200)
    pt2List=(300 470 600 800 1000 1400 1800 2400 3200 Inf )
    
    # Loop over all models
    for ((imodel=1; imodel<10; imodel++)); do
        echo ""

        # Get model parameters
        pt1=${pt1List[${imodel}]}
        pt2=${pt2List[${imodel}]}

        # Define model name
        model=QCD_pt_${pt1}to${pt2}

        # Submit jobs
        submit_one_job_per_file ${sample} ${model} ${file_list_directory} ${xsec_directory} ${xrdcp_files} ${use_custom_slurmjob_template}
    done
fi


###############
###  ttbar  ###
###############

if [ "${sample}" == "ttbar" ]; then
    file_list_directory=${PWD}/samplesInformation/ttbar_106X_v02/lite
    xsec_directory=${PWD}/samplesInformation/ttbar_106X_v02
    xrdcp_files=true
    use_custom_slurmjob_template=false

    # Define model name
    model=ttbar

    # Submit jobs
    submit_one_job_per_file ${sample} ${model} ${file_list_directory} ${xsec_directory} ${xrdcp_files} ${use_custom_slurmjob_template}
fi

