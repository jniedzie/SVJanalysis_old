#!/bin/bash

if [ -z $SLURM_JOB_UTILITIES_SOURCED ]; then
SLURM_JOB_UTILITIES_SOURCED=true


THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source ${THIS_DIRECTORY}/../xrootdUtilities.sh
THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source ${THIS_DIRECTORY}/../bashUtilities.sh


get_file_list() {

  local file_argument=$1

  if [ "$(get_last_characters ${file_argument} 4)" == ".txt" ]; then
     echo $(cat ${file_argument})
  else
     echo ${file_argument}
  fi
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
  echo $1 | cut -d/ -f-4
}


get_path() {
  echo $1 | cut -d/ -f4-
}


file_arg_exist() {
  # Return 0 if file exist, else 54

  local file_argument=$1

  # The file argument could countain a comma separated list of files
  if [ $(echo ${file_argument} | grep ".root,") ]; then
      local returned_value=0
      file_names=${file_argument//'.root,'/'.root '}
      for file_name in ${file_names}; do
          local returned_value=$((${returned_value} + $(file_exists ${file_name})))
      done
  else
      local returned_value=$(file_exists ${file_argument})
  fi
 
  echo ${returned_value}
}


make_flags() {

  local flag_string=$1

  echo ${flag_string} | tr : ' '
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
  python ${THIS_DIRECTORY}/make_slurmjob_file.py $@
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
    local slurmjob_template_name=$1
    local partition=$2
    local time_limit=$3
    local memory=$4
    local model=$5
    local file_argument=$6
    local output_directory=$7
    local data_tier=$8
    local file_prefix=$9
    local xrdcp_files=${10}
    local force_recreate=${11}
    local first_part_number=${12}
    local last_part_number=${13}
    local specific_flags=$(make_flags ${14})

    # Get sample cross-section, file list and number of files
    local file_list=$(get_file_list ${file_argument})
    local number_of_files=$(get_number_of_files ${file_list})

    # Define first and last part numbers
    if [ "${first_part_number}" == "-1" ]; then
        local first_part_number=1
    fi
    if [ "${last_part_number}" == "-1" ]; then
        local last_part_number=${number_of_files}
    fi

    # Submit one job per file
    for ((part=${first_part_number}; part<=${last_part_number}; part++)); do
    #for ((part=1; part<2; part++)); do  # for tests
        echo ""

        # Make the slurm job file
        local ifile=$((${part}-1))
        local id=${model}_part-${part}
        local input_file=$(get_file_number_from_file_list ${ifile} ${file_list[@]})
        local output_file=${output_directory}/${datatier}/${file_prefix}${id}.root
        local slurmjob_file=${id}

        # Check if input files exist
        if [ "$(file_arg_exist ${input_file})" == "true" ]; then

            # Check if output file already exists
            if [ "$(file_arg_exist ${output_file})" == "false" ] || [ "${force_recreate}" == "true" ]; then

                make_slurmjob_file -s ${slurmjob_template_name} -p ${partition} -t ${time_limit} -m ${memory} -j ${id} -l ${log_dir} -i ${input_file} -o ${output_file} -xrdcp ${xrdcp_files} ${specific_flags} -f ${slurmjob_dir}/${slurmjob_file}

                # Submit job
                submit_job ${slurmjob_dir} ${slurmjob_file}
                #cat ${slurmjob_dir}/${slurmjob_file}

            else
                echo "File ${output_file} already exists and force_recreate is false. Skipping part ${part}."

            fi

        else
            echo "Input file(s) ${input_file} do not exist. Skipping part ${part}."
        fi
    done
}

fi
