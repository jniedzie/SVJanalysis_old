#!/bin/bash

source ../../make_files_list.sh

pf_nano_skim_file_list_directory=/t3home/fleble/SVJ/SVJanalysis/branchesProducer/batch/samplesInformation/ttbar_106X_v02/files/
file_name=ttbar_files.txt
pf_nano_skim_file_full_path=${pf_nano_skim_file_list_directory}/${file_name}

$(make_files_list ${pf_nano_skim_file_full_path} ${file_name})
