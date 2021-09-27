#!/bin/bash

make_files_list() {

    local pf_nano_skim_file_full_path=$1
    local file_name=$2

    > ${file_name}
    
    for pf_nano_skim_file_name in $(cat ${pf_nano_skim_file_full_path}); do
        echo ${pf_nano_skim_file_name//PFNANOAODSKIM/PFNANOAODEXTENSION},${pf_nano_skim_file_name} >> ${file_name}
    done
}
