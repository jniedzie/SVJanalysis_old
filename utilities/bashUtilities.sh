#!/bin/bash

if [ -z $BASH_UTILITIES_SOURCED ]; then
BASH_UTILITIES_SOURCED=true

THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source ${THIS_DIRECTORY}/xrootdUtilities.sh


file_exists() {

    local file_name=$1
  
    # Check only when there is no file name expansion
    if [ $(echo ${file_name} | grep -Ev "{[0-9]+\.\.[0-9]+}") ]; then
        if [ "$(has_redirector ${file_name})" == "true" ]; then
            local redirector=$(get_redirector ${file_name})
            local file_path=$(get_path ${file_name})
      
            xrdfs ${redirector} ls ${file_path} > /dev/null 2>&1
            local returned_value=$?
          
            # When running CMSSW, xrdfs returns code 50 with the following error message:
            # [ERROR] Internal error
            # Try to locate the file on T3 PSI storage element instead in that case
            # if the redirector is the T3 PSI redirector
            if [ ${returned_value} -eq 50 ] && [ "${redirector}" == "root://t3se01.psi.ch:1094/" ]; then
                local path_to_storage_element=/pnfs/psi.ch/cms/trivcat
                ls ${path_to_storage_element}/${file_path} > /dev/null 2>&1
                local returned_value=$?
            fi
      
        else
            ls ${file_name} > /dev/null 2>&1
            local returned_value=$?
        fi
    else
        local returned_value=0
    fi
  
    if [ ${returned_value} -eq 0 ]; then echo true; else echo false; fi
}


make_dir() {
    ### Make directory given as argument if it does not exist ###

    directory=$1

    if [ "$(file_exists ${directory})" == "false" ]; then
        echo "Making directory ${directory}"
        if [ "$(has_redirector ${directory})" == "true" ]; then
            local redirector=$(get_redirector ${directory})
            local path=$(get_path ${directory})
            xrdfs ${redirector} mkdir -p ${path}
        else
            if [ ! -d ${directory} ]; then mkdir -p ${directory}; fi
        fi
    fi
}


get_first_characters() {

    local str=$1
    local n=$2
    local first_chars=${str:0:$n}
    echo ${first_chars}
}


get_last_characters() {

    local str=$1
    local n=$2
    local ichar=(${#str}-${n})
    local last_chars=${str:$ichar:$n}
    echo ${last_chars}
}

fi
