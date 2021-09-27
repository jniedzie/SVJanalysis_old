#!/bin/bash

if [ -z $XROOTD_UTILITIES_SOURCED ]; then
XROOTD_UTILITIES_SOURCED=true

THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source ${THIS_DIRECTORY}/bashUtilities.sh


has_redirector() {

  if [ "$(get_first_characters $1 7)" == "root://" ]; then
      echo "true"  
  else
      echo "false"
  fi
}


get_redirector() {
  echo $1 | cut -d/ -f-4
}


get_path() {
  echo $1 | cut -d/ -f4-
}

fi
