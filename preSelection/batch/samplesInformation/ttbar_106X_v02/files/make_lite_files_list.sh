#!/bin/bash

directory=./lite
if [ ! -d ${directory} ]; then mkdir -p ${directory}; fi

file=ttbar_files.txt
cathead -f ${file} -p 1 > ${directory}/${file}
