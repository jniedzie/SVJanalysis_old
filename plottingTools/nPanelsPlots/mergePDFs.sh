#!/bin/bash

##################################    README   #################################
#
# Merge PDFs in ${figuresDir} into ${mergedPDF}.
# PDFs can be merged in alphabetic order, order="alpha", or by decreasing auc,
# order="auc".
# If order="auc", then the path to the auc csv file, ${aucFile}, must be
# provided. The auc file must have the following form:
#   variable, auc
#   MET_pt,   0.95
#   Jet_pt,   0.8
#   jet_mass, 0.7
#   MET_phi,  0.5
#
################################################################################


figuresDir=./figures/
mergedPDF=merged.pdf
order="auc"  # choices: "alpha", "auc"
aucFile=auc.csv

## Make list of PDFs...
PDFs=""

# ... in alphabetic order
if [ "${order}" == "alpha" ]; then
  for file in `ls ${figuresDir}`; do
    a=${#file}-4
    if [ "${file:a:4}" = ".pdf" ]; then
      PDFs=${PDFs}' '${figuresDir}${file}
    fi
  done

# ... in decreasing auc order
elif [ "${order}" == "auc" ]; then
  n=$(($(cat ${aucFile} | wc -l)-1))
  variables=`cat ${aucFile} | cut -d, -f1 | tail -n ${n}`
  for variable in ${variables[@]}; do
    PDFs=${PDFs}' '${figuresDir}${variable}.pdf
  done
fi

## Merge PDF
gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=${mergedPDF} -dBATCH ${PDFs} > /dev/null 2>&1

## Print out if merging was successful
if [ "$?" = "0" ]; then
  echo "${mergedPDF} has been created"
else
  echo "Could not create ${mergedPDF}"
fi
