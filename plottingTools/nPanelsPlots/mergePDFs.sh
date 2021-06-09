#!/bin/bash

##################################    README   #################################
#
# Merge pdfs in ${figures_dir} into ${merged_pdf}.
# pdfs can be merged in alphabetic order, order="alpha", or by decreasing auc,
# order="auc".
# If order="auc", then the path to the auc csv file, ${auc_file}, must be
# provided. The auc file must have the following form:
#   variable, auc
#   MET_pt,   0.95
#   Jet_pt,   0.8
#   jet_mass, 0.7
#   MET_phi,  0.5
#
################################################################################


figures_dir=./figures/
merged_pdf=merged.pdf
order="auc"  # choices: "alpha", "auc"
auc_file=auc.csv

## Make list of pdfs...
pdfs=""

# ... in alphabetic order
if [ "${order}" == "alpha" ]; then
  for file in `ls ${figures_dir}`; do
    a=${#file}-4
    if [ "${file:a:4}" = ".pdf" ]; then
      pdfs=${pdfs}' '${figures_dir}${file}
    fi
  done

# ... in decreasing auc order
elif [ "${order}" == "auc" ]; then
  n=$(($(cat ${auc_file} | wc -l)-1))
  variables=`cat ${auc_file} | cut -d, -f1 | tail -n ${n}`
  for variable in ${variables[@]}; do
    pdfs=${pdfs}' '${figures_dir}${variable}.pdf
  done
fi

## Merge PDF
gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=${merged_pdf} -dBATCH ${pdfs} > /dev/null 2>&1

## Print out if merging was successful
if [ "$?" = "0" ]; then
  echo "${merged_pdf} has been created"
else
  echo "Could not create ${merged_pdf}"
fi
