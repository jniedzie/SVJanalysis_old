#!/bin/bash

figuresDir=./figures/
mergedPDF=merged.pdf

PDFs=""
for file in `ls ${figuresDir}`; do
  a=${#file}-4
  if [ "${file:a:4}" = ".pdf" ]; then
    PDFs=${PDFs}' '${figuresDir}${file}
  fi
done

gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=${mergedPDF} -dBATCH ${PDFs} > /dev/null 2>&1
if [ "$?" = "0" ]; then
  echo "${mergedPDF} has been created"
else
  echo "Could not creat ${mergedPDF}"
fi
