#!/bin/bash

directory=./lite
if [ ! -d ${directory} ]; then mkdir -p ${directory}; fi

pt_bins_edges=(170 300 470 600 800 1000 1400 1800 2400 3200 Inf)
nbins=$((${#pt_bins_edges[@]}-1))

for ((ibin=0; ibin<${nbins}; ibin++)); do
    pt1=${pt_bins_edges[${ibin}]}
    pt2=${pt_bins_edges[$((${ibin}+1))]}
    file=QCD_pt_${pt1}to${pt2}_files.txt
    echo "cathead -f ${file} -p 5 > ${directory}/${file}"
    cathead -f ${file} -p 5 > ${directory}/${file}
done
