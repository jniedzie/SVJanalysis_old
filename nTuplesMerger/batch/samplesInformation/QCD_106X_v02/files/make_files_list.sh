#!/bin/bash

source ../../make_files_list.sh

pf_nano_skim_file_list_directory=/t3home/fleble/SVJ/SVJanalysis/branchesProducer/batch/samplesInformation/QCD_106X_v02/files/

pt_bins_edges=(170 300 470 600 800 1000 1400 1800 2400 3200 Inf)
nbins=$((${#pt_bins_edges[@]}-1))

for ((ibin=0; ibin<${nbins}; ibin++)); do
    pt1=${pt_bins_edges[${ibin}]}
    pt2=${pt_bins_edges[$((${ibin}+1))]}
    file_name=QCD_pt_${pt1}to${pt2}_files.txt
    pf_nano_skim_file_full_path=${pf_nano_skim_file_list_directory}/${file_name}

    $(make_files_list ${pf_nano_skim_file_full_path} ${file_name})
done
