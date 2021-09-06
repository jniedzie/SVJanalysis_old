from coffea import processor
from coffea.nanoevents import NanoAODSchema
import uproot3
import numpy as np
import sys
import os
import time
import argparse
from collections import OrderedDict

sys.path.append("../utilities/")
import utilities as utl
import uproot3Utilities as uproot3utl
import processorPreSelection


def calculate_cut_efficiencies(cutflow):
    """Print cut efficiencies for cuts needed for defining some of the variables.

    For instance, to compute deltaR between leading 2 jets, events needs
    to have at least 2 jets.

    Args:
        cutflow (dict[coffea.processor.defaultdict_accumulator]):
            Keys are cut names
            Values are numbers of processed events (sum of gen weights)
 
    Returns:
        dict[float]
            Keys are cut names
            Values are absolute efficiencies
    """

    cut_names = []
    sum_gen_weights = []
    absolute_efficiencies = []
    relative_efficiencies = []

    sum_gen_weight_no_cut = cutflow["NoCut"]
    for cut_name, sum_gen_weight in cutflow.items():
        cut_names.append(cut_name)
        sum_gen_weights.append(sum_gen_weight)
        if cut_name == "NoCut":
            absolute_efficiency = 1.
            relative_efficiency = 1.
        else:
            absolute_efficiency = sum_gen_weight/sum_gen_weight_no_cut
            previous_cut_absolute_efficiency = absolute_efficiencies[-1]
            if previous_cut_absolute_efficiency > 0:
                relative_efficiency = absolute_efficiency / previous_cut_absolute_efficiency
            else:
                relative_efficiency = 0

        absolute_efficiencies.append(absolute_efficiency)
        relative_efficiencies.append(relative_efficiency)


    # Store information for all cuts (this is a repetition to access this information with convenience)
    cut_names.append("AllCuts")
    sum_gen_weights.append(sum_gen_weight)
    absolute_efficiencies.append(absolute_efficiency)
    relative_efficiencies.append(relative_efficiency)

    efficiencies = {
        "Names": cut_names,
        "SumGenWeights": sum_gen_weights,
        "Absolute": absolute_efficiencies,
        "Relative": relative_efficiencies,
    }

    return efficiencies


def print_cutflow(cut_names, absolute_efficiencies, relative_efficiencies):
    """Print cut efficiencies for cuts needed for defining some of the variables.

    For instance, to compute deltaR between leading 2 jets, events needs
    to have at least 2 jets.

    Args:
        cutflow (dict[coffea.processor.defaultdict_accumulator]):
            Keys are cut names
            Values are numbers of processed events (sum of gen weights)
 
    Returns:
        dict[float]
            Keys are cut names
            Values are absolute efficiencies
    """

    len_column1 = max([ len(k) for k in cut_names])

    print("\nCutflow:")
    print("\tCut" + (len_column1-3)*" " + "  Abs. eff. [%]   Rel. eff. [%]")
    for cut_name, absolute_efficiency, relative_efficiency in zip(cut_names, absolute_efficiencies, relative_efficiencies):
        absolute_efficiency = 100*absolute_efficiency
        relative_efficiency = 100*relative_efficiency
        spaces = (len_column1 - len(cut_name) + (absolute_efficiency<10))*" "
        spaces2 = (11 - (absolute_efficiency==100))*" "
        print("\t%s%s  %.2f%s%.2f" %(cut_name, spaces, absolute_efficiency, spaces2, relative_efficiency))
    print("")
 
    return


def make_cutflow_info(efficiency_dict):

    cutflow_info = { cut_name: sum_gen_weight for cut_name, sum_gen_weight in zip(efficiency_dict["Names"], efficiency_dict["SumGenWeights"]) }
    return cutflow_info


def make_metadata_info(total_cut_efficiency, cross_section, input_file_type):

    file_type_to_float = {
        "PFNanoAOD_106X_v02": 0.10602,
        "Delphes": 1.,
    }

    metadata = {
        "TotalCutEfficiency": total_cut_efficiency,
        "GenCrossSection": cross_section,
        "OriginalFormat": file_type_to_float[input_file_type],
    }

    return metadata


def write_root_file(output_file_name, events, cutflow, metadata, debug):
    """Write ROOT file with Events and Efficiencies TTrees.

    Args:
        output_file_name (str)
        events (coffea.processor.accumulator.dict_accumulator)
            Keys are Events branch names
            Values are branches (coffea.processor.accumulator.column_accumulator)
        efficiencies (dict[float])
            Keys are Efficicies branch names
            Values are efficiencies (float)
        debug (bool)

    Returns:
        None
    """

    # Making branches to write to Events tree
    branches = {}
    branches_init = {}
    branches["Events"], branches_init["Events"] = uproot3utl.make_branches(events, debug)
    branches["Cutflow"], branches_init["Cutflow"] = uproot3utl.make_branches(cutflow, debug)
    branches["Metadata"], branches_init["Metadata"] = uproot3utl.make_branches(metadata, debug)

    if debug:
        print("\nEvents branches keys:")
        print(branches["Events"].keys())
        print("\nEvents branchesInit keys:")
        print(branches_init["Events"].keys())


    print("\nWriting down output ROOT file %s" %output_file_name)
    with uproot3.recreate(output_file_name) as output_file:
        for tree_name in branches.keys():
            # Save tree to output ROOT file
            uproot3utl.write_tree_to_root_file(output_file, tree_name, branches[tree_name], branches_init[tree_name])
            print("TTree %s saved to output file" %tree_name)

        ## Save cutflow to output ROOT file
        #tree_name = "Cutflow"
        #uproot3utl.write_tree_to_root_file(root_file, tree_name, branches_cutflow, branches_init_cutflow)
        #print("TTree %s saved to output file %s" %(tree_name, output_file))

        ## Save some metadata to output ROOT file
        #tree_name = "Metadata"
        #uproot3utl.write_tree_to_root_file(root_file, tree_name, branches_metadata, branches_init_metadata)
        #print("TTree %s saved to output file %s" %(tree_name, output_file))

    return


def main(input_files, file_type, cross_section, output_file, processor_name, chunksize, maxchunks, nworkers, debug):
    """Produce ROOT file with pre-selected events.
    
    Args:
        input_files (list[str])
        output_file (str)
        file_type (str)
        processor_name (str)
        chunksize (int)
        maxchunks (int or None)
        nworkers (int)
        debug (bool)

    Returns:
        None
    """

    print("Input files:")
    print(input_files)

    # Fileset
    fileset = { "fileset": input_files }

    # Make pre-selections
    events = processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = getattr(processorPreSelection, processor_name)(file_type),
        executor = processor.iterative_executor,
        executor_args = {"schema": NanoAODSchema, "workers": nworkers},
        chunksize = chunksize,
        maxchunks = maxchunks
        )

    # Print out cutflow and make cutflow info
    cutflow_accumulator = events.pop("cutflow")
    efficiency_dict = calculate_cut_efficiencies(cutflow_accumulator)
    print_cutflow(efficiency_dict["Names"], efficiency_dict["Absolute"], efficiency_dict["Relative"])
    cutflow = make_cutflow_info(efficiency_dict)
    
    # Make matadata info
    total_cut_efficiency = efficiency_dict["Absolute"][efficiency_dict["Names"].index("AllCuts")]
    metadata = make_metadata_info(total_cut_efficiency, cross_section, file_type)

    # Write output ROOT file
    if total_cut_efficiency == 0.:
        print("\nWARNING: 0 event passed pre-selections. Will not write an empty ROOT file.")
    else:
        write_root_file(output_file, events, cutflow, metadata, debug)
    return


if __name__ == "__main__":

    tstart = time.time()

    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputFiles",
        help="Input ROOT files to skim. Format: Comma separated list of files "\
             +"or txt file with list of file (1 file name per line)",
        # Example for the format:
        # Ex1: file1.root,file2.root
        # Ex2: $ cat input.txt
        #      file1.root
        #      file2.root
        required=True
        )
    parser.add_argument(
        "-o", "--output",
        help="Output ROOT file",
        required=True
        )
    parser.add_argument(
        "-t", "--fileType",
        help="Input file type, mandatory argument",
        choices=["PFNanoAOD_102X", "PFNanoAOD_106X_v01", "PFNanoAOD_106X_v02"],
        required=True
        )
    parser.add_argument(
        "-xsec", "--genCrossSection",
        help="Sample cross section before pre-selection, mandatory argument",
        type=float,
        required=True
        )
    parser.add_argument(
        "-p", "--processor",
        help="Coffea processor to be used, as defined in processorPreSelection.py",
        required=True
        )
    parser.add_argument(
        "-c", "--chunksize",
        help="Size of the data chunks (default=100000)",
        default=100000, type=int
        )
    parser.add_argument(
        "-m", "--maxchunks",
        help="Maximum number of chunks to process, no flag means no maximum",
        type=int
        )
    parser.add_argument(
        "-n", "--nworkers",
        help="Number of worker nodes (default=4)",
        default=4, type=int
        )
    parser.add_argument(
        "-debug", "--debug",
        help="Debug mode",
        action="store_true",
        )
    
    
    args = parser.parse_args()

    ## Create output directory if it does not exist
    split = args.output.split("/")
    if len(split) > 1:
        if len(split) == 2:
            outputDirectory = (args.output[0]=="/")*"/" + split[0]
        elif len(split) > 1:
            outputDirectory = (args.output[0]=="/")*"/" + os.path.join(*args.output.split("/")[:-1])
        if not os.path.exists(outputDirectory):
            os.makedirs(outputDirectory)

    input_files = utl.make_file_list(args.inputFiles)

    ## Make pre-selection
    main(input_files, args.fileType, args.genCrossSection, args.output, args.processor, args.chunksize, args.maxchunks, args.nworkers, args.debug)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


