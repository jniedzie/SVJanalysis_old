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


def print_cutflow(cutflow):
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

    len_column1 = max([ len(k) for k in cutflow.keys() ])
    efficiencies = {}

    print("\nCutflow:")
    print("\tCut" + (len_column1-3)*" " + "  Abs. eff. [%]   Rel. eff. [%]")
    n_all = cutflow["all"]
    for cut, n in cutflow.items():
        if cut != "all":
            absolute_efficiency = 100*n/n_all
            if n_previous_cut > 0.:
                relative_efficiency = 100*n/n_previous_cut
            else:
                relative_efficiency = np.nan
            spaces = (len_column1-len(cut)+(absolute_efficiency<10))*" "
            print("\t%s%s  %.2f           %.2f" %(cut, spaces, absolute_efficiency, relative_efficiency))
            efficiencies[cut] = absolute_efficiency/100
        n_previous_cut = n
    
    efficiencies["totalEfficiency"] = absolute_efficiency
 
    return efficiencies


def write_root_file(output_file, events, efficiencies, debug):
    """Write ROOT file with Events and Efficiencies TTrees.

    Args:
        output_file (str)
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
    branches_events, branches_init_events = uproot3utl.make_branches(events, debug)
    branches_efficiencies, branches_init_efficiencies = uproot3utl.make_branches(efficiencies, debug)

    if debug:
        print("\nEvents branches keys:")
        print(branches_events.keys())
        print("\nEvents branchesInit keys:")
        print(branches_init_events.keys())


    with uproot3.recreate(output_file) as root_file:
        # Save events to output ROOT file
        uproot3utl.write_tree_to_root_file(root_file, "Events", branches_events, branches_init_events)
        print("\nTTree Events saved to output file %s" %output_file)

        # Save cut efficiencies to output ROOT file
        uproot3utl.write_tree_to_root_file(root_file, "Efficiencies", branches_efficiencies, branches_init_efficiencies)
        print("TTree Efficiencies saved to output file %s" %output_file)

    return


def main(input_files, output_file, file_type, processor_name, chunksize, maxchunks, nworkers, debug):
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

    ## Fileset
    fileset = { "fileset": input_files }

    ## Make pre-selections
    accumulator = processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = getattr(processorPreSelection, processor_name)(file_type),
        executor = processor.iterative_executor,
        executor_args = {"schema": NanoAODSchema, "workers": nworkers},
        chunksize = chunksize,
        maxchunks = maxchunks
        )

    ## Print out cutflow
    cutflow = accumulator.pop("cutflow")
    efficiencies = print_cutflow(cutflow)

    ## Making output ROOT file
    if efficiencies["totalEfficiency"] == 0.:
        print("\nWARNING: 0 event passed pre-selections. Will not write an empty ROOT file.")
    else:
        write_root_file(output_file, accumulator, efficiencies, debug)

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
    main(input_files, args.output, args.fileType, args.processor, args.chunksize, args.maxchunks, args.nworkers, args.debug)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


