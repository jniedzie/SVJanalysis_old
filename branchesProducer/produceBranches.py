import sys
import time
import argparse

import numpy as np
import awkward as ak
import uproot3
from coffea import processor
from coffea.nanoevents import PFNanoAODSchema

sys.path.append("../utilities/")
import uproot3Utilities as uproot3utl
from BranchesProducer import BranchesProducer


def write_root_file(accumulator, output_file_name, debug, tree_name="Events"):
    """Write histograms, stored in a coffea accumulator, to a ROOT file.

    Args:
        accumulator (coffea.processor.dict_accumulator)
        output_file_name (str): Name of the ROOT file to create
        debug (bool)

    Returns:
        None
    """

    # Need to use compression=None because there is a bug with EFPs when the file is compressed:
    # More info at https://github.com/scikit-hep/uproot3/issues/506
    with uproot3.recreate(output_file_name, compression=None) as root_file:
        # Making branches to write to tree
        branches, branches_init = uproot3utl.make_branches(accumulator, debug)

        if debug:
            print("\n%s branches keys:" %tree_name)
            print(branches.keys())
            print("\n%s branches_init keys:" %tree_name)
            print(branches_init.keys())

        # Save branches to ROOT file
        uproot3utl.write_tree_to_root_file(root_file, tree_name, branches, branches_init)
        print("\nTTree %s saved to output file %s" %(tree_name, output_file_name))

    return


def main(input_file, output_file, chunksize, maxchunks, nworkers, debug):

    print("Input file:")
    print(input_file)

    ## Fileset
    fileset = { "fileset": [input_file] }

    ## Make pre-selections
    output = processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = BranchesProducer(),
        executor = processor.iterative_executor,
        executor_args = {"schema": PFNanoAODSchema, "workers": nworkers},
        chunksize = chunksize,
        maxchunks = maxchunks,
        )

    ## Making output ROOT file
    # if output not empty:
    if output_file is None:
        output_file = input_file[:-5] + "_extension.root"
    write_root_file(output, output_file, debug)



if __name__ == "__main__":

    tstart = time.time()

    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputfile",
        help="Input ROOT file used to make new branches",
        required=True,
        )
    parser.add_argument(
        "-o", "--outputfile",
        help="Output ROOT file name",
        )
    parser.add_argument(
        "-c", "--chunksize",
        help="Size of the data chunks (default=%(default)s)",
        default=100000,
        type=int,
        )
    parser.add_argument(
        "-m", "--maxchunks",
        help="Maximum number of chunks to process, no flag means no maximum",
        type=int,
        )
    parser.add_argument(
        "-n", "--nworkers",
        help="Number of worker nodes (default=%(default)s)",
        default=4, 
        type=int,
        )
    parser.add_argument(
        "-debug", "--debug",
        help="Debug mode",
        action="store_true",
        )
    
    
    args = parser.parse_args()

    main(args.inputfile, args.outputfile, args.chunksize, args.maxchunks, args.nworkers, args.debug)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


