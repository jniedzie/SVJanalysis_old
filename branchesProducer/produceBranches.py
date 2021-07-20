from coffea import processor
from coffea.nanoevents import PFNanoAODSchema
import awkward as ak
import uproot3
import numpy as np
import time
import argparse

from BranchesProducer import BranchesProducer


# Note: This function could be put in a separate file to be used by other scripts as well
def make_events_branches(accumulator, debug):
    """Make branches for Events tree."""

    branches = {}
    branches_init = {}
    len_keys = []

    # Finding keys giving the length of jagged arrays
    for k, v in accumulator.items():
        nkey = "n"+str(k.split("_")[0])
        if (not k.startswith("n")) and (nkey in accumulator.keys()) and (nkey not in len_keys):
            len_keys.append(nkey)

    if debug:
        print("Len keys: ", len_keys)

    # Making branches
    # Need to use ak0 because it is not yet implemented in uproot4 i.e. ak1 (Feb. 2021)
    for k, v in accumulator.items():
        if debug:
            print("%s: %s" %(k, ak.type(ak.Array(v.value))))
        if not k.startswith("n"):
            nkey = "n"+str(k.split("_")[0])
            if nkey in len_keys:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
                # Case distinction for type of the jagged-array collections, treated as "object" in the processor
                # _candIdx and _jetIdx must be saved as integers ("i4") to use array at once syntax like 
                #     jetIdx = (events.JetPFCandsAK4_jetIdx == ijet)
                #     candIdx = events.JetPFCandsAK4_candIdx[jetIdx]
                #     PFcands_eta = events.JetPFCands_eta[candIdx]
                if k.endswith("_candIdx") or k.endswith("_jetIdx"):
                    branches_init[k] = uproot3.newbranch(np.dtype("i4"), size=nkey)
                else:
                    branches_init[k] = uproot3.newbranch(np.dtype("f8"), size=nkey)
            else:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
                branches_init[k] = uproot3.newbranch(v.value.dtype)
        else:
            branches[k] = ak.to_awkward0(ak.Array(v.value))
            if k not in len_keys:
                branches_init[k] = uproot3.newbranch(v.value.dtype)

    return branches, branches_init


def write_root_file(accumulator, output_file_name, debug):
    """Write histograms, stored in a coffea accumulator, to a ROOT file.

    Args:
        accumulator (coffea.processor.dict_accumulator)
        output_file_name (str): Name of the ROOT file to create
        debug (bool)

    Returns:
        None
    """

    # Making branches to write to Events tree
    branches, branches_init = make_events_branches(accumulator, debug)

    if debug:
        print("\nEvents branches keys:")
        print(branches.keys())
        print("\nEvents branches_init keys:")
        print(branches_init.keys())

    # Save branches to ROOT file
    # Need to use uproot3 because it is not implemented yet in uproot4 (Feb. 2021)
    # Need to use compression=None because there is a bug with EFPs when the file is compressed:
    # More info at https://github.com/scikit-hep/uproot3/issues/506
    with uproot3.recreate(output_file_name, compression=None) as f:
        f["Events"] = uproot3.newtree(branches_init)
        f["Events"].extend(branches)
        print("\nTTree Events saved to output file %s" %output_file_name)

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
        help="Output ROOT file name (default=inputfile_extension.root)",
        )
    parser.add_argument(
        "-c", "--chunksize",
        help="Size of the data chunks (default=100000)",
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
        help="Number of worker nodes (default=4)",
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


