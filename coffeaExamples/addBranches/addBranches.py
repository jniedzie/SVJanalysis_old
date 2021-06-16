from coffea import processor
from coffea.nanoevents import BaseSchema, NanoAODSchema, PFNanoAODSchema
import awkward as ak
import uproot3
import numpy as np
import time
import argparse
import processorAddBranches


def make_input_files_list(input_files):

    # If coma separated list of ROOT files (which ends with the .root of the last ROOT file)
    if input_files.endswith(".root"):
        input_files = input_files.split(",")
    # Else list of ROOT files in txt file
    else:
        with open (input_files, "r") as txt_file:
            input_files = txt_file.readlines()
        input_files = [ x.replace("\n", "") for x in input_files ]
        input_files = [ x for x in input_files if x.endswith(".root") ]

    return input_files


def make_events_branches(accumulator, debug):
    """Make branches for Events tree."""

    branches = {}
    branches_init = {}
    len_keys = []

    # Finding keys giving the length of jagged arrays
    for k, v in accumulator.items():
        nkey = "n"+str(k.split("_")[0])
        if (not k.startswith("n")) and (nkey in accumulator.keys()):
            len_keys.append(nkey)

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
            if k in len_keys:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
            else:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
                branches_init[k] = uproot3.newbranch(v.value.dtype)

    return branches, branches_init


def write_root_file(accumulator, output_file, debug):
    """Write a test.root ROOT file with the new branches."""

    output_file = "test.root" # For test purposes

    # Making branches to write to Events tree
    branches, branches_init = make_events_branches(accumulator, debug)

    if debug:
        print("\nEvents branches keys:")
        print(branches.keys())
        print("\nEvents branches_init keys:")
        print(branches_init.keys())


    # Save branches to ROOT file
    # Need to use uproot3 because it is not implemented yet in uproot4 (Feb. 2021)
    with uproot3.recreate(output_file) as f:
        f["Events"] = uproot3.newtree(branches_init)
        f["Events"].extend(branches)
        print("\nTTree Events saved to output file %s" %output_file)

    return


def main(input_file, processor_name, schema, pf_nano_aod_version, chunksize, maxchunks, nworkers, debug):

    print("Input file:")
    print(input_file)

    ## Fileset
    fileset = { "fileset": [input_file] }

    print(processorAddBranches)

    ## Make pre-selections
    output = processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = getattr(processorAddBranches, processor_name)(schema, pf_nano_aod_version),
        executor = processor.iterative_executor,
        executor_args = {"schema": eval(schema+"Schema"), "workers": nworkers},
        chunksize = chunksize,
        maxchunks = maxchunks,
        )

    ## Making output ROOT file
    write_root_file(output, input_file, debug)



if __name__ == "__main__":

    tstart = time.time()

    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputfiles",
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
        "-p", "--processor",
        help="Coffea processor to be used, as defined in processor.py",
        required=True
        )
    parser.add_argument(
        "-s", "--schema",
        help="Coffea nanoevents schema",
        choices=["Base", "PFNanoAOD"],
        required=True
        )
    parser.add_argument(
        "-pv", "--pfnanoaodversion",
        help="PF nano AOD version",
        choices=["102X", "106X"],
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

    ## Get list of input ROOT files
    input_files = make_input_files_list(args.inputfiles)

    ## Make pre-selection
    for input_file in input_files:
        main(input_file, args.processor, args.schema, args.pfnanoaodversion, args.chunksize, args.maxchunks, args.nworkers, args.debug)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


