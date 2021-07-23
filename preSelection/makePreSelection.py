from coffea import processor
#from coffea.nanoevents import PFNanoAODSchema
from coffea.nanoevents import NanoAODSchema
import awkward as ak
import awkward0 as ak0
import uproot3
import numpy as np
import os
import time
import argparse
from collections import OrderedDict

import processorPreSelection


def print_cutflow(cutflow):
    """
    Print cut efficiencies for cuts needed for defining some of the variables.
    E.g. for deltaR between leading 2 jets, events needs to have at least 2 jets.
    """

    lenCol1 = max([ len(k) for k in cutflow.keys() ])
    efficiencies = OrderedDict()

    print("\nCutflow:")
    print("\tCut" + (lenCol1-3)*" " + "  Abs. eff. [%]   Rel. eff. [%]")

    nAll = cutflow["all"]
    for cut, n in cutflow.items():
        if cut != "all":
            absoluteEfficiency = 100*n/nAll
            if nPreviousCut > 0.:
                relativeEfficiency = 100*n/nPreviousCut
            else:
                relativeEfficiency = np.nan
            spaces = (lenCol1-len(cut)+(absoluteEfficiency<10))*" "
            print("\t%s%s  %.2f           %.2f" %(cut, spaces, absoluteEfficiency, relativeEfficiency))
            efficiencies[cut] = absoluteEfficiency/100
        nPreviousCut = n
    
    efficiencies["totalEfficiency"] = absoluteEfficiency
    
    return efficiencies


def delete_empty_fields(accumulator, debug):
    """ ... """

    deletion = False
    if debug: print("\nDeleting empty fields...")

    keys_to_delete = []
    for key, column in accumulator.items():
        if column.value.shape[0] == 0:
            keys_to_delete.append(key)
            if debug:
                deletion = True
                print("%s is empty. Deleting from accumulator." %key)

    if debug:
        if not deletion:
            print("No empty field in accumulator.")
        print("")
        
    for key in keys_to_delete: accumulator.pop(key)

    return
      

def make_events_branches(accumulator, debug):
    """Make branches for Events tree."""

    branches = {}
    branchesInit = {}
    lenKeys = []

    # Finding keys giving the length of jagged arrays
    for k, v in accumulator.items():
        nKey = "n"+str(k.split("_")[0])
        if (not k.startswith("n")) and (nKey in accumulator.keys()):
            lenKeys.append(nKey)

    # Making branches
    # Need to use ak0 because it is not yet implemented in uproot4 i.e. ak1 (Feb. 2021)
    for k, v in accumulator.items():
        if debug:
            print("%s: %s" %(k, ak.type(ak.Array(v.value))))
        if not k.startswith("n"):
            nKey = "n"+str(k.split("_")[0])
            if nKey in lenKeys:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
                # Case distinction for type of the jagged-array collections, treated as "object" in the processor
                # _candIdx and _jetIdx must be saved as integers ("i4") to use array at once syntax like 
                #     jetIdx = (events.JetPFCandsAK4_jetIdx == ijet)
                #     candIdx = events.JetPFCandsAK4_candIdx[jetIdx]
                #     PFcands_eta = events.JetPFCands_eta[candIdx]
                if k.endswith("_pFCandsIdx") or k.endswith("_jetIdx"):
                    branchesInit[k] = uproot3.newbranch(np.dtype("i4"), size=nKey)
                else:
                    branchesInit[k] = uproot3.newbranch(np.dtype("f8"), size=nKey)
            else:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
                branchesInit[k] = uproot3.newbranch(v.value.dtype)
        else:
            if k in lenKeys:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
            else:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
                branchesInit[k] = uproot3.newbranch(v.value.dtype)

    return branches, branchesInit


def make_efficiencies_branches(efficiencies):
    """Make branches for Efficiencies tree."""

    branches = {}
    branchesInit = {}

    for k, v in efficiencies.items():
        branchesInit[k] = np.dtype("f8")
        branches[k] = np.array([v])

    return branches, branchesInit


def write_root_file(accumulator, outputFile, efficiencies, debug):
    """
    2 trees are written:
       * Events: standard NTuple events tree
       * Cuts: a tree with one branch `Efficiency`, having only one leaf, representing pre-selection efficiency
    """

    # Making branches to write to Events tree
    branchesEvents, branchesInitEvents = make_events_branches(accumulator, debug)
    branchesEfficiencies, branchesInitEfficiencies = make_efficiencies_branches(efficiencies)

    if debug:
        print("\nEvents branches keys:")
        print(branchesEvents.keys())
        print("\nEvents branchesInit keys:")
        print(branchesInitEvents.keys())


    # Save branches to ROOT file
    # Need to use uproot3 because it is not implemented yet in uproot4 (Feb. 2021)
    with uproot3.recreate(outputFile) as f:
        f["Events"] = uproot3.newtree(branchesInitEvents)
        f["Events"].extend(branchesEvents)
        print("\nTTree Events saved to output file %s" %outputFile)

        # Save cut efficiency to ROOT file
        f["Efficiencies"] = uproot3.newtree(branchesInitEfficiencies)
        f["Efficiencies"].extend(branchesEfficiencies)
        print("TTree Efficiencies saved to output file %s" %outputFile)

    return


def main(inputFiles, outputFile, fileType, processorName, chunksize, maxchunks, nworkers, debug):

    print("Input files:")
    print(inputFiles)

    ## Fileset
    fileset = { "fileset": inputFiles }

    ## Make pre-selections
    output = processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = getattr(processorPreSelection, processorName)(fileType),
        executor = processor.iterative_executor,
        executor_args = {"schema": NanoAODSchema, "workers": nworkers},
        chunksize = chunksize,
        maxchunks = maxchunks
        )

    ## Print out cutflow
    cutflow = output.pop("cutflow")
    efficiencies = print_cutflow(cutflow)


    ## Delete empty fields
    delete_empty_fields(output, debug)

    ## Making output ROOT file
    if efficiencies["totalEfficiency"] == 0.:
        print("\nWARNING: 0 event passed pre-selections. Will not write an empty ROOT file.")
    else:
        write_root_file(output, outputFile, efficiencies, debug)


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

    ## Get list of input ROOT files
    # If coma separated list of ROOT files (which ends with the .root of the last ROOT file)
    if args.inputFiles.endswith(".root"):
        inputFiles = args.inputFiles.split(",")
    # Else list of ROOT files in txt file
    else:
        with open (args.inputFiles, "r") as txtfile:
            inputFiles = txtfile.readlines()
        inputFiles = [ x.replace("\n", "") for x in inputFiles ]
        inputFiles = [ x for x in inputFiles if x.endswith(".root") ]

    if args.debug:
        debug = True
    else:
        debug = False

    ## Make pre-selection
    main(inputFiles, args.output, args.fileType, args.processor, args.chunksize, args.maxchunks, args.nworkers, debug)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


