from coffea import processor
from coffea.nanoevents import BaseSchema
import awkward as ak
import awkward0 as ak0
import uproot3
import numpy as np
import os
import time
import argparse

from processorPreSelection import Preselection
#from preSelectionProcessor_lite import Preselection


def main(inputFiles, outputFile, fileType, chunksize, maxchunks, nworkers):

    print("Input files:")
    print(inputFiles)

    ## Fileset
    fileset = { "fileset": inputFiles }


    ## Make pre-selections
    output = processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = Preselection(fileType),
        executor = processor.iterative_executor,
        executor_args = {"schema": BaseSchema, "workers": nworkers},
        chunksize = chunksize,
        maxchunks = maxchunks
        )


    ## Print out cutflow
    cutflow = output.pop("cutflow")
    lenCol1 = max([ len(k) for k in cutflow.keys() ])

    print("\nCutflow:")
    print("\tCut" + (lenCol1-3)*" " + "  Abs. eff. [%]   Rel. eff. [%]")
    nAll = cutflow["all"].value
    for cut, n in cutflow.items():
        if cut != "all":
            print("\t%s%s  %.2f           %.2f" %(cut, (lenCol1-len(cut))*" ", 100*n.value/nAll, 100*n.value/nPreviousCut))
        nPreviousCut = n.value
    totalEfficiency = n.value/nAll


    ## Making output ROOT file
    branches = {}
    branches_init = {}
    lenKeys = []

    # Finding keys giving the length of jagged arrays
    for k, v in output.items():
        nKey = "n"+str(k.split("_")[0])
        if (not k.startswith("n")) and (nKey in output.keys()):
            lenKeys.append(nKey)

    # Making branches
    # Need to use ak0 because it is not yet implemented in uproot4 i.e. ak1 (Feb. 2021)
    for k, v in output.items():
        if not k.startswith("n"):
            nKey = "n"+str(k.split("_")[0])
            if nKey in lenKeys:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
                # Case distinction for type of the jagged-array collections, treated as "object" in the processor
                # _candIdx and _jetIdx must be saved as integers ("i4") to use array at once syntax like 
                #     jetIdx = (events.JetPFCandsAK4_jetIdx == ijet)
                #     candIdx = events.JetPFCandsAK4_candIdx[jetIdx]
                #     PFcands_eta = events.JetPFCands_eta[candIdx]
                if k.endswith("_candIdx") or k.endswith("_jetIdx"):
                    branches_init[k] = uproot3.newbranch(np.dtype("i4"), size=nKey)
                else:
                    branches_init[k] = uproot3.newbranch(np.dtype("f8"), size=nKey)
            else:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
                branches_init[k] = uproot3.newbranch(v.value.dtype)
        else:
            if k in lenKeys:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
            else:
                branches[k] = ak.to_awkward0(ak.Array(v.value))
                branches_init[k] = uproot3.newbranch(v.value.dtype)

    # Save branches to ROOT file
    # Need to use uproot3 because it is not implemented yet in uproot4 (Feb. 2021)
    with uproot3.recreate(outputFile) as f:
        f["Events"] = uproot3.newtree(branches_init)
        f["Events"].extend(branches)
        print("\nTTree Events saved to output file %s" %outputFile)

        # Save cut efficiency to ROOT file
        f["Cuts"] = uproot3.newtree({"Efficiency": np.dtype("f8")})
        f["Cuts"].extend({"Efficiency": np.array([totalEfficiency])})
        print("TTree Cuts saved to output file %s" %outputFile)


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
        choices=["PFnano102X", "PFnano106X"],
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


    ## Make pre-selection
    main(inputFiles, args.output, args.fileType, args.chunksize, args.maxchunks, args.nworkers)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


