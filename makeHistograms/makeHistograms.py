import coffea as cf
import awkward as ak
import awkward0 as ak0
import uproot
import uproot3
import numpy as np
import time
import argparse
import os
import sys
import json

sys.path.append("../utilities/")
import utilities as utl
import processorHistograms as procHist


def print_cutflow(cutflow, histogramKeys):
    """
    Print cut efficiencies for cuts needed for defining some of the variables.
    E.g. for deltaR between leading 2 jets, events needs to have at least 2 jets.
    """

    cutflowKeys = []
    for key in cutflow.keys():
        if key == "noCut": continue
        if key not in histogramKeys:
            cutflowKeys.append(key)

    lenCol1 = max([ len(k) for k in cutflowKeys ])

    print("\nCutflow:")
    print("\tCut" + (lenCol1-3)*" " + "  Abs. eff. [%]")
    nAll = cutflow["noCut"]
    for cut in cutflowKeys:
        print("\t%s%s  %.2f" %(cut, (lenCol1-len(cut))*" ", 100*cutflow[cut]/nAll))

    return


def calculate_efficiency(fileset, efficiencyBranch, genWeightBranch="Events/genWeight"):
    """
    Calculate efficiency using the effienciency written in input ROOT files,
    weighted by the sum of generator weights for the events in this ROOT file.
    """

    efficiencies = []
    sumGenWeights = []

    for file_ in fileset:
        data = uproot.open(file_)
        efficiencies.append( data[efficiencyBranch].arrays()[efficiencyBranch.split("/")[-1]][0] )
        sumGenWeights.append( ak.sum(data[genWeightBranch].arrays()[genWeightBranch.split("/")[-1]]) )

    efficiency = sum([ eff*sgw for eff, sgw in zip(efficiencies, sumGenWeights)]) / sum(sumGenWeights)
    
    return efficiency


def write_root_file(accumulator, rootFileName, cutflow, lumi, xSection, efficiency):
    """
    Write histograms, stored in a coffea accumulator, to a ROOT file.
    Number of generated events (cutflow), pre-selection cut efficiency (efficiency), luminosity (lumi)
    and cross-section (xSection) are needed to normalize the histograms.
    """

    # Need to use uproot3 because it is not implemented yet in uproot4 (Feb. 2021)
    print("\nWill recreate ROOT file %s\n" %rootFileName)

    writtenHists = []
    with uproot3.recreate(rootFileName) as f:    # Create/update output file
        for variable, hist in accumulator.items():
            if isinstance(hist, cf.hist.Hist):
                hist.scale(lumi * xSection * efficiency / cutflow["noCut"])
                f[variable] = cf.hist.export1d(hist)
                writtenHists.append(variable)

    print("Histograms written to output ROOT file:")
    for writtenHist in writtenHists:
        print("\t%s" %writtenHist)

    return


def main(binning, sample, fileType, processor, chunksize, maxchunks, nworkers, outputDirectory, lumi, useEfficiencies, efficiency):

    ## Get sample name and cross-section
    sampleName = sample["name"]
    xSection = sample["XSection"]


    ## Fileset
    fileset = { sampleName: sample["fileset"] }


    ## Make histograms
    output = cf.processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = getattr(procHist, processor)(binning, fileType),
        executor = cf.processor.iterative_executor,
        executor_args = {"schema": cf.nanoevents.BaseSchema, "workers": nworkers},
        chunksize = chunksize,
        maxchunks = maxchunks,
    )

    ## Efficiencies
    if useEfficiencies:
        if isinstance(efficiency, str):
            efficiency = calculate_efficiency(sample["fileset"], efficiency)
        # else efficiency is a float that already corresponds to the efficiency
    else:
        efficiency = 1.


    ## Print out cutflow
    cutflow = output.pop("cutflow")
    print_cutflow(cutflow, output.keys())


    ## Save histograms to ROOT file
    rootFileName = outputDirectory + sampleName + ".root"
    write_root_file(output, rootFileName, cutflow, lumi, xSection, efficiency)

    return


if __name__ == "__main__":

    tstart = time.time()

    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-b", "--binning",
        help="json file describing binning of the histograms",
        required=True
        )
    parser.add_argument(
        "-sd", "--samplesDescription",
        help="json file describing samples location",
        required=True
        )
    parser.add_argument(
        "-s", "--samples",
        help="Comma separated list samples to pick up among the samples described in the sample file",
        required=True
        )
    parser.add_argument(
        "-t", "--fileType",
        help="Input file type. If not given, then the tyoe will be infered. Choices = PFnano102X, PFnano106",
        choices=["PFnano102X", "PFnano106X"],
        )
    parser.add_argument(
        "-p", "--processor",
        help="Coffea processor to be used, as defined in processorHistogram.py",
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
        "-o", "--outputDirectory",
        help="Path to the directory where to recreate/update ROOT file (default=./)",
        default="./"
        )
    parser.add_argument(
        "-l", "--lumi",
        help="Total luminosity for normalization of the histograms (default=59725 pb-1)",
        default=59725.0,   # 21071.0+38654.0
        )
    parser.add_argument(
        "-e", "--efficiency",
        help="Use efficiencies from the input files",
        nargs='?', const="Efficiencies/totalEfficiency", default="False",
        )

    args = parser.parse_args()

    
    ## Create output directory if it does not exist
    outputDirectory = args.outputDirectory
    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)

    ## All samples description
    samplesDescription = utl.makeJsonData(args.samplesDescription)

    ## Define the binning of the different variables to histogram
    with open(args.binning, 'r') as f:
        binning = json.load(f)["binning"]

    ## Get samples for which to make histograms
    if not args.samples:
        samplesNames = list(samplesDescription.keys())
    else:
        samplesNames = args.samples.split(",")

    ## Efficiencies
    if args.efficiency == "False":
        useEfficiencies = False
        efficiency = None
    else:
        useEfficiencies = True
        if args.efficiency.replace(".", "").isdigit():
            efficiency = float(args.efficiency)
        else:
            efficiency = args.efficiency

    ## Loop over all samples
    for sampleName in samplesNames:
        sample = samplesDescription[sampleName]
        if "name" not in sample.keys():
            sample["name"] = sampleName

        # Get / infer file type
        if not args.fileType:
            # Open 1st file and check some branches to infer the type
            file_ = uproot3.open(sample["fileset"][0])
            events = file_["Events"]  # Assuming the events TTree is called Events
            branches = [ branch.decode("utf-8") for branch in events.keys() ]

            if "JetPFCandsAK4_jetIdx" in branches:
                fileType = "PFnano106X"
            elif "JetPFCands_jetIdx" in branches:
                fileType = "PFnano102X"
            else:
                print("ERROR: file type cannot be inferred!")
                sys.exit()
            print("Inferred file type: %s" %fileType)

        else:
            fileType = args.fileType

        # Make histograms
        main(binning, sample, fileType, args.processor, args.chunksize, args.maxchunks, args.nworkers, outputDirectory, float(args.lumi), useEfficiencies, efficiency)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


