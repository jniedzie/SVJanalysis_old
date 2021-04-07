import coffea as cf
import awkward as ak
import awkward0 as ak0
import uproot3
import numpy as np
import time
import argparse
import os
import sys
import json

sys.path.append("../pythonUtils/")
import utilities as utl
import processorHistogram as procHist


def print_cutflow(cutflow, histogramKeys):
    """
    Print cut efficiencies for cuts needed for defining some of the variables.
    E.g. for deltaR between leading 2 jets, events needs to have at least 2 jets.
    """

    cutflowKeys = []
    for key in cutflow.keys():
        if key == "all": continue
        if key not in histogramKeys:
            cutflowKeys.append(key)

    lenCol1 = max([ len(k) for k in cutflowKeys ])

    print("\nCutflow:")
    print("\tCut" + (lenCol1-3)*" " + "  Abs. eff. [%]")
    nAll = cutflow["noCut"]
    for cut in cutflowKeys:
        print("\t%s%s  %.2f" %(cut, (lenCol1-len(cut))*" ", 100*cutflow[cut]/nAll))

    return


def write_root_file(accumulator, rootFileName, mode, cutflow, lumi, xSection, processor):
    """
    Write histograms, stored in a coffea accumulator, to a ROOT file.
    Cut efficiencies for the variables histogrammed (cutflow), luminosity (lumi)
    and cross-section (xSection) are needed to normalize the histograms.
    """

    # Need to use uproot3 because it is not implemented yet in uproot4 (Feb. 2021)
    print("\nWill %s ROOT file %s.\n" %(mode.lower(), rootFileName))

    writtenHists = []
    with getattr(uproot3, mode)(rootFileName) as f:    # Create/update output file
        for variable, hist in accumulator.items():
            if isinstance(hist, cf.hist.Hist):
                if variable not in cutflow.keys():
                    print("WARNING: Sum of gen weights not available for %s." %variable)
                    print("         Histogram cannot to be written into ROOT file.")
                    print("         Check processor %s." %processor)
                else:
                    if cutflow[variable] > 0.:
                        hist.scale(lumi * xSection / cutflow[variable])
                        f[variable] = cf.hist.export1d(hist)
                        writtenHists.append(variable)
                    else:
                        print("WARNING: Sum of gen weights is 0 for %s." %variable)
                        print("         Histogram cannot to be written into ROOT file.")

    print("Histograms written to output ROOT file:")
    for writtenHist in writtenHists:
        print("\t%s" %writtenHist)

    return


def main(mode, binning, sample, processor, outputDirectory, lumi):

    ## Get sample name and cross-section
    sampleName = sample["name"]
    xSection = sample["XSection"]


    ## Fileset
    fileset = { sampleName: sample["fileset"] }


    ## Make histograms
    output = cf.processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = getattr(procHist, processor)(binning),
        executor = cf.processor.iterative_executor,
        executor_args = {"schema": cf.nanoevents.BaseSchema, "workers": 4},
        chunksize = 100000,
        maxchunks = None,
    )


    ## Print out cutflow
    cutflow = output.pop("cutflow")
    print_cutflow(cutflow, output.keys())


    ## Save histograms to ROOT file
    rootFileName = outputDirectory + sampleName + ".root"
    write_root_file(output, rootFileName, mode, cutflow, lumi, xSection, processor)

    return


if __name__ == "__main__":

    tstart = time.time()

    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--mode",
        help="Mode in which to open the output ROOT file. Choices: recreate (default), update",
        choices=["recreate", "update"], default="recreate"
        )
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
        "-p", "--processor",
        help="Coffea processor to be used, as defined in processorHistogram.py",
        required=True
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

    ## Loop over all samples
    for sampleName in samplesNames:
        sample = samplesDescription[sampleName]
        if "name" not in sample.keys():
            sample["name"] = sampleName

        # and make histograms
        main(args.mode, binning, sample, args.processor, outputDirectory, float(args.lumi))


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


