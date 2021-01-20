import sys
import json
import argparse
import time
import numpy as np

sys.path.append("../pythonUtils/")
import utilities as utl

import ROOT
ROOT.gROOT.SetBatch(True)


## Custom C++ function to be used in Define
ROOT.gInterpreter.Declare("""
Float_t DeltaPhiMinN(int NjetsMax, ROOT::VecOps::RVec<Float_t>& phi, Float_t& phi2) {
    ROOT::VecOps::RVec<Float_t> dPhi;
    int size = (int)phi.size();
    auto idxMax = std::min(size, NjetsMax);

    for (auto idx=0; idx<idxMax; idx++) {
        dPhi.push_back(ROOT::VecOps::DeltaPhi(phi[idx], phi2));
    }

    return(ROOT::VecOps::Min(dPhi));
}
""")



## Useful functions

# Book a histogram for a specific variable
def bookHistogram(df, variable, range_, weight):
    return df.Histo1D(ROOT.ROOT.RDF.TH1DModel(variable, variable, range_[0], range_[1], range_[2]), variable, weight)

# Write a histogram with a given name to the output ROOT file
def writeHistogram(h, name):
    h.SetName(name)
    h.Write()


## Main function making histograms

def main(MODE, variables, binning, sample, outputDirectory, LUMI):

    # Set up multi-threading capability of ROOT
    ROOT.ROOT.EnableImplicitMT()

    # Get sample name and cross-section
    sampleName = sample["name"]
    XSection = sample["XSection"]

    # Create/update output file
    ROOTfileName = outputDirectory + sampleName + ".root"
    print("\nWill %s ROOT file %s." %(MODE.lower(), ROOTfileName))
    tfile = ROOT.TFile(ROOTfileName, MODE)

    # Initialise event counter
    nGenEvts = 0

    # Make batches of files with a total of less than X million events
    # Need to do that because RDataFrame efficient features seems to break down
    # for too many events at once
    batches = [[]]
    nEvts = 0
    for file_ in sample["fileset"]:

        # If file is on eos, add global redirector
        if file_.startswith("/store/mc/") or file_.startswith("/store/user/"):
            file_ = "root://cms-xrd-global.cern.ch/" + file_

        if nEvts < 2e6:
            batches[-1].append(file_)
            f = ROOT.TFile.Open(file_, "READ")
            t = f.Get("Events")
            nEvts += t.GetEntries()
        else:
            nEvts = 0
            batches.append([file_])

    print("%d batches of files made" %(len(batches)))

    # Object to store histograms from the different batches of file sets
    hists = []


    # Loop over all batches of files
    for fileset in batches:
        print("")

        # Object to store histograms from the different files in this batch
        histsBatch = []

        # Loop over all files
        nFiles = len(fileset)
        for iFile, fileName in enumerate(fileset):

            tstart = time.time()
            print("Reading file %d out of %d" %(iFile+1, nFiles))

            histsBatch.append({})

            # If file is on eos, add global redirector
            if fileName.startswith("/store/mc/") or fileName.startswith("/store/user/"):
                fileName = "root://cms-xrd-global.cern.ch/" + fileName

            # Lists to store pointers to different RDataFrames (different filters)
            # and to store variables that have been defined
            dfList = []
            definedVars = []


            ## Make RDataFrames with all requested variables
            for idf, dataframe in enumerate(dataframes):

                # Check ordering of dataframes read
                dfIdx = dataframe["idx"]
                if dfIdx != idf:
                    print("ERROR: Mismatch between dataframe index and the number of dataframes defined so far.")
                    break

                # Define dataframe
                # Either read ROOT file (1st dataframe)
                if dfIdx == 0:
                    dfList.append(ROOT.ROOT.RDataFrame("Events", fileName))
                    # Increment the number of events processed
                    nGenEvts += dfList[0].Sum("genWeight").GetValue()  
                # Or build it from filtering an already defined dataframe
                else:
                    idxBase = dataframe["idxBase"]
                    f1 = dataframe["filter"][0]
                    f2 = dataframe["filter"][1]
                    dfList.append(dfList[idxBase].Filter(f1 , f2))

                # Define variables in this dataframe
                definedVars.append([])
                for define in dataframe["defines"]:
                    variable = define[0]
                    definition = define[1]
                    dfList[dfIdx] = dfList[dfIdx].Define(variable, definition)
                    definedVars[dfIdx].append(variable)

           
            ## Book histograms
            for variable in variables:

                # Check if binning is defined for this variable
                if variable not in binning["noregex"]:
                    regexes = list(binning["regex"].keys())
                    indices = utl.inregex(variable, regexes)
                    if len(indices) == 0:
                        print("Binning of %s is not defined. Skipping." %variable)
                        continue
                    elif len(indices) > 1:
                        print("%s matches several regexes. Binning cannot be defined. Skipping." %variable)
                        continue
                    else: 
                        binning_ = binning["regex"][regexes[indices[0]]]
                else:
                    binning_ = binning["noregex"][variable]

                # Find from which RDataFrame the variable has been defined
                for idx in range(len(definedVars)):
                    if variable in definedVars[idx]:
                        break
                    else:
                        # If variable is not defined, then it's taken from the uncut dataframe
                        if idx == len(definedVars)-1:
                            idx = 0

                # Check if variable present in dataframe
                if not variable in dfList[idx].GetColumnNames():
                    print("WARNING: Variable %s is not in dataframe (index %d). Will not be saved to output ROOT file." %(variable, idx))
                    continue
                
                # Book histogram
                # Histograms should NOT be added together yet as RDataFrame would not proceed in one loop in an efficient way
                # Should be done at the very end
                weight = "genWeight"
                histsBatch[iFile][variable] = bookHistogram(dfList[idx], variable, binning_, weight)


            ## For sanity check, print defined variables not asked to be saved in histogram
            definedVarsFlat = []
            for idx in range(len(definedVars)):
                definedVarsFlat = definedVarsFlat + definedVars[idx]
            unsavedVars = [ variable for variable in definedVarsFlat if variable not in variables ]
            if len(unsavedVars) > 0:
                print("INFO: The following variable were defined but not instructed to be saved in ROOT file:")
                for x in unsavedVars: print("\t%s" %x)
                print("")

            elapsed = time.time() - tstart
            print("Elapsed time: %d s" %elapsed)


        ## Adding histograms of the batch together
        print("Adding together histograms of batch %d..." %(len(hists)+1))
        tstrat = time.time()
        hists.append({})
        for variable in histsBatch[0].keys():
            hists[-1][variable] = histsBatch[0][variable]
            for iFile in range(1, len(histsBatch)):
                hists[-1][variable].Add(histsBatch[iFile][variable].GetValue())
        elapsed = time.time() - tstart
        print("Elapsed time: %d s" %elapsed)


    ## Normalize histograms and write to output ROOT file
    print("\nNormalizing histograms and writing to output ROOT file...")
    tstart = time.time()
    tfile.cd()
    for variable in hists[0].keys():
        hist = hists[0][variable]
        for iFile in range(1, len(hists)):
            hist.Add(hists[iFile][variable].GetValue())
        hist.Scale( XSection*LUMI/nGenEvts )
        writeHistogram(hist, "{}".format(variable))
        print("%s histogram saved" %variable)
    elapsed = time.time() - tstart
    print("Elapsed time: %d s" %elapsed)

    tfile.Close()



if __name__ == "__main__":

    tstart = time.time()

    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--mode",
        choices=["RECREATE", "UPDATE"],
        help="Mode in which to open the output ROOT file"
        )
    parser.add_argument(
        "-v", "--variables",
        help="json file listing variables to save"
        )
    parser.add_argument(
        "-d", "--dataframes",
        help="json file describing filters and defines of the RDataFrames needed to comoute the required variables"
        )
    parser.add_argument(
        "-b", "--binning",
        help="json file describing binning of the histograms"
        )
    parser.add_argument(
        "-sd", "--samplesDescription",
        help="json file describing samples location"
        )
    parser.add_argument(
        "-s", "--samples",
        help="Comma separated list samples to pick up among the samples described in the sample file"
        )
    parser.add_argument(
        "-o", "--outputDirectory",
        nargs="?",
        const="./",
        help="Path to the directory where to recreate/update ROOT file"
        )
    parser.add_argument(
        "-l", "--lumi",
        nargs="?",
        default=59725.0,   # 21071.0+38654.0
        help="Total luminosity for normalization of the histograms"
        )


    args = parser.parse_args()
    

    # All samples description
    with open(args.samplesDescription, 'r') as f:
        samplesDescription = json.load(f)

    # Variables to histogram
    with open(args.variables, 'r') as f:
        variables = json.load(f)["variables"]

    # Define the binning of the different variables to histogram
    with open(args.binning, 'r') as f:
        binning = json.load(f)["binning"]

    # Read filters and defines instructions to build up sequentially RDataFrames
    with open(args.dataframes, 'r') as f:
        dataframes = json.load(f)["dataframes"]

    # Get samples for which to make histograms
    if not args.samples:
        samplesNames = list(samplesDescription.keys())
    else:
        samplesNames = args.samples.split(",")

    # Loop over all samples
    for sampleName in samplesNames:
        sample = samplesDescription[sampleName]
        if "name" not in sample.keys():
            sample["name"] = sampleName
        # and make histograms
        main(args.mode, variables, binning, sample, args.outputDirectory, args.lumi)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)
