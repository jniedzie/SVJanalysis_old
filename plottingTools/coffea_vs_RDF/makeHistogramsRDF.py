
from collections import OrderedDict

import ROOT
ROOT.gROOT.SetBatch(True)

import time



LUMI = 21071.0+38654.0
XSection = 0.01891


# Used to define root files to read
samplesPath = "/eos/home-f/fleble/SVJ/data/production/102X/tchannel/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR/"
filesNumber = [ i for i in range(1,101) ]


# Define the binning of the different variables to histogram
ranges = {
#    "nGoodFatJet"          : (20 , 0 , 20  ),
    "GoodFatJet_pt"        : (200, 0 , 2000),
}

# Define the objects and variables to build up the variables to histogram
objDefinitions = OrderedDict()
o = objDefinitions
o["GoodFatJet"]       = "abs(FatJet_eta) < 2.4 && FatJet_pt > 200"
#o["nGoodFatJet"]      = "Sum(GoodFatJet)"
o["GoodFatJet_pt"]    = "FatJet_pt[GoodFatJet]"



# Book a histogram for a specific variable
def bookHistogram(df, variable, range_, weight):
    return df.Histo1D(ROOT.ROOT.RDF.TH1DModel(variable, variable, range_[0], range_[1], range_[2]), variable, weight)


# Write a histogram with a given name to the output ROOT file
def writeHistogram(h, name):
    h.SetName(name)
    h.Write()



def main():

    # Set up multi-threading capability of ROOT
    ROOT.ROOT.EnableImplicitMT(8)

    # Variables to histogram
    variables = ranges.keys()

    # Create output file
    tfile = ROOT.TFile("./histograms.root", "RECREATE")

    # Initialise event counter
    nGenEvts = 0

    # Dict to store histograms
    hists = {}


    # Loop through datasets and produce histograms of variables
    firstPass = True
    for fn in filesNumber:

        fileName = samplesPath + str(fn) + ".root"

        # Load dataset
        df = ROOT.ROOT.RDataFrame("Events", fileName)
        dfRun = ROOT.ROOT.RDataFrame("Runs", fileName)

        for obj, definition in objDefinitions.items():
            df = df.Define(obj, definition)

        # Increment the number of event processed
        nGenEvts += dfRun.Mean("genEventSumw_").GetValue()
        # nGenEvts += df.Sum("genWeight").GetValue()

        # Book histograms
        for variable in variables:
            weight = "genWeight"
            if firstPass:
                hists[variable] = bookHistogram(df, variable, ranges[variable], weight)
                firstPass = False
            else:
                h = bookHistogram(df, variable, ranges[variable], weight)
                hists[variable].Add(h.GetValue())


    # Normalise histograms
    for variable in variables:
        hists[variable].Scale( XSection*LUMI/nGenEvts )

    # Write histograms to output file
    for variable in variables:
        tfile.cd()
        writeHistogram(hists[variable], "{}".format(variable))

    tfile.Close()


if __name__ == "__main__":
    tstart = time.time()
    main()
    elapsed = time.time() - tstart
    print("Elapsed time: %d s" %elapsed)
