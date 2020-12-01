import sys
import os
import ROOT
import DataFormats.FWLite as fwlite #import Events, Handle
import numpy as np

sys.path.append("../../pythonUtils/")
import utilities as ut
from plotInfo import getPlotInfo 


ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)


## Some global constants
AOD_FORMATS = [ "MiniAOD", "NanoAOD", "PFNanoAOD" ]
NBINS = 30


## Path the mini, nano, PFnano AOD file to compare
#fileNameMiniAOD = "root://cmseos.fnal.gov//store/user/keanet/tchannel/SVJP_08272020/MINIAOD/step4_MINIAOD_t-channel_mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000_part-1.root"
fileNameMiniAOD = "/eos/home-f/fleble/SVJ/data/test/signalSamples/102X/t-channel/MiniAOD/step4_MINIAOD_t-channel_mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000_part-1.root"
fileNamePFNanoAOD = "/eos/user/f/fleble/SVJ/data/test/signalSamples/102X/t-channel/PFnanoAOD/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000_part-1.root"
fileNameNanoAOD = ""


def buildVarInfo(plotInfo):
    '''
    Build dictionary with information about the variable to plot
    Create handle and label outside of loop

    '''
    # Create a pointer to the subdictionary plotInfo["varInfo"]
    varInfo = plotInfo["varInfo"]

    varList0 = plotInfo["varInfo"].keys()
    plotInfo["varList"] = []     # List to store actual variables for which to plot distribution
    for var in varList0:
        if (var[-6:] != "common" and var != "defaults"): plotInfo["varList"].append(var)
        else: continue

	varInfoTmp = {}

	# If common variables are present, copy content into the plotInfo dict
	# First xxx_common
	varc = var.split("_")[0]+"_common"
	if (varc in varList0):
	    for key in varInfo[varc].keys(): varInfoTmp[key] = varInfo[varc][key]
	# Then xxx_xxx_common
	varc = ut.list2str(var.split("_")[:2], "_")+"_common"
	if (varc in varList0):
	    for key in varInfo[varc].keys(): varInfoTmp[key] = varInfo[varc][key]

	# Finally get the input from the variable (overide common if need be)
	for key in varInfo[var].keys(): varInfoTmp[key] = varInfo[var][key]

	# Create array for values read from the mini and nano AOD files
	varInfoTmp["arrayMiniAOD"] = []
	varInfoTmp["arrayNanoAOD"] = []
	varInfoTmp["arrayPFNanoAOD"] = []
	
	# Define optional parameters
	for par in varInfo["defaults"].keys():
	    if par not in varInfoTmp.keys():
		varInfoTmp[par] = varInfo["defaults"][par]

	# Update dict
	del varInfo[var]
	varInfo[var] = varInfoTmp.copy()

	# Get handle, label for non composite variables
	if (var[0] != "!"):
	    varInfo[var]["handle"] = fwlite.Handle(varInfo[var]["handleName"])
	    varInfo[var]["label"] = (varInfo[var]["labelName"])

    return 0


def getMiniAODdata(plotInfo, fileName):

    # Create a pointer to the subdictionary plotInfo["varInfo"]
    varInfo = plotInfo["varInfo"]

    ## Open miniAOD file
    events = fwlite.Events(fileName)

    nevt=0
    ## Loop over events
    for event in events:
        compositeVarInfoAvail = True
        nevt+=1
        baseVar = {}
        for var in plotInfo["varList"]:
            if (var[0] != "!"):
                event.getByLabel(varInfo[var]["label"], varInfo[var]["handle"])
                objs = varInfo[var]["handle"].product()
            # Jet
            if var[:3] == "Jet" or var[:6] == "FatJet":
                if (varInfo[var]["leadingIdx"]<objs.size()):
                    varValue = eval("objs["+str(varInfo[var]["leadingIdx"])+"]"+varInfo[var]["accessor"])
                    if varInfo[var]["plot"]: varInfo[var]["arrayMiniAOD"].append(varValue)
                    baseVar[var] = varValue
                else:
                    compositeVarInfoAvail = False
            # MET
            elif var[:3] == "MET":
                varValue = eval("objs.front()"+varInfo[var]["accessor"])
                if varInfo[var]["plot"]: varInfo[var]["arrayMiniAOD"].append(varValue)
                baseVar[var] = varValue
            # Composite variable
            elif var[0] == "!":
                if compositeVarInfoAvail:
		    varValue = eval(varInfo[var]["expr"].replace("#","baseVar"))
		    if varInfo[var]["plot"]: varInfo[var]["arrayMiniAOD"].append(varValue)
		    baseVar[var] = varValue
            # Other simple variable
            else:
		varValue = eval("objs["+str(varInfo[var]["leadingIdx"])+"]"+varInfo[var]["accessor"])
		if varInfo[var]["plot"]: varInfo[var]["arrayMiniAOD"].append(varValue)
		baseVar[var] = varValue
                
    return nevt


def getNanoAODdata(plotInfo, fileName):

    # Create a pointer to the subdictionary plotInfo["varInfo"]
    varInfo = plotInfo["varInfo"]

    ## Open nanoAOD
    f = ROOT.TFile.Open(fileName)
    tree = f.Get("Events")
    nEntries = tree.GetEntries()

    ## Loop over events
    for ie in range(nEntries):
        compositeVarInfoAvail = True
        tree.GetEntry(ie)
        baseVar = {}
        for var in plotInfo["varList"]:
            # Composite variable
            if var[0] == "!":
                if compositeVarInfoAvail:
		    varValue = eval(varInfo[var]["expr"].replace("#","baseVar"))
		    if varInfo[var]["plot"]: varInfo[var]["arrayPFNanoAOD"].append(varValue)
		    baseVar[var] = varValue
            # Jet
            elif var[:3] == "Jet" or var[:6] == "FatJet":
		if varInfo[var]["leadingIdx"]<tree.GetLeaf("nJet").GetValue(0):
		    varValue = tree.GetLeaf(varInfo[var]["leaf"]).GetValue(varInfo[var]["leadingIdx"])
		    if varInfo[var]["plot"]: varInfo[var]["arrayPFNanoAOD"].append(varValue)
		    baseVar[var] = varValue
		else:
		    compositeVarInfoAvail = False
            # Other simple variables
            else:
		varValue = tree.GetLeaf(varInfo[var]["leaf"]).GetValue(varInfo[var]["leadingIdx"])
		if varInfo[var]["plot"]: varInfo[var]["arrayPFNanoAOD"].append(varValue)
		baseVar[var] = varValue

    f.Close()   # Close file
    return nEntries


def getBinning(varInfo, method="minMaxHistograms"):

    def getMinBin(varInfo):
        minArray = []
        for AODformat in AOD_FORMATS:
            if len(varInfo["array"+AODformat]) > 0:
                minArray.append(min(varInfo["array"+AODformat]))
        minBin = min(minArray)
        minBin = minBin - 0.1*abs(minBin)
        return(minBin)

    def getMaxBin(varInfo):
        maxArray = []
        for AODformat in AOD_FORMATS:
            if len(varInfo["array"+AODformat]) > 0:
                maxArray.append(max(varInfo["array"+AODformat]))
        maxBin = max(maxArray)
        maxBin = maxBin + 0.1*abs(maxBin)
        return(maxBin)

    if method == "minMaxHistograms":
        return(getMinBin(varInfo), getMaxBin(varInfo))

    elif method == "zeroMaxHistograms":
        return(0, getMaxBin(varInfo))

    elif method == "minMaxKeys":
        return(varInfo["min"], varInfo["max"])

    else:
        print("WARNING: Unknown binning method for var, using min max of histograms" %varInfo["name"])
        return getMinBin(varInfo), getMaxBin(varInfo)



def makeHistograms(plotInfo):

    # Create a pointer to the subdictionary plotInfo["varInfo"]
    varInfo = plotInfo["varInfo"]

    plotInfo["varPlotList"] = []
    for var in plotInfo["varList"]:
        if varInfo[var]["plot"]:
            plotInfo["varPlotList"].append(var)
            minBin, maxBin = getBinning(varInfo[var], plotInfo["plotInfo"]["minMaxBinning"])

            for AODformat in AOD_FORMATS:
                if len(varInfo[var]["array"+AODformat]) > 0:
                    hTitle = varInfo[var]["name"]+AODformat
                    varInfo[var]["hist"+AODformat] = ROOT.TH1D(hTitle, hTitle, NBINS, minBin, maxBin)
                    varInfo[var]["hist"+AODformat].SetTitle(varInfo[var]["name"])
                    for x in varInfo[var]["array"+AODformat]: varInfo[var]["hist"+AODformat].Fill(x)
                else:
                    varInfo[var]["hist"+AODformat] = None


def makePlots(plotInfo, AOD1, AOD2):

    # Create a pointer to the subdictionary plotInfo["varInfo"]
    varInfo = plotInfo["varInfo"]

    for var in plotInfo["varPlotList"]:

        canvas = ROOT.TCanvas(varInfo[var]["name"], varInfo[var]["name"])
        if "logscale" in varInfo[var].keys() and varInfo[var]["logscale"]: canvas.SetLogy()

        varInfo[var]["hist"+AOD1].SetLineColor(2)
        varInfo[var]["hist"+AOD2].SetLineColor(3)
        varInfo[var]["hist"+AOD1].Draw()
        varInfo[var]["hist"+AOD2].Draw("SAMEHIST")

        # Update and draw legend
        legend = ROOT.TLegend(0.7,0.7,0.9,0.9)
        legend.SetFillStyle(0)
        legend.AddEntry(varInfo[var]["hist"+AOD1], AOD1, "lep")
        legend.AddEntry(varInfo[var]["hist"+AOD2], AOD2, "lep")
        legend.Draw()

        dirpath="figures/"+AOD1+"vs"+AOD2+"/"
        if not os.path.exists(dirpath): os.makedirs(dirpath)
        canvas.SaveAs(dirpath+var.replace("!","")+".pdf")
        del canvas



if __name__ == "__main__":

    print("\nReading plot info...")
    plotInfo = getPlotInfo()
    buildVarInfo(plotInfo)


    print("\nReading files...")
    x = plotInfo["comparisons"]
    if x["MiniAODvsPFNanoAOD"] or x["MiniAODvsNanoAOD"]:
        print("Reading %s" %fileNameMiniAOD)
        nEvt = getMiniAODdata(plotInfo, fileNameMiniAOD)
        print("%d events read" %nEvt)

    if x["MiniAODvsPFNanoAOD"] or x["NanoAODvsPFNanoAOD"]:
        print("Reading %s" %fileNamePFNanoAOD)
        nEvt = getNanoAODdata(plotInfo, fileNamePFNanoAOD)
        print("%d events read" %nEvt)

    if x["MiniAODvsNanoAOD"] or x["NanoAODvsPFNanoAOD"]:
        print("Reading %s" %fileNameNanoAOD)
        nEvt = getNanoAODdata(plotInfo, fileNameNanoAOD)
        print("%d events read" %nEvt)


    print("\nMaking histograms...")
    makeHistograms(plotInfo)
  
    print("\nMaking plots...")
    for comparison in plotInfo["comparisons"].keys():
        if plotInfo["comparisons"][comparison]:
            [ AOD1, AOD2 ] = comparison.split("vs")
            makePlots(plotInfo, AOD1, AOD2)
