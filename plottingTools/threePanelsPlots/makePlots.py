import ROOT
import math
import os
import sys
import json
from array import array

sys.path.append("../common/")
import histogramUtils as hst


ROOT.gROOT.SetBatch(True)

# Output path for plots
OUTPUT_PATH = "./figures"

# Maximum number of bins for histogram
# Also maximum number of point of middle and bottom panels TGraphs
NBINS_MAX = 60

# Normalize histograms to 1
NORMALIZE_TO_1 = True


## Define variables, processes, colors

# Variables 
with open("varLabels.json", 'r') as f:
    varLabels = json.load(f)


# Processes to plot
pathToHistograms = "/eos/user/f/fleble/SVJ/data/histograms/"
processes = {
    # varying mMed
    #"tchannel_mMed-1000_mDark-20_rinv-0p3": [ pathToHistograms + "tchannel_mMed-1000_mDark-20_rinv-0p3.root" ],
    #"tchannel_mMed-3000_mDark-20_rinv-0p3": [ pathToHistograms + "tchannel_mMed-3000_mDark-20_rinv-0p3.root" ],
    #"tchannel_mMed-4000_mDark-20_rinv-0p3": [ pathToHistograms + "tchannel_mMed-4000_mDark-20_rinv-0p3.root" ],

    # varying rinv
    "tchannel_mMed-3000_mDark-20_rinv-0p1": [ pathToHistograms + "tchannel_mMed-3000_mDark-20_rinv-0p1.root" ],
    "tchannel_mMed-3000_mDark-20_rinv-0p3": [ pathToHistograms + "tchannel_mMed-3000_mDark-20_rinv-0p3.root" ],
    "tchannel_mMed-3000_mDark-20_rinv-0p7": [ pathToHistograms + "tchannel_mMed-3000_mDark-20_rinv-0p7.root" ],
    "QCD": [ 
         pathToHistograms + "QCD_Pt_170to300.root",
         pathToHistograms + "QCD_Pt_300to470.root",
         pathToHistograms + "QCD_Pt_470to600.root",
         pathToHistograms + "QCD_Pt_600to800.root",
         pathToHistograms + "QCD_Pt_800to1000.root",
         pathToHistograms + "QCD_Pt_1000to1400.root",
         pathToHistograms + "QCD_Pt_1400to1800.root",
         pathToHistograms + "QCD_Pt_1800to2400.root",
         pathToHistograms + "QCD_Pt_2400to3200.root",
         pathToHistograms + "QCD_Pt_3200toInf.root"
         ]
    }

# Label for the processes (legend)
processLabels = {
    "tchannel_mMed-3000_mDark-20_rinv-0p3": "mMed=3000 mDark=20 rinv=0.3",
    "tchannel_mMed-1000_mDark-20_rinv-0p3": "mMed=1000 mDark=20 rinv=0.3",
    "tchannel_mMed-6000_mDark-20_rinv-0p3": "mMed=6000 mDark=20 rinv=0.3",
    "tchannel_mMed-3000_mDark-20_rinv-0p1": "mMed=3000 mDark=20 rinv=0.1",
    "tchannel_mMed-3000_mDark-20_rinv-0p7": "mMed=3000 mDark=20 rinv=0.7",
    "QCD": "QCD"
}

# Processes to stack or not stack
stack = {
    "tchannel_mMed-3000_mDark-20_rinv-0p3": False,
    "tchannel_mMed-1000_mDark-20_rinv-0p3": False,
    "tchannel_mMed-6000_mDark-20_rinv-0p3": False,
    "tchannel_mMed-3000_mDark-20_rinv-0p1": False,
    "tchannel_mMed-3000_mDark-20_rinv-0p7": False,
    "QCD": True
}

# Define signal and background
# Background is a list of processes to add together (2 bkgs = 1 bkg curve representing the sum of the 2 bkgs)
# Signal is a list a signals that will be treated separately (2 signals = 2 graphs/histograms in plot)
processClass = {
    "signal": [
        # varying mMed
        #"tchannel_mMed-1000_mDark-20_rinv-0p3",
        #"tchannel_mMed-3000_mDark-20_rinv-0p3",
        #"tchannel_mMed-6000_mDark-20_rinv-0p3",

        # varying rinv
        "tchannel_mMed-3000_mDark-20_rinv-0p1",
        "tchannel_mMed-3000_mDark-20_rinv-0p3",
        "tchannel_mMed-3000_mDark-20_rinv-0p7",

        ],
    "background": ["QCD"]
}

# Color of the different processes
colors = {
    # varying mMed
    #"tchannel_mMed-1000_mDark-20_rinv-0p3": ROOT.kBlue+1,
    #"tchannel_mMed-3000_mDark-20_rinv-0p3": ROOT.kRed+1,
    #"tchannel_mMed-6000_mDark-20_rinv-0p3": ROOT.kGreen+2,

    # varying rinv
    "tchannel_mMed-3000_mDark-20_rinv-0p1": ROOT.kBlue+1,
    "tchannel_mMed-3000_mDark-20_rinv-0p3": ROOT.kRed+1,
    "tchannel_mMed-3000_mDark-20_rinv-0p7": ROOT.kGreen+2,

    # Background
    "QCD": ROOT.TColor.GetColor(155, 152, 204),
    }

colors["background"] = ROOT.TColor.GetColor(155, 152, 204)

colors["signal"] = []
for signal in processClass["signal"]:
    colors["signal"].append(colors[signal])


# Function to compute significance 

def calcSignificance(cumulativeSumSig, cumulativeSumBkg):
    '''
    Return significance given cumulative sum of signal and background histograms.
    Input:  list, list (or list-like objects)
    Output: list

    '''

    nbins = len(cumulativeSumSig) - 1
    
    # If cumulative sums start with zeros
    if cumulativeSumSig[0] + cumulativeSumBkg[0] == 0:
        ibin0 = 1
        while cumulativeSumSig[ibin0] + cumulativeSumBkg[ibin0] == 0:
            ibin0 += 1
        part = ibin0 * [0.]
        range_ = range(ibin0, nbins+1)
        significance = part + [ cumulativeSumSig[ibin] / math.sqrt(cumulativeSumSig[ibin] + cumulativeSumBkg[ibin]) for ibin in range_ ]

    # If cumulative sums ends with zeros
    elif cumulativeSumSig[nbins] + cumulativeSumBkg[nbins] == 0:
        ibin0 = nbins-1
        while cumulativeSumSig[ibin0] + cumulativeSumBkg[ibin0] == 0:
            ibin0 -= 1
        part = (nbins-ibin0) * [0.]
        range_ = range(0, ibin0+1)
        significance = [ cumulativeSumSig[ibin] / math.sqrt(cumulativeSumSig[ibin] + cumulativeSumBkg[ibin]) for ibin in range_ ] + part
    
    else:
        significance = [ cumulativeSumSig[ibin] / math.sqrt(cumulativeSumSig[ibin] + cumulativeSumBkg[ibin]) for ibin in range(nbins+1) ]

    return significance
    


## Main function making 3 panels plots

def makePlot(variable, rootFiles, outputPath):

    # Define plot style
    ROOT.gStyle.SetOptStat(0)

    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(600)
    ROOT.gStyle.SetCanvasDefW(600)
    ROOT.gStyle.SetCanvasDefX(0)
    ROOT.gStyle.SetCanvasDefY(0)

    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.15)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.02)

    ROOT.gStyle.SetHistLineColor(1)
    ROOT.gStyle.SetHistLineStyle(0)
    ROOT.gStyle.SetHistLineWidth(1)
    ROOT.gStyle.SetEndErrorSize(2)
    ROOT.gStyle.SetMarkerStyle(20)
    ROOT.gStyle.SetMarkerSize(0.9)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetTitleTextColor(1)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFontSize(0.05)

    ROOT.gStyle.SetTitleColor(1, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.07, "XYZ")
    ROOT.gStyle.SetTitleXOffset(1.00)
    ROOT.gStyle.SetTitleYOffset(0.90)

    ROOT.gStyle.SetLabelColor(1, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
    ROOT.gStyle.SetLabelSize(0.06, "XYZ")

    ROOT.gStyle.SetAxisColor(1, "XYZ")
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetTickLength(0.03, "XYZ")
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    ROOT.gStyle.SetPaperSize(20., 20.)
    ROOT.gStyle.SetHatchesLineWidth(5)
    ROOT.gStyle.SetHatchesSpacing(0.05)

    ROOT.TGaxis.SetExponentOffset(-0.08, 0.01, "Y")


    # Get histograms from files
    # Merge histograms from the different files corresponding to 1 process
    # Histograms are assumed to be correctly normalized to xs*lumi
    histograms = {}
    for process in processes.keys():
        firstPass = True
        for file_ in rootFiles[process]:
            if firstPass:
                histograms[process] = hst.getHistogram(file_, variable)
                firstPass = False
            else:
                histograms[process].Add(hst.getHistogram(file_, variable))

    # Sanity check: all histograms must have the same numbers of bins
    # Histograms should have the same binning
    process0 = next(iter(histograms))
    nbins0 = histograms[process0].GetXaxis().GetNbins()
    for process in histograms.keys():
        nbins = histograms[process].GetXaxis().GetNbins()
        if nbins != nbins0:
            print("ERROR: All histograms for the different processes must have \
                   the same number of bins. Not producing plots.")
            return
    
    # Normalize histogram to 1 if asked
    if NORMALIZE_TO_1:
        for process in histograms.keys():
            histograms[process].Scale(1./histograms[process].Integral())

    # If histograms have too many bins, rebin it
    if nbins0 > NBINS_MAX:
        rebinFound = False
        ngroup = 2
        while not rebinFound:
            if nbins0%ngroup == 0 and nbins0/ngroup < NBINS_MAX:
                rebinFound = True
            else:
                ngroup += 1
        print("Rebin histograms %s with %d groups. New numbers of bins = %d"\
               %(variable, ngroup, nbins/ngroup))
        for process in histograms.keys():
            histograms[process].Rebin(ngroup)
    

    # Make histograms for signal and background by adding togther histograms
    # of the relevant processes
    histogramSigBkg = {"signal": []}
    for process in histograms.keys():
        if process in processClass["background"]:
            if "background" in histogramSigBkg.keys():
                histogramSigBkg["background"].Add(histograms[process])
            else:
                histogramSigBkg["background"] = histograms[process].Clone()

    N_SIGNAL = len(processClass["signal"])
    for iSignal in range(N_SIGNAL):
        histogramSigBkg["signal"].append(None)
        for process in histograms.keys():
            if process in processClass["signal"][iSignal]:
                histogramSigBkg["signal"][iSignal] = histograms[process].Clone()


    # Stack histograms
    tstack = ROOT.THStack("", "")
    for process in histograms.keys():
        if stack[process]:
            hist = histograms[process]
            hist.SetLineWidth(0)
            hist.SetFillColor(colors[process])
            tstack.Add(hist)
            
    
    # Make X axis label
    name = variable
    if name in varLabels:
        XLABEL = varLabels[name]
    else:
        XLABEL = name


    # Make canvas
    canvas = ROOT.TCanvas("", "", 600, 700)
    canvas.Divide(1, 3, 0, 0, 0)   # Divide canvas without margins between subcanvases


    ###################
    ###  Top panel  ###
    ###################

    canvas.cd(1)
    ROOT.gPad.SetRightMargin(0.03)

    # Draw histograms
    tstack.Draw("HIST")
    for process in histograms.keys():
        if not stack[process]:
            hist = histograms[process]
            hist.SetLineWidth(3)
            hist.SetLineColor(colors[process])
            hist.Draw("HIST SAME")


    # Add y axis labels
    tstack.GetYaxis().SetTitle("N_{Events} ")
    ROOT.gPad.SetLogy()


    # Adapt maximum
    tstack.SetMaximum(max(tstack.GetMaximum(), max([histograms[process].GetMaximum() for process in histograms.keys() if not stack[process]])) * 100)
    if NORMALIZE_TO_1:
        tstack.SetMinimum(1e-6)
    else:
        tstack.SetMinimum(1.0)


    # Add legend
    legend = ROOT.TLegend(0.3, 0.82-0.1*N_SIGNAL//2, 0.93, 0.97)
    legend.SetNColumns(2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    for process in histograms.keys():
        processLabel = processLabels[process]
        if stack[process]:
            legend.AddEntry(histograms[process], processLabel, "f")
        else:
            legend.AddEntry(histograms[process], processLabel, "lep")
    legend.Draw()


    ######################
    ###  Middle panel  ###
    ######################

    canvas.cd(2)
    ROOT.gPad.SetRightMargin(0.03)

    # Get number and edges of bins
    xaxis = histogramSigBkg["background"].GetXaxis()
    nbins = xaxis.GetNbins()
    binEdges = [ xaxis.GetBinLowEdge(i) for i in range(1, nbins+2) ]

    # Calculate cumulative sum of the background histogram and background efficiencies
    # lower for < cut
    # upper for > cut
    cumulativeSumLowerBkg, cumulativeSumUpperBkg, histSumBkg = hst.calcCumulativeSums(histogramSigBkg["background"])
    bkgEffLower = [ 100 * cumulativeSumLowerBkg[i]/histSumBkg for i in range(nbins+1) ]
    bkgEffUpper = [ 100 * cumulativeSumUpperBkg[i]/histSumBkg for i in range(nbins+1) ]

    sigEffLower = []
    sigEffUpper = []

    # Loop over all signals
    maxSignificance = 0
    tgraphSigLower = []
    tgraphSigUpper = []
    for iSignal in range(len(histogramSigBkg["signal"])):
        # Compute cumulative histogram
        cumulativeSumLowerSig, cumulativeSumUpperSig, histSumSig = hst.calcCumulativeSums(histogramSigBkg["signal"][iSignal])

        # Compute efficiencies
        sigEffLower.append([ 100 * cumulativeSumLowerSig[i]/histSumSig for i in range(nbins+1) ])
        sigEffUpper.append([ 100 * cumulativeSumUpperSig[i]/histSumSig for i in range(nbins+1) ])
        
        # Calculate significance for values defined by bin limits
        significanceLower = calcSignificance(cumulativeSumLowerSig, cumulativeSumLowerBkg)
        significanceUpper = calcSignificance(cumulativeSumUpperSig, cumulativeSumUpperBkg)

        # Draw significance ...
        # ... for lower cut
        if iSignal==0:
            tgraphSigLower0 = ROOT.TGraph(len(significanceLower), array("d", binEdges), array("d", significanceLower))
            tgraphSigLower0.SetMarkerColor(colors["signal"][iSignal])
            tgraphSigLower0.SetMarkerStyle(20)
            tgraphSigLower0.GetXaxis().SetRangeUser(binEdges[0], binEdges[-1])
            tgraphSigLower0.Draw("AP")
            maxSig = max(significanceLower)
            if maxSignificance < maxSig: maxSignificance = maxSig
        else:
            tgraphSigLower.append(ROOT.TGraph(len(significanceLower), array("d", binEdges), array("d", significanceLower)))
            
            tgraphSigLower[iSignal-1].SetMarkerColor(colors["signal"][iSignal])
            tgraphSigLower[iSignal-1].SetMarkerStyle(20)
            tgraphSigLower[iSignal-1].SetMarkerSize(1-0.15*iSignal)
            tgraphSigLower[iSignal-1].Draw("P SAME")
            maxSig = max(significanceLower)
            if maxSignificance < maxSig: maxSignificance = maxSig

        # ... for upper cut
        tgraphSigUpper.append(ROOT.TGraph(len(significanceUpper), array("d", binEdges), array("d", significanceUpper)))
        tgraphSigUpper[iSignal].SetMarkerColor(colors["signal"][iSignal])
        tgraphSigUpper[iSignal].SetMarkerStyle(21)
        tgraphSigUpper[iSignal].SetMarkerSize(1-0.15*iSignal)
        tgraphSigUpper[iSignal].Draw("P SAME")
        maxSig = max(significanceUpper)
        if maxSignificance < maxSig: maxSignificance = maxSig


    # Y axis range and label
    SET_LOG_Y = False
    if SET_LOG_Y:
        ROOT.gPad.SetLogy()
        ymax = 5*maxSignificance
        ymin = 8*10**-6
    else:
        ymax = 1.2*maxSignificance
        ymin = - ymax * 0.03
    tgraphSigLower0.GetYaxis().SetRangeUser(ymin, ymax)
    tgraphSigLower0.GetYaxis().SetTitle("Significance ")


    # Add legend
    tgraphDummy1 = ROOT.TGraph(1, array("d", [0.]), array("d", [0.]))
    tgraphDummy1.SetMarkerStyle(24)
    tgraphDummy1.SetMarkerColor(ROOT.kBlack)
    tgraphDummy2 = ROOT.TGraph(1, array("d", [0.]), array("d", [0.]))
    tgraphDummy2.SetMarkerStyle(25)
    tgraphDummy2.SetMarkerColor(ROOT.kBlack)
    
    legend2 = ROOT.TLegend(0.74, 0.85, 0.94, 0.98)
    legend2.SetNColumns(2)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0)
    legend2.AddEntry(tgraphDummy1, "< cut", "p")
    legend2.AddEntry(tgraphDummy2, "#geq cut", "p")
    legend2.Draw()


    ######################
    ###  Bottom panel  ###
    ######################

    canvas.cd(3)
    ROOT.gPad.SetRightMargin(0.03)

    # Make legend
    legend3 = ROOT.TLegend(0.3, 0.82-0.1*N_SIGNAL//2, 0.93, 0.97)
    legend3.SetNColumns(2)
    legend3.SetBorderSize(0)
    legend3.SetFillStyle(0)

    # Draw signal efficiencies
    tgraphSigEff = []
    for iSignal in range(len(sigEffUpper)):

        if iSignal == 0:
            tgraphSigEff0 = ROOT.TGraph(len(sigEffUpper[iSignal]), array("d", binEdges), array("d", sigEffUpper[iSignal]))
            tgraphSigEff0.SetMarkerColor(colors["signal"][iSignal])
            tgraphSigEff0.SetMarkerStyle(21)
            tgraphSigEff0.Draw("AP")
            legend3.AddEntry(tgraphSigEff0, processLabels[processClass["signal"][iSignal]], "p")
        else:
            tgraphSigEff.append(ROOT.TGraph(len(sigEffUpper[iSignal]), array("d", binEdges), array("d", sigEffUpper[iSignal])))
            tgraphSigEff[iSignal-1].SetMarkerColor(colors["signal"][iSignal])
            tgraphSigEff[iSignal-1].SetMarkerStyle(21)
            tgraphSigEff[iSignal-1].SetMarkerSize(1-0.15*iSignal)
            tgraphSigEff[iSignal-1].Draw("P SAME")
            legend3.AddEntry(tgraphSigEff[iSignal-1], processLabels[processClass["signal"][iSignal]], "p")


    # Make x and y labels, ranges
    tgraphSigEff0.GetXaxis().SetTitle(XLABEL)
    tgraphSigEff0.GetYaxis().SetTitle("Efficiency #geq cut [%] ")
    tgraphSigEff0.GetXaxis().SetRangeUser(binEdges[0], binEdges[-1])
    tgraphSigEff0.GetYaxis().SetRangeUser(0, 130)

    # Draw background efficiency
    tgraphBkgEff = ROOT.TGraph(len(bkgEffUpper), array("d", binEdges), array("d", bkgEffUpper))
    tgraphBkgEff.SetMarkerColor(colors["background"])
    tgraphBkgEff.SetMarkerStyle(21)
    tgraphBkgEff.Draw("P SAME")
    legend3.AddEntry(tgraphBkgEff, "Background", "p")

    # Draw legend
    legend3.Draw()


    # Save the whole canvas
    canvas.SaveAs(outputPath + "/{}.pdf".format(variable))
    #canvas.SaveAs(outputPath + "/{}.png".format(variable))


if __name__ == "__main__":

    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)

    # Open all ROOT files beforehand
    rootFiles = {}
    for process in processes.keys():
        rootFiles[process] = []
        for fileName in processes[process]:
            rootFiles[process].append(ROOT.TFile(fileName, "READ"))

    # Make plots
    for variable in varLabels.keys():
        makePlot(variable, rootFiles, OUTPUT_PATH)

