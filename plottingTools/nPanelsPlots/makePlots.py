import ROOT
import math
import pandas as pd
import os
import sys
import json
from array import array

sys.path.append("../common/")
import histogramUtils as hst

ROOT.gROOT.SetBatch(True)


def make_x_label(varLabels, name):
    """Return xlabel name."""

    if name in varLabels:
        xlabel = varLabels[name]
    else:
        xlabel = name
    return xlabel


def calc_cumulative_sums(histogramSigBkg):
    """
    Calculate cumulative sum of the signal and background histograms.

    Parameters
    ----------
    histogramSigBkg: dict
        Form of the dict: { "background": ROOT.TH1D, "signal": list < ROOT.TH1D > }
    
    Returns
    -------
    cumulatives: dict
        form of the dict:
        { "lowerBkg": list < float >, "upperBkg": list < float >, "sumBkg": float,
          "lowerSig": list < list <float> >, "upperSig": list < list <float> >, "sumSig": list < float >
        }
    """

    cumulatives = {}

    # Calculate cumulative sum of the background histogram
    # lower for < cut
    # upper for > cut
    cumulatives["lowerBkg"], cumulatives["upperBkg"], cumulatives["sumBkg"] = hst.calc_cumulative_sums(histogramSigBkg["background"])

    # Loop over all signals
    nSignals = len(histogramSigBkg["signal"])
    cumulatives["lowerSig"] = nSignals * [0.]
    cumulatives["upperSig"] = nSignals * [0.]
    cumulatives["sumSig"]   = nSignals * [0.]
    for iSignal in range(nSignals):
        # Compute cumulative histogram
        cumulatives["lowerSig"][iSignal], cumulatives["upperSig"][iSignal], cumulatives["sumSig"][iSignal] = hst.calc_cumulative_sums(histogramSigBkg["signal"][iSignal])

    return cumulatives


def calc_efficiencies(cumulatives):
    """
    Calculate efficiencies from cumulative sum.

    Parameters
    ----------
    cumulatives: dict
        form of the dict:
        { "lowerBkg": list < float >, "upperBkg": list < float >, "sumBkg": float,
          "lowerSig": list < list <float> >, "upperSig": list < list <float> >, "sumSig": list < float >
        }

    Returns
    -------
    effiencies: dict
        form of the dict:
        { "lowerBkg": list < float >, "upperBkg": list < float >,
          "lowerSig": list < list < float > >,  "upperBkg": list < list < float > >
        }

    """

    efficiencies = {}
    efficiencies["lowerSig"] = []
    efficiencies["upperSig"] = []
    for iSignal in range(len(cumulatives["lowerSig"])):
        efficiencies["lowerSig"].append(hst.calc_efficiencies(cumulatives["lowerSig"][iSignal]))
        efficiencies["upperSig"].append(hst.calc_efficiencies(cumulatives["upperSig"][iSignal]))
 
    efficiencies["lowerBkg"] = hst.calc_efficiencies(cumulatives["lowerBkg"])
    efficiencies["upperBkg"] = hst.calc_efficiencies(cumulatives["upperBkg"])

    return efficiencies


def define_plot_style():
    """Define a CMS-style ROOT canvas."""

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

    return



def read_histograms(variable, rootFiles, processes, nBinsMax=50, normalizeTo1=True, plotOverflowBin=True):
    """
    Get histograms from input ROOT files.
    Merge histograms from the different files corresponding to 1 process.
    Histograms are assumed to be correctly normalized to xs*lumi.

    Parameters
    ----------
    variable: str
        Name of the TH1D to get in the NTuple
    rootFiles: dict < str, list < ROOT.TFile > >
        Keys are process names
        Values are list of ROOT.TFile from which to read the histograms
    processes: dict < str, list < str > >
    nBinsMax: int
    normalizeTo1: boolean
    plotOverflowBin: boolean

    Returns
    -------
    histograms: dict < str, ROOT.TH1D >
    histogramSigBkg: dict
        Form of the dict: { "background": ROOT.TH1D, "signal": list < ROOT.TH1D > }
    tstack: ROOT.THStack

    """

    histograms = {}
    for process in processes.keys():
        firstPass = True
        for file_ in rootFiles[process]:
            if firstPass:
                histograms[process] = hst.get_histogram(file_, variable)
                firstPass = False
            else:
                histograms[process].Add(hst.get_histogram(file_, variable))

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
    if normalizeTo1:
        for process in histograms.keys():
            integral = histograms[process].Integral()
            histograms[process].Scale(1./integral)

    # If histograms have too many bins, rebin it
    if nbins0 > nBinsMax:
        rebinFound = False
        ngroup = 2
        while not rebinFound:
            if nbins0%ngroup == 0 and nbins0/ngroup < nBinsMax:
                rebinFound = True
            else:
                ngroup += 1
        print("Rebin histograms %s with %d groups. New numbers of bins = %d"\
               %(variable, ngroup, nbins/ngroup))
        for process in histograms.keys():
            histograms[process].Rebin(ngroup)
    
        nbins0 = int(nbins0/ngroup)

    # Add overflow to the last bin
    if plotOverflowBin:
        for process in histograms.keys():
            overflow = histograms[process].GetBinContent(nbins+1)
            lastBinContent = histograms[process].GetBinContent(nbins0)
            histograms[process].SetBinContent(nbins0+1, 0.)
            histograms[process].SetBinContent(nbins0, lastBinContent+overflow)

    # Make histograms for signal and background by adding togther histograms
    # of the relevant processes
    histogramSigBkg = {"signal": []}
    for process in histograms.keys():
        if process in processClass["background"]:
            if "background" in histogramSigBkg.keys():
                histogramSigBkg["background"].Add(histograms[process])
            else:
                histogramSigBkg["background"] = histograms[process].Clone()

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
 
    return histograms, histogramSigBkg, tstack


def get_bin_info(histogramSigBkg):
    """Return number and edges of bins."""

    xaxis = histogramSigBkg["background"].GetXaxis()
    binInfo = {}
    binInfo["nbins"] = xaxis.GetNbins()
    binInfo["binEdges"] = [ xaxis.GetBinLowEdge(i) for i in range(1, binInfo["nbins"]+2) ]
    
    return binInfo


def draw_histograms(histograms, tstack, xlabel=None):
    """Draw signal (empty) and background (filled) histograms on canvas."""

    tstackClone = tstack.DrawClone("HIST")
    for process in histograms.keys():
        if not stack[process]:
            hist = histograms[process]
            hist.SetLineWidth(3)
            hist.SetLineColor(colors[process])
            hist.DrawClone("HIST SAME")

    # Adapt maximum
    tstackClone.SetMaximum(max(tstackClone.GetMaximum(), max([histograms[process].GetMaximum() for process in histograms.keys() if not stack[process]])) * 100)
    if NORMALIZE_TO_1:
        tstackClone.SetMinimum(1e-6)
    else:
        tstackClone.SetMinimum(1.0)

    # Axis labels
    tstackClone.GetYaxis().SetTitle("N_{Events} ")
    if xlabel is not None:
        tstackClone.GetXaxis().SetTitle(xlabel)

    ROOT.gPad.SetLogy()

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
    legend.DrawClone()

    return


def draw_significance(binInfo, cumulatives, histogramSigBkg, xlabel=None):
    """Draw significance as a function of cut value, for < and >= cuts."""

    ROOT.gPad.SetRightMargin(0.03)

    nbins = binInfo["nbins"]
    binEdges = binInfo["binEdges"]

    # Loop over all signals
    maxSignificance = 0
    tgraphSigLower = []
    tgraphSigUpper = []
    for iSignal in range(len(histogramSigBkg["signal"])):
        # Calculate significance for values defined by bin limits
        significanceLower = hst.calc_significance(cumulatives["lowerSig"][iSignal], cumulatives["lowerBkg"])
        significanceUpper = hst.calc_significance(cumulatives["upperSig"][iSignal], cumulatives["upperBkg"])

        # Draw significance ...
        # ... for lower cut
        if iSignal==0:
            tgraphSigLower0 = ROOT.TGraph(len(significanceLower), array("d", binEdges), array("d", significanceLower))
            tgraphSigLower0.SetMarkerColor(colors["signal"][iSignal])
            tgraphSigLower0.SetMarkerStyle(20)
            tgraphSigLower0.GetXaxis().SetRangeUser(binEdges[0], binEdges[-1])
            tgraphSigLower0Clone = tgraphSigLower0.DrawClone("AP")
            maxSig = max(significanceLower)
            if maxSignificance < maxSig: maxSignificance = maxSig
        else:
            tgraphSigLower.append(ROOT.TGraph(len(significanceLower), array("d", binEdges), array("d", significanceLower)))
            
            tgraphSigLower[iSignal-1].SetMarkerColor(colors["signal"][iSignal])
            tgraphSigLower[iSignal-1].SetMarkerStyle(20)
            tgraphSigLower[iSignal-1].SetMarkerSize(1-0.15*iSignal)
            tgraphSigLower[iSignal-1].DrawClone("P SAME")
            maxSig = max(significanceLower)
            if maxSignificance < maxSig: maxSignificance = maxSig

        # ... for upper cut
        tgraphSigUpper.append(ROOT.TGraph(len(significanceUpper), array("d", binEdges), array("d", significanceUpper)))
        tgraphSigUpper[iSignal].SetMarkerColor(colors["signal"][iSignal])
        tgraphSigUpper[iSignal].SetMarkerStyle(21)
        tgraphSigUpper[iSignal].SetMarkerSize(1-0.15*iSignal)
        tgraphSigUpper[iSignal].DrawClone("P SAME")
        maxSig = max(significanceUpper)
        if maxSignificance < maxSig: maxSignificance = maxSig


    # Y axis range
    SET_LOG_Y = False
    if SET_LOG_Y:
        ROOT.gPad.SetLogy()
        ymax = 5*maxSignificance
        ymin = 8*10**-6
    else:
        ymax = 1.2*maxSignificance
        ymin = - ymax * 0.03
    tgraphSigLower0Clone.GetYaxis().SetRangeUser(ymin, ymax)

    # Axis labels
    tgraphSigLower0Clone.GetYaxis().SetTitle("Significance ")
    if xlabel is not None:
        tgraphSigLower0Clone.GetXaxis().SetTitle(xlabel)

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
    legend2.DrawClone()

    return


def draw_efficiency(binInfo, efficiencies, xlabel):
    """Draw efficiency as a function of cut value, for >= cut."""

    nbins = binInfo["nbins"]
    binEdges = binInfo["binEdges"]

    # Make legend
    legend3 = ROOT.TLegend(0.3, 0.82-0.1*N_SIGNAL//2, 0.93, 0.97)
    legend3.SetNColumns(2)
    legend3.SetBorderSize(0)
    legend3.SetFillStyle(0)

    # Draw signal efficiencies
    tgraphSigEff = []
    for iSignal in range(len(efficiencies["upperSig"])):

        sigEffUpper = efficiencies["upperSig"][iSignal]
        if iSignal == 0:
            tgraphSigEff0 = ROOT.TGraph(len(sigEffUpper), array("d", binEdges), array("d", sigEffUpper))
            tgraphSigEff0.SetMarkerColor(colors["signal"][iSignal])
            tgraphSigEff0.SetMarkerStyle(21)
            tgraphSigEff0Clone = tgraphSigEff0.DrawClone("AP")
            legend3.AddEntry(tgraphSigEff0, processLabels[processClass["signal"][iSignal]], "p")
        else:
            tgraphSigEff.append(ROOT.TGraph(len(sigEffUpper), array("d", binEdges), array("d", sigEffUpper)))
            tgraphSigEff[iSignal-1].SetMarkerColor(colors["signal"][iSignal])
            tgraphSigEff[iSignal-1].SetMarkerStyle(21)
            tgraphSigEff[iSignal-1].SetMarkerSize(1-0.15*iSignal)
            tgraphSigEff[iSignal-1].DrawClone("P SAME")
            legend3.AddEntry(tgraphSigEff[iSignal-1], processLabels[processClass["signal"][iSignal]], "p")


    # Make x and y labels, ranges
    if xlabel is not None:
        tgraphSigEff0Clone.GetXaxis().SetTitle(xlabel)
    tgraphSigEff0Clone.GetYaxis().SetTitle("Efficiency #geq cut [%] ")
    tgraphSigEff0Clone.GetXaxis().SetRangeUser(binEdges[0], binEdges[-1])
    tgraphSigEff0Clone.GetYaxis().SetRangeUser(0, 130)

    # Draw background efficiency
    bkgEffUpper = efficiencies["upperBkg"]
    tgraphBkgEff = ROOT.TGraph(len(bkgEffUpper), array("d", binEdges), array("d", bkgEffUpper))
    tgraphBkgEff.SetMarkerColor(colors["background"])
    tgraphBkgEff.SetMarkerStyle(21)
    tgraphBkgEff.DrawClone("P SAME")
    legend3.AddEntry(tgraphBkgEff, "Background", "p")

    # Draw legend
    legend3.DrawClone()

    return 


## Main function making 3 panels plots

def makePlot(variable, rootFiles, outputPath, panels):
    """
    Make an N panel plot and returns a ROC AUC.

    Parameters
    ----------
    variable: str
    rootFiles: dict < str, list < ROOT.TFile > >
        Keys are process names
        Values are list of ROOT.TFile from which to read the histograms
    outputPath: str
    panels: list < str >
        Defines, from top to bottom, the subplots to draw

    Returns
    -------
    auc: float
        The ROC AUC for the specified signal and all the backgrounds combined

    """

    histograms, histogramSigBkg, tstack = read_histograms(variable, rootFiles, processes, nBinsMax=NBINS_MAX, normalizeTo1=NORMALIZE_TO_1, plotOverflowBin=PLOT_OVERFLOW_BIN)
    binInfo = get_bin_info(histogramSigBkg)
    
    # Make X axis label
    xlabel = make_x_label(varLabels, variable)

    # Make canvas
    canvas = ROOT.TCanvas("", "", 600, 700)
    canvas.Divide(1, len(panels), 0, 0, 0)   # Divide canvas without margins between subcanvases


    cumulatives = None
    efficiencies = None

    xlabel_ = None
    for ipanel, panel in enumerate(panels):
        canvas.cd(ipanel+1)
        ROOT.gPad.SetRightMargin(0.03)
        if ipanel == len(panels)-1:
            xlabel_ = xlabel
  
        if (panel == "efficiency" or panel == "significance") and cumulatives is None:
            cumulatives = calc_cumulative_sums(histogramSigBkg)
        
        if panel == "histograms":
            draw_histograms(histograms, tstack, xlabel_)

        elif panel == "efficiency":
            if efficiencies is None:
                efficiencies = calc_efficiencies(cumulatives)
            draw_efficiency(binInfo, efficiencies, xlabel_)

        elif panel == "significance":
            draw_significance(binInfo, cumulatives, histogramSigBkg, xlabel_)
        
    # Calculate auc
    if efficiencies is None:
        if cumulatives is None:
            cumulatives = calc_cumulative_sums(histogramSigBkg)
        efficiencies = calc_efficiencies(cumulatives)
    auc1 = hst.calc_auc(efficiencies["lowerBkg"]      , efficiencies["lowerSig"][iSignalForAuc]      , percentage=True)
    auc2 = hst.calc_auc(efficiencies["upperBkg"][::-1], efficiencies["upperSig"][iSignalForAuc][::-1], percentage=True)
    if auc1 > 0.5: auc = auc1
    else: auc = auc2


    # Save the whole canvas
    canvas.SaveAs(outputPath + "/{}.pdf".format(variable))
    #canvas.SaveAs(outputPath + "/{}.png".format(variable))

    return auc


if __name__ == "__main__":

    # Import variables from config file
    from config import *

    define_plot_style()

    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)

    # Dataframe to store AUCs
    aucDataFrame = pd.DataFrame({"variable": [], "auc": []})

    # Open all ROOT files beforehand
    rootFiles = {}
    for process in processes.keys():
        rootFiles[process] = []
        for fileName in processes[process]:
            rootFiles[process].append(ROOT.TFile(fileName, "READ"))

    # Make plots
    for variable in varLabels.keys():
        auc = makePlot(variable, rootFiles, OUTPUT_PATH, panels)
        aucDataFrame = aucDataFrame.append({"variable": variable, "auc": auc}, ignore_index=True)

    # Sort dataframe by descending AUCs
    aucDataFrame.sort_values("auc", ascending=False, inplace=True)

    # Print aucs
    print("\nAUCs:")
    print(aucDataFrame)
    aucDataFrame.to_csv("auc.csv", index=False)

