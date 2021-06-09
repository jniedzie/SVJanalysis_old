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


def make_x_label(var_labels, name):
    """Return xlabel name."""

    if name in var_labels:
        xlabel = var_labels[name]
    else:
        xlabel = name
    return xlabel


def calc_cumulative_sums(histogram_sig_bkg):
    """Calculate cumulative sum of the signal and background histograms.

    Args:
        histogram_sig_bkg (dict): Signals and background histograms:
            { "background": ROOT.TH1D, "signal": list[ROOT.TH1D] }
        
    Returns:
        dict: cumulative sums for <, >= cuts, for signals and background
            form of the dict:
            { "lowerBkg": list[float], "upperBkg": list[float],
              "lowerSig": list[list[float]], "upperSig": list[list[float]],
              "sumBkg": float", sumSig": list[float] }
    """

    cumulatives = {}

    # Calculate cumulative sum of the background histogram
    # lower for < cut
    # upper for > cut
    cumulatives["lowerBkg"], cumulatives["upperBkg"], cumulatives["sumBkg"] = hst.calc_cumulative_sums(histogram_sig_bkg["background"])

    # Loop over all signals
    nsignals = len(histogram_sig_bkg["signal"])
    cumulatives["lowerSig"] = nsignals * [0.]
    cumulatives["upperSig"] = nsignals * [0.]
    cumulatives["sumSig"]   = nsignals * [0.]
    for isignal in range(nsignals):
        # Compute cumulative histogram
        cumulatives["lowerSig"][isignal], cumulatives["upperSig"][isignal], cumulatives["sumSig"][isignal] = hst.calc_cumulative_sums(histogram_sig_bkg["signal"][isignal])

    return cumulatives


def calc_efficiencies(cumulatives):
    """Calculate efficiencies from cumulative sum.

    Args:
        cumulatives (dict): Cumulative sums for <, >= cuts, for signals and background
            form of the dict:
            { "lowerBkg": list[float], "upperBkg": list[float],
              "lowerSig": list[list[float]], "upperSig": list[list[float]],
              "sumBkg": float, "sumSig": list[float] }

    Returns:
        dict: Efficiencies for <, >= cuts, for signals and background
            form of the dict:
            { "lowerBkg": list[float], "upperBkg": list[float],
              "lowerSig": list[list[float]],  "upperBkg": list[list[float]] }
    """

    efficiencies = {}
    efficiencies["lowerSig"] = []
    efficiencies["upperSig"] = []
    for isignal in range(len(cumulatives["lowerSig"])):
        efficiencies["lowerSig"].append(hst.calc_efficiencies(cumulatives["lowerSig"][isignal]))
        efficiencies["upperSig"].append(hst.calc_efficiencies(cumulatives["upperSig"][isignal]))
 
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


def make_histogram_sig_bkg(histograms, process_class):
    """Make histograms for signal and background.
 
    Make histograms for the different signal points and and for the sum of all
    background by adding together histograms of the relevant processes.

    Args:
        histograms (dict[str, ROOT.TH1D]])
            Keys are process names
            Values are histograms
        process_class (dict(str, list[str]))

    Returns:
        dict: { "background": ROOT.TH1D, "signal": list[ROOT.TH1D] }
    """

    histogram_sig_bkg = {"signal": []}
    for process in histograms.keys():
        if process in process_class["background"]:
            if "background" in histogram_sig_bkg.keys():
                histogram_sig_bkg["background"].Add(histograms[process])
            else:
                histogram_sig_bkg["background"] = histograms[process].Clone()

    for isignal in range(N_SIGNAL):
        histogram_sig_bkg["signal"].append(None)
        for process in histograms.keys():
            if process in process_class["signal"][isignal]:
                histogram_sig_bkg["signal"][isignal] = histograms[process].Clone()

    return histogram_sig_bkg


def make_thstack(histograms, stack):
    """Make a ROOT.THStack histogram from histograms.

    Args:
        histograms (dict[str, ROOT.TH1D]]):
            Keys are process names
            Values are histograms
        stack (dict[str, bool]):
            Keys are process names
            Values are booleans for whether the histogram for this process
            is to be stacked or not

    Returns:
        ROOT.THStack
    """

    hstack = ROOT.THStack("", "")
    for process in histograms.keys():
        if stack[process]:
            hist = histograms[process]
            hist.SetLineWidth(0)
            hist.SetFillColor(colors[process])
            hstack.Add(hist)

    return hstack
 

def read_histograms(variable, root_files, processes, nbins_max=50, normalize_to_1=True, plot_overflow_bin=True):
    """
    Get histograms from input ROOT files.
    Merge histograms from the different files corresponding to 1 process.
    Histograms are assumed to be correctly normalized to xs*lumi.

    Args:
        variable (str): Name of the TH1D to get in the NTuple
        root_files (dict[str, list[ROOT.TFile]]):
            Keys are process names
            Values are list of ROOT.TFile from which to read the histograms
        processes (dict[str, list[str]]): List of ROOT files for each process
            Keys are process names
            Values are list of str path to ROOT files from which to read the histograms
        nbins_max (int)
        normalize_to_1 (bool)
        plot_overflow_bin (boolean)

    Returns:
        tuple: tuple containing:
            histograms (dict[str, ROOT.TH1D]])
            histogram_sig_bkg (dict): { "background": ROOT.TH1D, "signal": list[ROOT.TH1D] }
            hstack (ROOT.THStack)

    """

    histograms = {}
    for process in processes.keys():
        first_pass = True
        for file_ in root_files[process]:
            if first_pass:
                histograms[process] = hst.get_histogram(file_, variable)
                first_pass = False
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
    if normalize_to_1:
        for process in histograms.keys():
            hst.normalize(histograms[process], norm=1.)

    # If histograms have too many bins, rebin it
    if nbins0 > nbins_max:
        for iprocess, process in enumerate(histograms.keys()):
           verbose = True if iprocess == 0 else False
           hst.rebin_histogram(histograms[process], nbins_max, verbose=verbose)

    # Add overflow to the last bin
    if plot_overflow_bin:
        for process in histograms.keys():
            hst.add_overflow_bin_to_last_bin(histograms[process])

    # Make histograms for signal and background by adding together histograms
    # of the relevant processes
    histogram_sig_bkg = make_histogram_sig_bkg(histograms, process_class)

    # Stack histograms
    hstack = make_thstack(histograms, stack)

    return histograms, histogram_sig_bkg, hstack


def get_bin_info(histogram_sig_bkg):
    """Return number and edges of bins."""

    xaxis = histogram_sig_bkg["background"].GetXaxis()
    bin_info = {}
    bin_info["nbins"] = xaxis.GetNbins()
    bin_info["bin_edges"] = [ xaxis.GetBinLowEdge(i) for i in range(1, bin_info["nbins"]+2) ]
    
    return bin_info


def draw_histograms(histograms, hstack, xlabel=None):
    """Draw signal (empty) and background (filled) histograms on canvas."""

    hstack_clone = hstack.DrawClone("HIST")
    for process in histograms.keys():
        if not stack[process]:
            hist = histograms[process]
            hist.SetLineWidth(3)
            hist.SetLineColor(colors[process])
            hist.DrawClone("HIST SAME")

    # Adapt maximum
    hstack_clone.SetMaximum(max(hstack_clone.GetMaximum(), max([histograms[process].GetMaximum() for process in histograms.keys() if not stack[process]])) * 100)
    if NORMALIZE_TO_1:
        hstack_clone.SetMinimum(1e-6)
    else:
        hstack_clone.SetMinimum(1.0)

    # Axis labels
    hstack_clone.GetYaxis().SetTitle("N_{Events} ")
    if xlabel is not None:
        hstack_clone.GetXaxis().SetTitle(xlabel)

    ROOT.gPad.SetLogy()

    # Add legend
    legend = ROOT.TLegend(0.3, 0.82-0.1*N_SIGNAL//2, 0.93, 0.97)
    legend.SetNColumns(2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    for process in histograms.keys():
        process_label = process_labels[process]
        if stack[process]:
            legend.AddEntry(histograms[process], process_label, "f")
        else:
            legend.AddEntry(histograms[process], process_label, "lep")
    legend.DrawClone()

    return


def draw_significance(bin_info, cumulatives, histogram_sig_bkg, xlabel=None):
    """Draw significance as a function of cut value, for < and >= cuts."""

    ROOT.gPad.SetRightMargin(0.03)

    nbins = bin_info["nbins"]
    bin_edges = bin_info["bin_edges"]

    # Loop over all signals
    max_significance = 0
    tgraph_sig_lower = []
    tgraph_sig_upper = []
    for isignal in range(len(histogram_sig_bkg["signal"])):
        # Calculate significance for values defined by bin limits
        significance_lower = hst.calc_significance(cumulatives["lowerSig"][isignal], cumulatives["lowerBkg"])
        significance_upper = hst.calc_significance(cumulatives["upperSig"][isignal], cumulatives["upperBkg"])

        # Draw significance ...
        # ... for lower cut
        if isignal==0:
            tgraph_sig_lower_0 = ROOT.TGraph(len(significance_lower), array("d", bin_edges), array("d", significance_lower))
            tgraph_sig_lower_0.SetMarkerColor(colors["signal"][isignal])
            tgraph_sig_lower_0.SetMarkerStyle(20)
            tgraph_sig_lower_0.GetXaxis().SetRangeUser(bin_edges[0], bin_edges[-1])
            tgraph_sig_lower_0_clone = tgraph_sig_lower_0.DrawClone("AP")
            maxSig = max(significance_lower)
            if max_significance < maxSig: max_significance = maxSig
        else:
            tgraph_sig_lower.append(ROOT.TGraph(len(significance_lower), array("d", bin_edges), array("d", significance_lower)))
            
            tgraph_sig_lower[isignal-1].SetMarkerColor(colors["signal"][isignal])
            tgraph_sig_lower[isignal-1].SetMarkerStyle(20)
            tgraph_sig_lower[isignal-1].SetMarkerSize(1-0.15*isignal)
            tgraph_sig_lower[isignal-1].DrawClone("P SAME")
            maxSig = max(significance_lower)
            if max_significance < maxSig: max_significance = maxSig

        # ... for upper cut
        tgraph_sig_upper.append(ROOT.TGraph(len(significance_upper), array("d", bin_edges), array("d", significance_upper)))
        tgraph_sig_upper[isignal].SetMarkerColor(colors["signal"][isignal])
        tgraph_sig_upper[isignal].SetMarkerStyle(21)
        tgraph_sig_upper[isignal].SetMarkerSize(1-0.15*isignal)
        tgraph_sig_upper[isignal].DrawClone("P SAME")
        maxSig = max(significance_upper)
        if max_significance < maxSig: max_significance = maxSig


    # Y axis range
    SET_LOG_Y = False
    if SET_LOG_Y:
        ROOT.gPad.SetLogy()
        ymax = 5*max_significance
        ymin = 8*10**-6
    else:
        ymax = 1.2*max_significance
        ymin = - ymax * 0.03
    tgraph_sig_lower_0_clone.GetYaxis().SetRangeUser(ymin, ymax)

    # Axis labels
    tgraph_sig_lower_0_clone.GetYaxis().SetTitle("Significance ")
    if xlabel is not None:
        tgraph_sig_lower_0_clone.GetXaxis().SetTitle(xlabel)

    # Add legend
    tgraph_dummy_1 = ROOT.TGraph(1, array("d", [0.]), array("d", [0.]))
    tgraph_dummy_1.SetMarkerStyle(24)
    tgraph_dummy_1.SetMarkerColor(ROOT.kBlack)
    tgraph_dummy_2 = ROOT.TGraph(1, array("d", [0.]), array("d", [0.]))
    tgraph_dummy_2.SetMarkerStyle(25)
    tgraph_dummy_2.SetMarkerColor(ROOT.kBlack)
    
    legend = ROOT.TLegend(0.74, 0.85, 0.94, 0.98)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(tgraph_dummy_1, "< cut", "p")
    legend.AddEntry(tgraph_dummy_2, "#geq cut", "p")
    legend.DrawClone()

    return


def draw_efficiency(bin_info, efficiencies, xlabel):
    """Draw efficiency as a function of cut value, for >= cut."""

    nbins = bin_info["nbins"]
    bin_edges = bin_info["bin_edges"]

    # Make legend
    legend = ROOT.TLegend(0.3, 0.82-0.1*N_SIGNAL//2, 0.93, 0.97)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Draw signal efficiencies
    tgraph_sig_eff = []
    for isignal in range(len(efficiencies["upperSig"])):

        sig_eff_upper = efficiencies["upperSig"][isignal]
        if isignal == 0:
            tgraph_sig_eff_0 = ROOT.TGraph(len(sig_eff_upper), array("d", bin_edges), array("d", sig_eff_upper))
            tgraph_sig_eff_0.SetMarkerColor(colors["signal"][isignal])
            tgraph_sig_eff_0.SetMarkerStyle(21)
            tgraph_sig_eff_0_clone = tgraph_sig_eff_0.DrawClone("AP")
            legend.AddEntry(tgraph_sig_eff_0, process_labels[process_class["signal"][isignal]], "p")
        else:
            tgraph_sig_eff.append(ROOT.TGraph(len(sig_eff_upper), array("d", bin_edges), array("d", sig_eff_upper)))
            tgraph_sig_eff[isignal-1].SetMarkerColor(colors["signal"][isignal])
            tgraph_sig_eff[isignal-1].SetMarkerStyle(21)
            tgraph_sig_eff[isignal-1].SetMarkerSize(1-0.15*isignal)
            tgraph_sig_eff[isignal-1].DrawClone("P SAME")
            legend.AddEntry(tgraph_sig_eff[isignal-1], process_labels[process_class["signal"][isignal]], "p")


    # Make x and y labels, ranges
    if xlabel is not None:
        tgraph_sig_eff_0_clone.GetXaxis().SetTitle(xlabel)
    tgraph_sig_eff_0_clone.GetYaxis().SetTitle("Efficiency #geq cut [%] ")
    tgraph_sig_eff_0_clone.GetXaxis().SetRangeUser(bin_edges[0], bin_edges[-1])
    tgraph_sig_eff_0_clone.GetYaxis().SetRangeUser(0, 130)

    # Draw background efficiency
    bkg_eff_upper = efficiencies["upperBkg"]
    tgraph_bkg_eff = ROOT.TGraph(len(bkg_eff_upper), array("d", bin_edges), array("d", bkg_eff_upper))
    tgraph_bkg_eff.SetMarkerColor(colors["background"])
    tgraph_bkg_eff.SetMarkerStyle(21)
    tgraph_bkg_eff.DrawClone("P SAME")
    legend.AddEntry(tgraph_bkg_eff, "Background", "p")

    # Draw legend
    legend.DrawClone()

    return 


## Main function making 3 panels plots

def make_plot(variable, root_files, output_path, panels):
    """Make an N panel plot and returns a ROC AUC.

    Args:
        variable (str)
        root_files (dict[str, list[ROOT.TFile]]):
            Keys are process names
            Values are list of ROOT.TFile from which to read the histograms
        output_path (str)
        panels (list[str]): Defines, from top to bottom, the subplots to draw

    Returns:
        float: ROC AUC for the specified signal and all the backgrounds combined
    """

    histograms, histogram_sig_bkg, hstack = read_histograms(variable, root_files, processes, nbins_max=NBINS_MAX, normalize_to_1=NORMALIZE_TO_1, plot_overflow_bin=PLOT_OVERFLOW_BIN)
    bin_info = get_bin_info(histogram_sig_bkg)
    
    # Make X axis label
    xlabel = make_x_label(var_labels, variable)

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
            cumulatives = calc_cumulative_sums(histogram_sig_bkg)
        
        if panel == "histograms":
            draw_histograms(histograms, hstack, xlabel_)

        elif panel == "efficiency":
            if efficiencies is None:
                efficiencies = calc_efficiencies(cumulatives)
            draw_efficiency(bin_info, efficiencies, xlabel_)

        elif panel == "significance":
            draw_significance(bin_info, cumulatives, histogram_sig_bkg, xlabel_)
        
    # Calculate auc
    if efficiencies is None:
        if cumulatives is None:
            cumulatives = calc_cumulative_sums(histogram_sig_bkg)
        efficiencies = calc_efficiencies(cumulatives)
    auc1 = hst.calc_auc(efficiencies["lowerBkg"]      , efficiencies["lowerSig"][ISIGNAL_FOR_AUC]      , percentage=True)
    auc2 = hst.calc_auc(efficiencies["upperBkg"][::-1], efficiencies["upperSig"][ISIGNAL_FOR_AUC][::-1], percentage=True)
    if auc1 > 0.5: auc = auc1
    else: auc = auc2


    # Save the whole canvas
    canvas.SaveAs(output_path + "/{}.pdf".format(variable))
    #canvas.SaveAs(output_path + "/{}.png".format(variable))

    return auc


if __name__ == "__main__":

    # Import variables from config file
    from config import *

    define_plot_style()

    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)

    # Dataframe to store AUCs
    auc_dataframe = pd.DataFrame({"variable": [], "auc": []})

    # Open all ROOT files beforehand
    root_files = {}
    for process in processes.keys():
        root_files[process] = []
        for file_name in processes[process]:
            root_files[process].append(ROOT.TFile(file_name, "READ"))

    # Make plots
    for variable in var_labels.keys():
        auc = make_plot(variable, root_files, OUTPUT_PATH, panels)
        auc_dataframe = auc_dataframe.append({"variable": variable, "auc": auc}, ignore_index=True)

    # Sort dataframe by descending AUCs
    auc_dataframe.sort_values("auc", ascending=False, inplace=True)

    # Print aucs
    print("\nAUCs:")
    print(auc_dataframe)
    auc_dataframe.to_csv("auc.csv", index=False)

