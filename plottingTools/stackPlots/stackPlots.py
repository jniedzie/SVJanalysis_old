import ROOT
import os

ROOT.gROOT.SetBatch(True)

# Output path for plots
OUTPUT_PATH = "./figures"


## Define variables, processes, colors

# Label for each variable (x-axis)
# This also defines variables that will be plotted
varLabels = {
    "nGoodFatJet": "Number of Fat Jets",
    "nGoodJet": "Number of Fat Jets",
    "J1_pt": "Fat Jet 1 p_{T} [GeV]",
    "GoodFatJet_pt": "Fat Jet p_{T} [GeV]",
     }

# Processes to plot
pathToHistograms = "/eos/user/f/fleble/SVJ/data/histograms/"
processes = {
    "tchannel": [ pathToHistograms + "tchannel_mMed-3000_mDark-20_rinv-0p3.root" ],
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
    "tchannel": "tchannel",
    "QCD": "QCD"
}

# Color of the different processes
colors = {
    "tchannel": ROOT.TColor.GetColor("#BF2229"),
    "QCD": ROOT.TColor.GetColor(155, 152, 204)
    }

# Processes to stack or not stack
stack = {
    "tchannel": False,
    "QCD": True
}


## Useful functions

# Retrieve a histogram from the input file based on the process and the variable
# name
def getHistogram(tfile, variable):
    h = tfile.Get(variable)
    if not h:
        raise Exception("Failed to load histogram {}.".format(variable))
    return h


## Main function making plots

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
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.05)

    ROOT.gStyle.SetHistLineColor(1)
    ROOT.gStyle.SetHistLineStyle(0)
    ROOT.gStyle.SetHistLineWidth(1)
    ROOT.gStyle.SetEndErrorSize(2)
    ROOT.gStyle.SetMarkerStyle(20)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetTitleTextColor(1)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFontSize(0.05)

    ROOT.gStyle.SetTitleColor(1, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleXOffset(1.00)
    ROOT.gStyle.SetTitleYOffset(1.60)

    ROOT.gStyle.SetLabelColor(1, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")

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
    # Histograms are assumed to be correctly normalized
    histograms = {}
    for process in processes.keys():
        firstPass = True
        for file_ in rootFiles[process]:
            if firstPass:
                histograms[process] = getHistogram(file_, variable)
                firstPass = False
            else:
                histograms[process].Add(getHistogram(file_, variable))

    # Stack histograms
    tstack = ROOT.THStack("", "")
    for process in histograms.keys():
        if stack[process]:
            hist = histograms[process]
            hist.SetLineWidth(0)
            hist.SetFillColor(colors[process])
            tstack.Add(hist)
            
    # Make canvas
    c = ROOT.TCanvas("", "", 600, 600)

    # Draw histograms
    tstack.Draw("HIST")
    for process in histograms.keys():
        if not stack[process]:
            hist = histograms[process]
            hist.SetLineWidth(3)
            hist.SetLineColor(colors[process])
            hist.Draw("HIST SAME")


    # Add axis labels
    name = variable
    if name in varLabels:
        title = varLabels[name]
    else:
        title = name
    tstack.GetXaxis().SetTitle(title)
    tstack.GetYaxis().SetTitle("N_{Events}")
    ROOT.gPad.SetLogy()


    # Adapt maximum
    tstack.SetMaximum(max(tstack.GetMaximum(), max([histograms[process].GetMaximum() for process in histograms.keys() if not stack[process]])) * 1.4)
    tstack.SetMinimum(1.0)


    # Add legend
    legend = ROOT.TLegend(0.4, 0.73, 0.90, 0.88)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    for process in histograms.keys():
        processLabel = processLabels[process]
        if stack[process]:
            legend.AddEntry(histograms[process], processLabel, "f")
        else:
            legend.AddEntry(histograms[process], processLabel, "lep")
    legend.Draw()

    # Save
    c.SaveAs(outputPath + "/{}.pdf".format(variable))
    #c.SaveAs(outputPath + "/{}.png".format(variable))


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

