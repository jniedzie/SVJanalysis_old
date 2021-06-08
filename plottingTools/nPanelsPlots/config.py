import json
import ROOT


# Output path for plots
OUTPUT_PATH = "./figures"

# Maximum number of bins for histogram
# Also maximum number of point of middle and bottom panels TGraphs
NBINS_MAX = 60

# Normalize histograms to 1
NORMALIZE_TO_1 = True

# PLOT overflow bin
PLOT_OVERFLOW_BIN = True

## Define variables, processes, colors

# Variables 
with open("varLabels_quick.json", 'r') as f:
    varLabels = json.load(f)
    # Quick hack
    varLabels = { x: x for x in sorted(varLabels) }


# Processes to plot
pathToHistogramsTChannel = "/work/fleble/t_channel_histograms/skim1/"
pathToHistogramsQCD = "/work/fleble/QCD_histograms/skim1/"
processes = {
    "tchannel_mMed-1000_mDark-20_rinv-0p3_alpha-peak": [ pathToHistogramsTChannel + "tchannel_mMed-1000_mDark-20_rinv-0.3_alpha-peak.root" ],
    "tchannel_mMed-3000_mDark-20_rinv-0p3_alpha-peak": [ pathToHistogramsTChannel + "tchannel_mMed-3000_mDark-20_rinv-0.3_alpha-peak.root" ],
    "tchannel_mMed-4000_mDark-20_rinv-0p3_alpha-peak": [ pathToHistogramsTChannel + "tchannel_mMed-4000_mDark-20_rinv-0.3_alpha-peak.root" ],

    "QCD": [ 
         pathToHistogramsQCD + "QCD_170_300.root",
         pathToHistogramsQCD + "QCD_300_470.root",
         pathToHistogramsQCD + "QCD_470_600.root",
         pathToHistogramsQCD + "QCD_600_800.root",
         pathToHistogramsQCD + "QCD_800_1000.root",
         pathToHistogramsQCD + "QCD_1000_1400.root",
         pathToHistogramsQCD + "QCD_1400_1800.root",
         pathToHistogramsQCD + "QCD_1800_2400.root",
         pathToHistogramsQCD + "QCD_2400_3200.root",
         pathToHistogramsQCD + "QCD_3200_Inf.root"
         ]
}

# Label for the processes (legend)
processLabels = {
    "tchannel_mMed-1000_mDark-20_rinv-0p3_alpha-peak": "mMed=1000 mDark=20 rinv=0.3 alpha=peak",
    "tchannel_mMed-3000_mDark-20_rinv-0p3_alpha-peak": "mMed=3000 mDark=20 rinv=0.3 alpha=peak",
    "tchannel_mMed-4000_mDark-20_rinv-0p3_alpha-peak": "mMed=4000 mDark=20 rinv=0.3 alpha=peak",
    "QCD": "QCD",
}

# Processes to stack or not stack
stack = {
    "tchannel_mMed-1000_mDark-20_rinv-0p3_alpha-peak": False,
    "tchannel_mMed-3000_mDark-20_rinv-0p3_alpha-peak": False,
    "tchannel_mMed-4000_mDark-20_rinv-0p3_alpha-peak": False,
    "QCD": True,
}

# Define signal and background
# Background is a list of processes to add together (2 bkgs = 1 bkg curve representing the sum of the 2 bkgs)
# Signal is a list a signals that will be treated separately (2 signals = 2 graphs/histograms in plot)
processClass = {
    "signal": [
        # varying mMed
        "tchannel_mMed-1000_mDark-20_rinv-0p3_alpha-peak",
        "tchannel_mMed-3000_mDark-20_rinv-0p3_alpha-peak",
        "tchannel_mMed-4000_mDark-20_rinv-0p3_alpha-peak",
        ],
    "background": ["QCD"]
}
N_SIGNAL = len(processClass["signal"])

signalForAuc = "tchannel_mMed-3000_mDark-20_rinv-0p3_alpha-peak"
iSignalForAuc = processClass["signal"].index(signalForAuc)

# Color of the different processes
colors = {
    # varying mMed
    "tchannel_mMed-1000_mDark-20_rinv-0p3_alpha-peak": ROOT.kBlue+1,
    "tchannel_mMed-3000_mDark-20_rinv-0p3_alpha-peak": ROOT.kRed+1,
    "tchannel_mMed-4000_mDark-20_rinv-0p3_alpha-peak": ROOT.kGreen+2,

    # Background
    "QCD": ROOT.TColor.GetColor(155, 152, 204),
}

colors["background"] = ROOT.TColor.GetColor(155, 152, 204)

colors["signal"] = []
for signal in processClass["signal"]:
    colors["signal"].append(colors[signal])


# Subplots to draw
# Available subplots: "histograms", "efficiency", "significance"
panels = ["histograms", "efficiency", "significance"]
