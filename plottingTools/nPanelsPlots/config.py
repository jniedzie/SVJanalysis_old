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
    var_labels = json.load(f)
    # Quick hack
    var_labels = { x: x for x in sorted(var_labels) }


# Processes to plot
path_to_histograms_tchannel = "/work/fleble/t_channel_histograms/skim1/"
path_to_histograms_QCD = "/work/fleble/QCD_histograms/skim1/"
processes = {
    "tchannel_mMed-1000_mDark-20_rinv-0p3_alpha-peak": [ path_to_histograms_tchannel + "tchannel_mMed-1000_mDark-20_rinv-0.3_alpha-peak.root" ],
    "tchannel_mMed-3000_mDark-20_rinv-0p3_alpha-peak": [ path_to_histograms_tchannel + "tchannel_mMed-3000_mDark-20_rinv-0.3_alpha-peak.root" ],
    "tchannel_mMed-4000_mDark-20_rinv-0p3_alpha-peak": [ path_to_histograms_tchannel + "tchannel_mMed-4000_mDark-20_rinv-0.3_alpha-peak.root" ],

    "QCD": [ 
         path_to_histograms_QCD + "QCD_170_300.root",
         path_to_histograms_QCD + "QCD_300_470.root",
         path_to_histograms_QCD + "QCD_470_600.root",
         path_to_histograms_QCD + "QCD_600_800.root",
         path_to_histograms_QCD + "QCD_800_1000.root",
         path_to_histograms_QCD + "QCD_1000_1400.root",
         path_to_histograms_QCD + "QCD_1400_1800.root",
         path_to_histograms_QCD + "QCD_1800_2400.root",
         path_to_histograms_QCD + "QCD_2400_3200.root",
         path_to_histograms_QCD + "QCD_3200_Inf.root"
         ]
}

# Label for the processes (legend)
process_labels = {
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
process_class = {
    "signal": [
        # varying mMed
        "tchannel_mMed-1000_mDark-20_rinv-0p3_alpha-peak",
        "tchannel_mMed-3000_mDark-20_rinv-0p3_alpha-peak",
        "tchannel_mMed-4000_mDark-20_rinv-0p3_alpha-peak",
        ],
    "background": ["QCD"]
}
N_SIGNAL = len(process_class["signal"])

SIGNAL_FOR_AUC = "tchannel_mMed-3000_mDark-20_rinv-0p3_alpha-peak"
ISIGNAL_FOR_AUC = process_class["signal"].index(SIGNAL_FOR_AUC)

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
for signal in process_class["signal"]:
    colors["signal"].append(colors[signal])


# Subplots to draw
# Available subplots: "histograms", "efficiency", "significance"
panels = ["histograms", "efficiency", "significance"]
