from coffea import processor
from coffea import hist
import awkward as ak
import numpy as np
from collections import OrderedDict
import sys

sys.path.append("../pythonUtils/")
import utilities as utl



def fillHistogram(output, variable, axis, weight):
    """Fill histogram if the variable has been defined."""
    if variable in output.keys():
        output[variable].fill(axis=axis, weight=weight)
    else:
        pass
        # A warning has been sent earlier when defining coffea histograms


class Histogram1(processor.ProcessorABC):
    """
    Make histogram of selected variables in input ROOT file.
    If variable do not exist in ROOT file, they are computed on the fly (e.g. HT, MET/HT, ...).
    """

    def __init__(self, binning, VERBOSE=1):
        """Define the variables to histogram, their binning and their coffea histogram object."""

        ## Variables to histogram and and their x label
        variables = {
            "nFatJet"          : r"Number of Fat Jet",
            "FatJet_pt"        : r"FatJet $p_{T}$ [GeV]",
            "FatJet_eta"       : r"FatJet $\eta$",
            "FatJet_mass"      : r"FatJet mass [GeV]",
            "FatJet_msoftdrop" : r"FatJet softdrop mass [GeV]",
            "FatJet_tau1"      : r"FatJet $\tau_{1}$",
            "FatJet_tau2"      : r"FatJet $\tau_{2}$",
            "FatJet_tau3"      : r"FatJet $\tau_{3}$",

            "nJet"             : r"Number of AK4 jets",
            "Jet_pt"           : r"Jet $p_{T}$ [GeV]",
            "Jet_eta"          : r"Jet $\eta$",
            "Jet_mass"         : r"Jet mass [GeV]",

            "MET_pt"           : r"MET [GeV]",
            "HT_AK8"           : r"H_{T}^{AK8}",
            "ST_AK8"           : r"S_{T}^{AK8}",
            "HT_AK4"           : r"H_{T}^{AK4}",
            "ST_AK4"           : r"S_{T}^{AK4}",

            "FatJet_tau21"     : r"FatJet $\tau_{21}$",
            "FatJet_tau32"     : r"FatJet $\tau_{32}$",
            "METrHT_AK8"       : r"MET/H_{T}^{AK8}",
            "METrST_AK8"       : r"MET/H_{T}^{AK8}",
            "METrHT_AK4"       : r"MET/H_{T}^{AK4}",
            "METrST_AK4"       : r"MET/H_{T}^{AK4}",
        }

        
        ## Make dict filled with coffea histograms
        histograms = {}

        # List of regexes for binning definition
        regexes = list(binning["regex"].keys())

        # Loop over all varaibles to histogram
        for variable, label in variables.items():
            # Get binning of the variable
            if variable not in binning["noregex"]:
                indices = utl.inregex(variable, regexes)
                if len(indices) == 0:
                    if VERBOSE>0: print("WARNING: Binning of %s is not defined. Skipping." %variable)
                    continue
                elif len(indices) > 1:
                    if VERBOSE>0: print("WARNING: %s matches several regexes. Binning cannot be defined. Skipping." %variable)
                    continue
                else: 
                    binning_ = binning["regex"][regexes[indices[0]]]
            else:
                binning_ = binning["noregex"][variable]

            # Define coffea histogra for the variable
            histograms[variable] = hist.Hist(variable, hist.Bin("axis", label, binning_[0], binning_[1], binning_[2]) )

        ## Define accumulator
        self._accumulator = processor.dict_accumulator({
            **histograms,
            **{
                # Cutflow in case cuts are applied (e.g. for delta phi between leding jets, at least 2 jets are needed!)
                'cutflow': processor.defaultdict_accumulator(float)
            }
            })


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        """Fill histograms."""

        output = self.accumulator.identity()


        ## Compute variables defined in constructor
        # No cut
        nFatJet = events["nFatJet"]
        FatJet_pt = ak.flatten(events["FatJet_pt"])
        FatJet_eta = ak.flatten(events["FatJet_eta"])
        FatJet_mass = ak.flatten(events["FatJet_mass"])
        FatJet_msoftdrop = ak.flatten(events["FatJet_msoftdrop"])
        FatJet_tau1 = ak.flatten(events["FatJet_tau1"])
        FatJet_tau2 = ak.flatten(events["FatJet_tau2"])
        FatJet_tau3 = ak.flatten(events["FatJet_tau3"])

        nJet = events["nJet"]
        Jet_pt = ak.flatten(events["Jet_pt"])
        Jet_eta = ak.flatten(events["Jet_eta"])
        Jet_mass = ak.flatten(events["Jet_mass"])

        MET_pt = events["MET_pt"]

        HT_AK8 = ak.sum(events["FatJet_pt"], axis=1)
        ST_AK8 = HT_AK8 + events["MET_pt"]

        HT_AK4 = ak.sum(events["Jet_pt"], axis=1)
        ST_AK4 = HT_AK4 + events["MET_pt"]

        # At least 1 FatJet
        mask_ge1J = (nFatJet>0)
        mask_ge1J_bc = ak.flatten(ak.broadcast_arrays(mask_ge1J, events["FatJet_pt"])[0])
        mask_ge1J_tau21_bc = mask_ge1J_bc & (FatJet_tau1[mask_ge1J_bc] != 0.)
        mask_ge1J_tau32_bc = mask_ge1J_bc & (FatJet_tau2[mask_ge1J_bc] != 0.)

        # The two lines below do not work if tau1 or tau2 is 0, which can apparently happen!
        #FatJet_tau21 = FatJet_tau2[mask_ge1J_bc]/FatJet_tau1[mask_ge1J_bc]
        #FatJet_tau32 = FatJet_tau3[mask_ge1J_bc]/FatJet_tau2[mask_ge1J_bc]
        FatJet_tau21 = FatJet_tau2[mask_ge1J_tau21_bc]/FatJet_tau1[mask_ge1J_tau21_bc]
        FatJet_tau32 = FatJet_tau3[mask_ge1J_tau32_bc]/FatJet_tau2[mask_ge1J_tau32_bc]

        METrHT_AK8 = MET_pt[mask_ge1J] / HT_AK8[mask_ge1J]
        METrST_AK8 = MET_pt[mask_ge1J] / ST_AK8[mask_ge1J]

        # At least 1 AK4 Jet
        mask_ge1j = (nJet>0)

        METrHT_AK4 = MET_pt[mask_ge1j] / HT_AK4[mask_ge1j]
        METrST_AK4 = MET_pt[mask_ge1j] / ST_AK4[mask_ge1j]


        ## Define gen weights
        genWeight = events["genWeight"]
        genWeight_J = ak.flatten(ak.broadcast_arrays(events["genWeight"], events["FatJet_pt"])[0])
        genWeight_j = ak.flatten(ak.broadcast_arrays(events["genWeight"], events["Jet_pt"])[0])
        genWeight_ge1J = genWeight[mask_ge1J]
        genWeight_ge1j = genWeight[mask_ge1j]
        genWeight_ge1J_bc = genWeight_J[mask_ge1J_bc]
        genWeight_ge1J_tau21_bc = genWeight_J[mask_ge1J_tau21_bc]
        genWeight_ge1J_tau32_bc = genWeight_J[mask_ge1J_tau32_bc]

        sumGenWeight = ak.sum(genWeight)
        sumGenWeight_ge1J = ak.sum(genWeight_ge1J)
        sumGenWeight_ge1j = ak.sum(genWeight_ge1j)


        ## Fill histograms
        fillHistogram(output, "nFatJet", nFatJet, genWeight)
        fillHistogram(output, "FatJet_pt", FatJet_pt, genWeight_J)
        fillHistogram(output, "FatJet_eta", FatJet_eta, genWeight_J)
        fillHistogram(output, "FatJet_mass", FatJet_mass, genWeight_J)
        fillHistogram(output, "FatJet_msoftdrop", FatJet_msoftdrop, genWeight_J)
        fillHistogram(output, "FatJet_tau1", FatJet_tau1, genWeight_J)
        fillHistogram(output, "FatJet_tau2", FatJet_tau2, genWeight_J)
        fillHistogram(output, "FatJet_tau3", FatJet_tau3, genWeight_J)

        fillHistogram(output, "nJet", nJet, genWeight)
        fillHistogram(output, "Jet_pt", Jet_pt, genWeight_j)
        fillHistogram(output, "Jet_eta", Jet_eta, genWeight_j)
        fillHistogram(output, "Jet_mass", Jet_mass, genWeight_j)

        fillHistogram(output, "MET_pt", MET_pt, genWeight)

        fillHistogram(output, "HT_AK8", HT_AK8, genWeight)
        fillHistogram(output, "ST_AK8", ST_AK8, genWeight)
        fillHistogram(output, "HT_AK4", HT_AK4, genWeight)
        fillHistogram(output, "ST_AK4", ST_AK4, genWeight)

        fillHistogram(output, "FatJet_tau21", FatJet_tau21, genWeight_ge1J_bc)
        fillHistogram(output, "FatJet_tau32", FatJet_tau32, genWeight_ge1J_tau32_bc)

        fillHistogram(output, "METrHT_AK8", METrHT_AK8, genWeight_ge1J)
        fillHistogram(output, "METrST_AK8", METrST_AK8, genWeight_ge1J)

        fillHistogram(output, "METrHT_AK4", METrHT_AK4, genWeight_ge1j)
        fillHistogram(output, "METrST_AK4", METrST_AK4, genWeight_ge1j)


        ## Cutflow information for histogram normalization

        # List of variables depending on cuts/weights
        variables_noCut = [
            "all",
            "nFatJet", "FatJet_pt", "FatJet_eta", "FatJet_mass", "FatJet_msoftdrop",
            "FatJet_tau1", "FatJet_tau2", "FatJet_tau3",
            "nJet", "Jet_pt", "Jet_eta", "Jet_mass",
            "MET_pt", "HT_AK8", "ST_AK8", "HT_AK4", "ST_AK4"
        ]
        variables_ge1j = ["ge1j", "METrHT_AK4", "METrST_AK4"]
        variables_ge1J = ["ge1J", "FatJet_tau21", "FatJet_tau32", "METrHT_AK8", "METrST_AK8"]

        for variable in variables_noCut:
            output["cutflow"][variable] += sumGenWeight

        for variable in variables_ge1j:
            output["cutflow"][variable] += sumGenWeight_ge1j

        for variable in variables_ge1J:
            output["cutflow"][variable] += sumGenWeight_ge1J

        return output


    def postprocess(self, accumulator):
        return accumulator

