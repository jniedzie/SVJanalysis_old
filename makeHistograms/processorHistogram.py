from coffea import processor
from coffea import hist
from coffea.nanoevents.methods import vector
import awkward as ak
import numpy as np
from collections import OrderedDict
import sys

# Needed so that ak.zip({"pt": [...], "eta": [...], "phi": [...], "mass": [...]},
#                         with_name="PtEtaPhiMLorentzVector")
# is understood as a PtEtaPhiMLorentzVector from coffea.nanoevents.methods.vector
ak.behavior.update(vector.behavior)


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

        ## Variables to histogram and their x label
        self.variables = {
            "nFatJet"             : r"Number of Fat Jet",
            "FatJet_pt"           : r"FatJet $p_{T}$ [GeV]",
            "FatJet_eta"          : r"FatJet $\eta$",
            "FatJet_mass"         : r"FatJet mass [GeV]",
            "FatJet_msoftdrop"    : r"FatJet softdrop mass [GeV]",
            "FatJet_tau1"         : r"FatJet $\tau_{1}$",
            "FatJet_tau2"         : r"FatJet $\tau_{2}$",
            "FatJet_tau3"         : r"FatJet $\tau_{3}$",

            "nJet"                : r"Number of AK4 jets",
            "Jet_pt"              : r"Jet $p_{T}$ [GeV]",
            "Jet_eta"             : r"Jet $\eta$",
            "Jet_mass"            : r"Jet mass [GeV]",

            "MET_pt"              : r"MET [GeV]",
            "HT_AK8"              : r"H_{T}^{AK8}",
            "ST_AK8"              : r"S_{T}^{AK8}",
            "HT_AK4"              : r"H_{T}^{AK4}",
            "ST_AK4"              : r"S_{T}^{AK4}",

            "FatJet_tau21"        : r"FatJet $\tau_{21}$",
            "FatJet_tau32"        : r"FatJet $\tau_{32}$",
            "METrHT_AK8"          : r"MET/H_{T}^{AK8}",
            "METrST_AK8"          : r"MET/H_{T}^{AK8}",
            "METrHT_AK4"          : r"MET/H_{T}^{AK4}",
            "METrST_AK4"          : r"MET/H_{T}^{AK4}",

            "deltaR_J1J2"         : r"$\Delta R$ J1 J2",
            "deltaPhi_J1J2"       : r"$\Delta \Phi$ J1 J2",
            "deltaR_J1_J1PFcands" : r"$\Delta R$ J1 J1PFCands",
        }


        ## In the process method, a genWeights dict is defined, storing the gen weights ak array to use for the different variables
        #  Here we define which gen weights to use for the different variables
        #  Explanation of the different dict values:
        #    * noCut         : no cut nor array broadcasting
        #    * noCut_J       : no cut and array broadcasting for FatJets
        #    * noCut_j       : no cut and array broadcasting for AK4 Jets
        #    * ge1J          : events with at least 1 FatJet
        #    * ge1J_tau21_bc : events with at least 1 FatJet, tau1 != 0, array broadcasting
        #    * ge1J_tau32_bc : events with at least 1 FatJet, tau2 != 0, array broadcasting
        #    * ge1J_J1PFcands: events with at least 1 FatJet, PF candidates for leading FatJet
        #    * ge2J          : events with at least 2 FatJet
        #    * ge1j          : events with at least 1 AK4 Jet
        self.genWeightsInfo = {
            "nFatJet"             : "noCut",
            "FatJet_pt"           : "noCut_J",
            "FatJet_eta"          : "noCut_J",
            "FatJet_mass"         : "noCut_J",
            "FatJet_msoftdrop"    : "noCut_J",
            "FatJet_tau1"         : "noCut_J",
            "FatJet_tau2"         : "noCut_J",
            "FatJet_tau3"         : "noCut_J",

            "nJet"                : "noCut",
            "Jet_pt"              : "noCut_j",
            "Jet_eta"             : "noCut_j",
            "Jet_mass"            : "noCut_j",

            "MET_pt"              : "noCut",
            "HT_AK8"              : "noCut",
            "ST_AK8"              : "noCut",
            "HT_AK4"              : "noCut",
            "ST_AK4"              : "noCut",

            "FatJet_tau21"        : "ge1J_tau21_bc",
            "FatJet_tau32"        : "ge1J_tau32_bc",
            "METrHT_AK8"          : "ge1J",
            "METrST_AK8"          : "ge1J",
            "METrHT_AK4"          : "ge1j",
            "METrST_AK4"          : "ge1j",

            "deltaR_J1J2"         : "ge2J",
            "deltaPhi_J1J2"       : "ge2J",

            "deltaR_J1_J1PFcands" : "ge1J_J1PFcands",
        }


        ## In the process method, a sumGenWeights dict is defined to store the sum of gen weights
        #  to keep track of the different event-level cuts for histogram normalization
        #  Note that the following convention is used:
        #  if cut is in self.genWeightsInfo[variable] then self.sumGenWeightsInfo[variable] = cut
        self.sumGenWeightsInfo = {}
        self.cuts = ["noCut", "ge1J", "ge2J", "ge1j"]
        for variable in self.genWeightsInfo.keys():
            for cut in self.cuts:
                if cut in self.genWeightsInfo[variable]:
                    self.sumGenWeightsInfo[variable] = cut
        

        ## Running some sanity checks
        if self.variables.keys() != self.genWeightsInfo.keys():
            if len(list(self.variables.keys())) > len(list(self.genWeightsInfo.keys())):
                print("ERROR: Some variables have been defined but do have a defined gen weight")
            else:
                print("ERROR: Some variables have a defined gen weight but have not been defined")
            print("Defined variables:")
            print(sorted(list(self.variables.keys())))
            print("Defined gen weights:")
            print(sorted(list(self.genWeightsInfo.keys())))
            sys.exit()

        if self.variables.keys() != self.sumGenWeightsInfo.keys():
            if len(list(self.variables.keys())) > len(list(self.sumGenWeightsInfo.keys())):
                print("ERROR: Some variables have been defined but do have a defined sum gen weight")
            else:
                print("ERROR: Some variables have a defined sum gen weight but have not been defined")
            print("Defined variables:")
            print(sorted(list(self.variables.keys())))
            print("Defined sun gen weights:")
            print(sorted(list(self.sumGenWeightsInfo.keys())))
            sys.exit()
        

        ## Make dict filled with coffea histograms to be used in the accumulator
        histograms = {}

        # List of regexes for binning definition
        regexes = list(binning["regex"].keys())

        # Loop over all varaibles to histogram
        for variable, label in self.variables.items():
            # Get binning of the variable
            if variable not in binning["noregex"]:
                indices = utl.inregex(variable, regexes)
                if len(indices) == 0:
                    print("ERROR: Binning of %s is not defined." %variable)
                    sys.exit()
                elif len(indices) > 1:
                    if VERBOSE>0: print("ERROR: %s matches several regexes. Binning cannot be defined." %variable)
                    sys.exit()
                else: 
                    binning_ = binning["regex"][regexes[indices[0]]]
            else:
                binning_ = binning["noregex"][variable]

            # Define coffea histogram for the variable
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

        ## Define accumulator
        output = self.accumulator.identity()


        ## Define dict storing the different ak arrays that will be used to fill the histograms
        akArrays = {}


        ## Compute variables to histogram / defined in constructor
        # No cut
        akArrays["nFatJet"] = events["nFatJet"]
        akArrays["FatJet_pt"] = ak.flatten(events["FatJet_pt"])
        akArrays["FatJet_eta"] = ak.flatten(events["FatJet_eta"])
        akArrays["FatJet_mass"] = ak.flatten(events["FatJet_mass"])
        akArrays["FatJet_msoftdrop"] = ak.flatten(events["FatJet_msoftdrop"])
        akArrays["FatJet_tau1"] = ak.flatten(events["FatJet_tau1"])
        akArrays["FatJet_tau2"] = ak.flatten(events["FatJet_tau2"])
        akArrays["FatJet_tau3"] = ak.flatten(events["FatJet_tau3"])

        akArrays["nJet"] = events["nJet"]
        akArrays["Jet_pt"] = ak.flatten(events["Jet_pt"])
        akArrays["Jet_eta"] = ak.flatten(events["Jet_eta"])
        akArrays["Jet_mass"] = ak.flatten(events["Jet_mass"])

        akArrays["MET_pt"] = events["MET_pt"]

        akArrays["HT_AK8"] = ak.sum(events["FatJet_pt"], axis=1)
        akArrays["ST_AK8"] = akArrays["HT_AK8"] + events["MET_pt"]

        akArrays["HT_AK4"] = ak.sum(events["Jet_pt"], axis=1)
        akArrays["ST_AK4"] = akArrays["HT_AK4"] + events["MET_pt"]


        # At least 1 FatJet
        mask_ge1J = (akArrays["nFatJet"]>0)
        mask_ge1J_bc = ak.flatten(ak.broadcast_arrays(mask_ge1J, events["FatJet_pt"])[0])
        mask_ge1J_tau21_bc = mask_ge1J_bc & (akArrays["FatJet_tau1"][mask_ge1J_bc] != 0.)
        mask_ge1J_tau32_bc = mask_ge1J_bc & (akArrays["FatJet_tau2"][mask_ge1J_bc] != 0.)

        # The two lines below do not work if tau1 or tau2 is 0, which can apparently happen!
        #akArrays["FatJet_tau21"] = akArrays["FatJet_tau2"][mask_ge1J_bc]/akArrays["FatJet_tau1"][mask_ge1J_bc]
        #akArrays["FatJet_tau32"] = akArrays["FatJet_tau3"][mask_ge1J_bc]/akArrays["FatJet_tau2"][mask_ge1J_bc]
        akArrays["FatJet_tau21"] = akArrays["FatJet_tau2"][mask_ge1J_tau21_bc]/akArrays["FatJet_tau1"][mask_ge1J_tau21_bc]
        akArrays["FatJet_tau32"] = akArrays["FatJet_tau3"][mask_ge1J_tau32_bc]/akArrays["FatJet_tau2"][mask_ge1J_tau32_bc]

        akArrays["METrHT_AK8"] = akArrays["MET_pt"][mask_ge1J] / akArrays["HT_AK8"][mask_ge1J]
        akArrays["METrST_AK8"] = akArrays["MET_pt"][mask_ge1J] / akArrays["ST_AK8"][mask_ge1J]

        FatJet_ge1J = ak.zip(
            {
                "pt":   events["FatJet_pt"][mask_ge1J],
                "eta":  events["FatJet_eta"][mask_ge1J],
                "phi":  events["FatJet_phi"][mask_ge1J],
                "mass": events["FatJet_mass"][mask_ge1J],
            },
            with_name="PtEtaPhiMLorentzVector",
        )
        FatJet_ge1J_J1 = FatJet_ge1J[:,0]

        # This only works for 102X - Will be improved in a next iteration
        mask_fatJetIdx0 = (events["FatJetPFCands_jetIdx"][mask_ge1J] == 0)
        J1PFCands = ak.zip(
            {
                "pt":   events["FatJetPFCands_pt"][mask_ge1J][mask_fatJetIdx0],
                "eta":  events["FatJetPFCands_eta"][mask_ge1J][mask_fatJetIdx0],
                "phi":  events["FatJetPFCands_phi"][mask_ge1J][mask_fatJetIdx0],
                "mass": events["FatJetPFCands_mass"][mask_ge1J][mask_fatJetIdx0],
            },
            with_name="PtEtaPhiMLorentzVector",
        )

        FatJet_ge1J_J1_bc = ak.broadcast_arrays(FatJet_ge1J_J1, J1PFCands)[0]
        deltaR_J1_J1PFcands = J1PFCands.delta_r(FatJet_ge1J_J1_bc)
        akArrays["deltaR_J1_J1PFcands"] = ak.flatten(deltaR_J1_J1PFcands)


        # At least 1 AK4 Jet
        mask_ge1j = (akArrays["nJet"]>0)

        akArrays["METrHT_AK4"] = akArrays["MET_pt"][mask_ge1j] / akArrays["HT_AK4"][mask_ge1j]
        akArrays["METrST_AK4"] = akArrays["MET_pt"][mask_ge1j] / akArrays["ST_AK4"][mask_ge1j]


        # At least 2 FatJets
        mask_ge2J = (akArrays["nFatJet"]>1)
        FatJet_ge2J = ak.zip(
            {
                "pt":   events["FatJet_pt"][mask_ge2J],
                "eta":  events["FatJet_eta"][mask_ge2J],
                "phi":  events["FatJet_phi"][mask_ge2J],
                "mass": events["FatJet_mass"][mask_ge2J],
            },
            with_name="PtEtaPhiMLorentzVector",
        )

        akArrays["deltaR_J1J2"] = FatJet_ge2J[:,0].delta_r(FatJet_ge2J[:,1])
        akArrays["deltaPhi_J1J2"] = FatJet_ge2J[:,0].delta_phi(FatJet_ge2J[:,1])


        ## Running some sanity checks
        if self.variables.keys() != akArrays.keys():
            if len(list(variables.keys())) > len(list(self.akArrays.keys())):
                print("ERROR: Some variables have been defined but no ak array has been made")
            else:
                print("ERROR: Some variables have a defined ak array but have not been defined")
            print("Defined variables:")
            print(sorted(list(self.variables.keys())))
            print("Defined akArrays variables:")
            print(sorted(list(akArrays.keys())))
            sys.exit()


        ## Make gen weights
        #  Here we make the gen weights defined in self.genWeights
        genWeights = {}
        genWeights["noCut"] = events["genWeight"]
        genWeights["noCut_J"] = ak.flatten(ak.broadcast_arrays(genWeights["noCut"], events["FatJet_pt"])[0])
        genWeights["noCut_j"] = ak.flatten(ak.broadcast_arrays(genWeights["noCut"], events["Jet_pt"])[0])
        genWeights["ge1J"] = genWeights["noCut"][mask_ge1J]
        genWeights["ge1J_bc"] = genWeights["noCut_J"][mask_ge1J_bc]
        genWeights["ge1J_tau21_bc"] = genWeights["noCut_J"][mask_ge1J_tau21_bc]
        genWeights["ge1J_tau32_bc"] = genWeights["noCut_J"][mask_ge1J_tau32_bc]
        genWeights["ge1J_J1PFcands"] = ak.flatten(ak.broadcast_arrays(genWeights["ge1J"], events["FatJetPFCands_pt"][mask_ge1J][mask_fatJetIdx0])[0])
        genWeights["ge2J"] = genWeights["noCut"][mask_ge2J]
        genWeights["ge1j"] = genWeights["noCut"][mask_ge1j]



        ## Make sum of gen weights
        #  Here we make the sum of gen weights defined in self.sumGenWeights
        sumGenWeights = {}
        for cut in self.cuts:
            sumGenWeights[cut] = ak.sum(genWeights[cut])


        ## Fill histograms
        for variable in akArrays.keys():
            fillHistogram(output, variable, akArrays[variable], genWeights[self.genWeightsInfo[variable]])


        ## Cutflow information for histogram normalization
        for variable in akArrays.keys():
            output["cutflow"][variable] += sumGenWeights[self.sumGenWeightsInfo[variable]]
        for cut in self.cuts:
            output["cutflow"][cut] += sumGenWeights[cut]


        return output


    def postprocess(self, accumulator):
        return accumulator

