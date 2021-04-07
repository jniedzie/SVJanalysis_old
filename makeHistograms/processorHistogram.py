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


def get_binning(variable, binning):
    """
    Return binning for variable in 1st argument providing the binning information in 2nd argument.
    binning = {
        "noregex": {
            "j1j2_mass"         :  [500  ,  0   , 5000 ],
        },
        "regex": {
            "n(Fat)?Jet"        :  [20   ,  0   , 20   ],
        }
    }
    """

    # List of regexes for binning definition
    regexes = list(binning["regex"].keys())

    if variable not in binning["noregex"]:
        indices = utl.inregex(variable, regexes)
        if len(indices) == 0:
            print("ERROR: Binning of %s is not defined." %variable)
            sys.exit()
        elif len(indices) > 1:
            print("ERROR: %s matches several regexes. Binning cannot be defined." %variable)
            sys.exit()
        else: 
            binning_ = binning["regex"][regexes[indices[0]]]
    else:
        binning_ = binning["noregex"][variable]

    return binning_


def make_PtEtaPhiMLorentzVector(pt, eta, phi, mass):
    """Take pt, eta, phi, mass awkward arrays and return the corresponding PtEtaPhiMLorentzVector."""

    vec = ak.zip(
        {
            "pt": pt,
            "eta": eta,
            "phi": phi,
            "mass": mass,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    return vec



def fillHistogram(output, variable, axis, weight):
    """Fill histogram if the variable has been defined."""

    if variable in output.keys():
        output[variable].fill(axis=axis, weight=weight)
    else:
        pass
        # A warning has been sent earlier when defining coffea histograms


def check_dict_keys(dict1, dict2, dict1Name, dict2Name):
    """Check if dict1 and dict2 have the same keys. Interrupt program if not the case."""

    if dict1.keys() != dict2.keys():
        if len(list(dict1.keys())) > len(list(dict2.keys())):
            print("ERROR: Some %s have been defined but do not have a defined %s" %(dict1Name, dict2Name))
        else:
            print("ERROR: Some %s have a defined %s but have not been defined" %(dict1Name, dict2Name))
        print("Defined %s:" %dict1Name)
        print(sorted(list(dict1.keys())))
        print("Defined %s:" %dict2Name)
        print(sorted(list(dict2.keys())))
        sys.exit()

    return
 

class Histogram1(processor.ProcessorABC):
    """
    Make histogram of selected variables in input ROOT file.
    If variable do not exist in ROOT file, they are computed on the fly (e.g. HT, MET/HT, ...).

    Conventions used in this class:
       * j denotes AK4 jets
       * J denotes AK8 jets
       * j1 (resp. J1) denotes leading AK4 (resp. AK8) jet
       * ge{n}j denotes a quantity with a cut selecting events with at least {n} AK4 jets
    """


    def __init__(self, binning, fileType, VERBOSE=1):
        """Define the variables to histogram, their binning and their coffea histogram object."""

        self.fileType = fileType

        ## Variables to histogram and their x label
        self.variables = {
            "nFatJet"                    : r"Number of Fat Jet",
            "FatJet_pt"                  : r"FatJet $p_{T}$ [GeV]",
            "FatJet_eta"                 : r"FatJet $\eta$",
            "FatJet_mass"                : r"FatJet mass [GeV]",
            "FatJet_msoftdrop"           : r"FatJet softdrop mass [GeV]",
            "FatJet_tau1"                : r"FatJet $\tau_{1}$",
            "FatJet_tau2"                : r"FatJet $\tau_{2}$",
            "FatJet_tau3"                : r"FatJet $\tau_{3}$",

            "nJet"                       : r"Number of AK4 jets",
            "Jet_pt"                     : r"Jet $p_{T}$ [GeV]",
            "Jet_eta"                    : r"Jet $\eta$",
            "Jet_mass"                   : r"Jet mass [GeV]",

            "MET_pt"                     : r"MET [GeV]",
            "HT_AK8"                     : r"H_{T}^{AK8}",
            "ST_AK8"                     : r"S_{T}^{AK8}",
            "HT_AK4"                     : r"H_{T}^{AK4}",
            "ST_AK4"                     : r"S_{T}^{AK4}",

            "FatJet_tau21"               : r"FatJet $\tau_{21}$",
            "FatJet_tau32"               : r"FatJet $\tau_{32}$",
            "METrHT_AK8"                 : r"MET/H_{T}^{AK8}",
            "METrST_AK8"                 : r"MET/H_{T}^{AK8}",
            "METrHT_AK4"                 : r"MET/H_{T}^{AK4}",
            "METrST_AK4"                 : r"MET/H_{T}^{AK4}",

            "deltaR_J1J2"                : r"$\Delta$ R J1 J2",
            "deltaPhi_J1J2"              : r"$\Delta$ $\Phi$ J1 J2",

            "deltaR_J1_J1PFcands"        : r"$\Delta$ R J1 J1PFCands",
            "sum_J1PFcandsPt_dRle0.1"    : r"Sum J1PFCandsPt $\Delta$ R $\leq$ 0.1",
            "sum_J1PFcandsPt_dR0.1To0.2" : r"0.1 $\lg$ Sum J1PFCandsPt $\Delta$ R $\leq$ 0.2",
            "sum_J1PFcandsPt_dR0.2To0.3" : r"0.2 $\lg$ Sum J1PFCandsPt $\Delta$ R $\leq$ 0.3",
            "sum_J1PFcandsPt_dR0.3To0.4" : r"0.3 $\lg$ Sum J1PFCandsPt $\Delta$ R $\leq$ 0.4",
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
            "nFatJet"                    : "noCut",
            "FatJet_pt"                  : "noCut_J",
            "FatJet_eta"                 : "noCut_J",
            "FatJet_mass"                : "noCut_J",
            "FatJet_msoftdrop"           : "noCut_J",
            "FatJet_tau1"                : "noCut_J",
            "FatJet_tau2"                : "noCut_J",
            "FatJet_tau3"                : "noCut_J",

            "nJet"                       : "noCut",
            "Jet_pt"                     : "noCut_j",
            "Jet_eta"                    : "noCut_j",
            "Jet_mass"                   : "noCut_j",

            "MET_pt"                     : "noCut",
            "HT_AK8"                     : "noCut",
            "ST_AK8"                     : "noCut",
            "HT_AK4"                     : "noCut",
            "ST_AK4"                     : "noCut",

            "FatJet_tau21"               : "ge1J_tau21_bc",
            "FatJet_tau32"               : "ge1J_tau32_bc",
            "METrHT_AK8"                 : "ge1J",
            "METrST_AK8"                 : "ge1J",
            "METrHT_AK4"                 : "ge1j",
            "METrST_AK4"                 : "ge1j",

            "deltaR_J1J2"                : "ge2J",
            "deltaPhi_J1J2"              : "ge2J",

            "deltaR_J1_J1PFcands"        : "ge1J_J1PFcands",
            "sum_J1PFcandsPt_dRle0.1"    : "ge1J",
            "sum_J1PFcandsPt_dR0.1To0.2" : "ge1J",
            "sum_J1PFcandsPt_dR0.2To0.3" : "ge1J",
            "sum_J1PFcandsPt_dR0.3To0.4" : "ge1J",
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
        check_dict_keys(self.variables, self.genWeightsInfo, "variable", "gen weight")
        check_dict_keys(self.variables, self.sumGenWeightsInfo, "variable", "sum gen weight")


        ## Make dict filled with coffea histograms to be used in the accumulator
        histograms = {}

        # Loop over all varaibles to histogram
        for variable, label in self.variables.items():
            # Get binning of the variable
            binning_ = get_binning(variable, binning)
            # Define coffea histogram for the variable
            histograms[variable] = hist.Hist(variable, hist.Bin("axis", label, binning_[0], binning_[1], binning_[2]))


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


    def initialize_ak_arrays(self, events):
        """Fill a dict with the ak arrays that will be used to fill the histograms."""

        ## Define dict storing the different ak arrays
        akArrays = {}

        ## Define dict storing the different masks used
        masks = {}


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
        masks["ge1J"] = (akArrays["nFatJet"]>0)
        masks["ge1J_bc"] = ak.flatten(ak.broadcast_arrays(masks["ge1J"], events["FatJet_pt"])[0])
        masks["ge1J_tau21_bc"] = masks["ge1J_bc"] & (akArrays["FatJet_tau1"][masks["ge1J_bc"]] != 0.)
        masks["ge1J_tau32_bc"] = masks["ge1J_bc"] & (akArrays["FatJet_tau2"][masks["ge1J_bc"]] != 0.)

        # The two lines below do not work if tau1 or tau2 is 0, which can apparently happen!
        #akArrays["FatJet_tau21"] = akArrays["FatJet_tau2"][masks["ge1J_bc]/akArrays["FatJet_tau1"][masks["ge1J_bc]
        #akArrays["FatJet_tau32"] = akArrays["FatJet_tau3"][masks["ge1J_bc]/akArrays["FatJet_tau2"][masks["ge1J_bc]
        akArrays["FatJet_tau21"] = akArrays["FatJet_tau2"][masks["ge1J_tau21_bc"]]/akArrays["FatJet_tau1"][masks["ge1J_tau21_bc"]]
        akArrays["FatJet_tau32"] = akArrays["FatJet_tau3"][masks["ge1J_tau32_bc"]]/akArrays["FatJet_tau2"][masks["ge1J_tau32_bc"]]

        akArrays["METrHT_AK8"] = akArrays["MET_pt"][masks["ge1J"]] / akArrays["HT_AK8"][masks["ge1J"]]
        akArrays["METrST_AK8"] = akArrays["MET_pt"][masks["ge1J"]] / akArrays["ST_AK8"][masks["ge1J"]]

        fatJet_ge1J = make_PtEtaPhiMLorentzVector(
            events["FatJet_pt"][masks["ge1J"]],
            events["FatJet_eta"][masks["ge1J"]],
            events["FatJet_phi"][masks["ge1J"]],
            events["FatJet_mass"][masks["ge1J"]],
        )
        fatJet_ge1J_J1 = fatJet_ge1J[:,0]


        if self.fileType == "PFnano102X":
            masks["fatJetIdx0"] = (events["FatJetPFCands_jetIdx"][masks["ge1J"]] == 0)
            J1PFCands = make_PtEtaPhiMLorentzVector(
                events["FatJetPFCands_pt"][masks["ge1J"]][masks["fatJetIdx0"]],
                events["FatJetPFCands_eta"][masks["ge1J"]][masks["fatJetIdx0"]],
                events["FatJetPFCands_phi"][masks["ge1J"]][masks["fatJetIdx0"]],
                events["FatJetPFCands_mass"][masks["ge1J"]][masks["fatJetIdx0"]],
            )

        elif self.fileType == "PFnano106X":
            masks["fatJetIdx0"] = (events["JetPFCandsAK8_jetIdx"][masks["ge1J"]] == 0)
            masks["fatJetCandIdx0"] = events["JetPFCandsAK8_candIdx"][masks["ge1J"]][masks["fatJetIdx0"]]
            J1PFCands = make_PtEtaPhiMLorentzVector(
                events["JetPFCands_pt"][masks["ge1J"]][masks["fatJetCandIdx0"]],
                events["JetPFCands_eta"][masks["ge1J"]][masks["fatJetCandIdx0"]],
                events["JetPFCands_phi"][masks["ge1J"]][masks["fatJetCandIdx0"]],
                events["JetPFCands_mass"][masks["ge1J"]][masks["fatJetCandIdx0"]],
            )

        # the else case cannot happen, it has already been tackled

        fatJet_ge1J_J1_bc = ak.broadcast_arrays(fatJet_ge1J_J1, J1PFCands)[0]
        deltaR_J1_J1PFcands = J1PFCands.delta_r(fatJet_ge1J_J1_bc)
        akArrays["deltaR_J1_J1PFcands"] = ak.flatten(deltaR_J1_J1PFcands)

        masks["deltaR_J1_J1PFcands_le0.1"] = (deltaR_J1_J1PFcands <= 0.1)
        masks["deltaR_J1_J1PFcands_0.1To0.2"] = (deltaR_J1_J1PFcands > 0.1) & (deltaR_J1_J1PFcands <= 0.2)
        masks["deltaR_J1_J1PFcands_0.2To0.3"] = (deltaR_J1_J1PFcands > 0.2) & (deltaR_J1_J1PFcands <= 0.3)
        masks["deltaR_J1_J1PFcands_0.3To0.4"] = (deltaR_J1_J1PFcands > 0.3) & (deltaR_J1_J1PFcands <= 0.4)
        akArrays["sum_J1PFcandsPt_dRle0.1"] = ak.sum(J1PFCands.pt[masks["deltaR_J1_J1PFcands_le0.1"]], axis=1)
        akArrays["sum_J1PFcandsPt_dR0.1To0.2"] = ak.sum(J1PFCands.pt[masks["deltaR_J1_J1PFcands_0.1To0.2"]], axis=1)
        akArrays["sum_J1PFcandsPt_dR0.2To0.3"] = ak.sum(J1PFCands.pt[masks["deltaR_J1_J1PFcands_0.2To0.3"]], axis=1)
        akArrays["sum_J1PFcandsPt_dR0.3To0.4"] = ak.sum(J1PFCands.pt[masks["deltaR_J1_J1PFcands_0.3To0.4"]], axis=1)


        # At least 1 AK4 Jet
        masks["ge1j"] = (akArrays["nJet"]>0)

        akArrays["METrHT_AK4"] = akArrays["MET_pt"][masks["ge1j"]] / akArrays["HT_AK4"][masks["ge1j"]]
        akArrays["METrST_AK4"] = akArrays["MET_pt"][masks["ge1j"]] / akArrays["ST_AK4"][masks["ge1j"]]


        # At least 2 FatJets
        masks["ge2J"] = (akArrays["nFatJet"]>1)
        FatJet_ge2J = make_PtEtaPhiMLorentzVector(
            events["FatJet_pt"][masks["ge2J"]],
            events["FatJet_eta"][masks["ge2J"]],
            events["FatJet_phi"][masks["ge2J"]],
            events["FatJet_mass"][masks["ge2J"]],
        )

        akArrays["deltaR_J1J2"] = FatJet_ge2J[:,0].delta_r(FatJet_ge2J[:,1])
        akArrays["deltaPhi_J1J2"] = FatJet_ge2J[:,0].delta_phi(FatJet_ge2J[:,1])

        return akArrays, masks


    def process(self, events):
        """Fill histograms."""

        ## Define accumulator
        output = self.accumulator.identity()

        ## Make akArrays that will be used to fill histograms
        #  Also return masks used, will be useful for making gen weights
        akArrays, masks = self.initialize_ak_arrays(events)

        ## Running some sanity checks
        check_dict_keys(self.variables, akArrays, "variable", "ak array")

        ## Make gen weights
        #  Here we make the gen weights defined in self.genWeights
        genWeights = {}
        genWeights["noCut"] = events["genWeight"]
        genWeights["noCut_J"] = ak.flatten(ak.broadcast_arrays(genWeights["noCut"], events["FatJet_pt"])[0])
        genWeights["noCut_j"] = ak.flatten(ak.broadcast_arrays(genWeights["noCut"], events["Jet_pt"])[0])
        genWeights["ge1J"] = genWeights["noCut"][masks["ge1J"]]
        genWeights["ge1J_tau21_bc"] = genWeights["noCut_J"][masks["ge1J_tau21_bc"]]
        genWeights["ge1J_tau32_bc"] = genWeights["noCut_J"][masks["ge1J_tau32_bc"]]
        if self.fileType == "PFnano102X":
            genWeights["ge1J_J1PFcands"] = ak.flatten(ak.broadcast_arrays(genWeights["ge1J"], events["FatJetPFCands_pt"][masks["ge1J"]][masks["fatJetIdx0"]])[0])
        elif self.fileType == "PFnano106X":
            genWeights["ge1J_J1PFcands"] = ak.flatten(ak.broadcast_arrays(genWeights["ge1J"], events["JetPFCands_pt"][masks["ge1J"]][masks["fatJetCandIdx0"]])[0])
        #genWeights["ge1J_J1PFcands_dRle0p1"] = ak.flatten(genWeights["ge1J_J1PFcands"][masks["deltaR_J1_J1PFcands_le0p1"]])
        genWeights["ge2J"] = genWeights["noCut"][masks["ge2J"]]
        genWeights["ge1j"] = genWeights["noCut"][masks["ge1j"]]

        ## Running some sanity checks
        check_dict_keys(genWeights, { k: "" for k in self.genWeightsInfo.values() }, "gen weight", "gen weight initialized in constructor")

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

