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


def divide_ak_arrays(akArray1, akArray2, divisionByZeroValue=1.):
    """
    Input: 2 ak arrays with floatting point values having the same jagged structure.
    Return a jagged list having same jagged structure corresponding to
           akArray1 / akArray2
           where the division by 0 is replaced by divisionByZeroValue
           e.g. akArray1 = [ [0, 3], [5], [2] ] akArray2 = [ [3, 3], [0], [1] ]
                division = [ [0, 1], [1], [0.5] ]
    """

    isNotZero = (akArray2!=0.)
    if (not ak.all(isNotZero)):
        print("The following warning about true_divide can be safely ignored.")

    rawDivision = akArray1/akArray2
    division = ak.where(isNotZero, rawDivision, divisionByZeroValue*ak.ones_like(akArray1))

    # This implementation seems slower:
    #division = ak.Array([ [ x1/x2 if x2 != 0. else divisionByZeroValue for x1, x2 in zip(y1, y2) ] for y1, y2 in zip(akArray1, akArray2) ])

    return division



def get_from_events(events, branchName):
    """Return the branch from the events TTree if it exists, else return None."""

    if branchName in events.fields:
        return events[branchName]
    else:
        return None


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



def fill_histogram(output, variable, axis, weight):
    """Fill histogram if the variable has been defined."""

    type1 = str(ak.type(axis))
    type2 = str(ak.type(weight))
    type1 = utl.list2str(type1.split(" * ")[:-1], " * ")
    type2 = utl.list2str(type2.split(" * ")[:-1], " * ")
    if variable in output.keys():
        if type1 == type2:
            output[variable].fill(axis=axis, weight=weight)
        else:
            print("ERROR: Array and weight for variable %s have different sizes:" %variable)
            print("array  size: %s" %(type1))
            print("weight size: %s" %(type2))
            sys.exit()
    else:
        pass
        # A warning has been sent earlier when defining coffea histograms


def make_variables_names(variablesDescriptions):
    """
    Make variables based on the following description from a collection of templates and list of
    params to be replaced by a given value:
    e.g. "variablesDescriptions = [ ( "{jet}Jet_n",  [ {"jet": "ak4"}, {"jet": "ak8"} ] ),
                                      "{jet}Jet_pt", [ {"jet": "ak4"}, {"jet": "ak8"} ] ),
                                  ]

    """

    variables = []

    for entry in variablesDescriptions:
        variableTemplate = entry[0]
        listOfParamsDict = entry[1]
        for paramsDict in listOfParamsDict:
            var = variableTemplate
            for k, v in paramsDict.items():
                var = var.replace("{"+k+"}", v)
            variables.append(var)

    return variables



def get_binning(variable, binning):
    """
    Return binning for variable in 1st argument providing the binning information in 2nd argument.
    binning = {
        "noregex": {
            "ak8J1_ak8J2_mass"  :  [500  ,  0   , 5000 ],
        },
        "regex": {
            "n(Ak8)?Jet"        :  [20   ,  0   , 20   ],
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
            print("regexes matched:")
            for idx in indices: print("\t%s" %(regexes[idx]))
            sys.exit()
        else: 
            binning_ = binning["regex"][regexes[indices[0]]]
    else:
        binning_ = binning["noregex"][variable]

    return binning_



def check_dict_keys(dict1, dict2, dict1Name, dict2Name, errorKeysDict1NotInKeysDict2=True, errorKeysDict2NotInKeysDict1=True):
    """Check if dict1 and dict2 have the same keys. Interrupt program if not the case."""

    if isinstance(dict1, list):
        dict1 = { k: "" for k in dict1 }
    if isinstance(dict2, list):
        dict1 = { k: "" for k in dict2 }

    if errorKeysDict1NotInKeysDict2:
        messageLevel1 = "ERROR"
    else:
        messageLevel1 = "WARNING"
    if errorKeysDict2NotInKeysDict1:
        messageLevel2 = "ERROR"
    else:
        messageLevel2 = "WARNING"

    error = False
    dict1NotNoneKeys = sorted([ k for k in dict1.keys() if dict1[k] is not None ])
    dict2NotNoneKeys = sorted([ k for k in dict2.keys() if dict2[k] is not None ])
    if dict1NotNoneKeys != dict2NotNoneKeys:
        dict1KeysNotInDict2 = [ k for k in dict1NotNoneKeys if k not in dict2NotNoneKeys ]
        dict2KeysNotInDict1 = [ k for k in dict2NotNoneKeys if k not in dict1NotNoneKeys ]
        if len(dict1KeysNotInDict2) > 0:
            if errorKeysDict2NotInKeysDict1: error = True
            print("\n%s: Some %s have been defined but do not have a defined %s" %(messageLevel2, dict1Name, dict2Name))
            print("Missing %s:" %dict2Name)
            print(dict1KeysNotInDict2)
        if len(dict2KeysNotInDict1) > 0:
            if errorKeysDict1NotInKeysDict2: error = True
            print("\n%s: Some %s have a defined %s but have not been defined" %(messageLevel1, dict1Name, dict2Name))
            print("Missing %s:" %dict1Name)
            print(dict2KeysNotInDict1)
        if error:
            sys.exit()
        else:
            return False

    return True
 

def make_pairs(n):
    """Return pairs of distinct integers from {0, 1, ...n}."""

    pairs = []
    for i1 in range(n-1):
        for i2 in range(i1+1, n):
            pairs.append([i1, i2])
    return pairs


class Histogram1(processor.ProcessorABC):
    """
    Make histograms of selected variables in input ROOT file.
    If variable does not exist in ROOT file, they are computed on the fly (e.g. HT, MET/HT, ...).
    """


    def __init__(self, binning, fileType, VERBOSE=1):
        """Define the variables to histogram, their binning and their coffea histogram object."""

        self.fileType = fileType

        self.nJetMax = 4

        ## Variables to histogram and their x label
        # Variables without cuts
        self.jets = ["ak4", "ak8"]

        ## Iterable over jets
        it1 = [ {"jet": "ak4"}, {"jet": "ak8"} ]
        it11 = [ {"jet": "ak8"} ]

        ## Iterable over jets, cuts and jet number
        it2 = []
        it3 = []
        for njet in range(1, self.nJetMax+1):
            for ijet in range(1, njet+1):
                for jet in self.jets:
                    it2.append({"jet": jet, "cut": "ge"+str(njet)+jet, "n": str(ijet)})
                for jet in ["ak8"]:
                    it3.append({"jet": jet, "cut": "ge"+str(njet)+jet, "n": str(ijet)})

        ## Iterable over cuts
        it4 = [ {"cut": ""} ]
        for njet in range(1, self.nJetMax+1):
            it4.append({"cut": "_ge"+str(njet)+jet})

        ## Iterable over cuts and jets
        it5 = []
        for njet in range(0, self.nJetMax+1):
            for jet in self.jets:
                if njet == 0:
                    it5.append({"jet": jet, "cut": ""})
                else:
                    it5.append({"jet": jet, "cut": "_ge"+str(njet)+jet})

        ## Iterable over cuts and jets
        it51 = []
        for njet in range(1, self.nJetMax+1):
            for jet in self.jets:
                    it51.append({"jet": jet, "cut": "_ge"+str(njet)+jet})

 
        ## Iterable over jet, cuts and jet pair numbers
        it6 = []
        for jet in self.jets:
           for njet in range(2, self.nJetMax+1):
               for ijet1, ijet2 in make_pairs(njet):
                   it6.append({"jet": jet, "cut": "ge"+str(njet)+jet, "n1": str(ijet1+1), "n2": str(ijet2+1)})

        ## Iterable over jet and cuts 
        it7 = []
        for jet in self.jets:
           for njet in range(2, self.nJetMax+1):
               it7.append({"jet": jet, "cut": "ge"+str(njet)+jet})

        
        variablesDescription = [
            ( "{jet}Jet_n"           , it1 ),
            ( "{jet}Jet_pt"          , it1 ),
            ( "{jet}Jet_eta"         , it1 ),
            ( "{jet}Jet_mass"        , it1 ),
            ( "{jet}Jet_msoftdrop"   , it11 ),
            ( "{jet}Jet_tau1"        , it11 ),
            ( "{jet}Jet_tau2"        , it11 ),
            ( "{jet}Jet_tau3"        , it11 ),
            ( "{jet}Jet_tau21"       , it11 ),
            ( "{jet}Jet_tau32"       , it11 ),

            ( "{jet}Jet{n}_pt_{cut}"        , it2 ),
            ( "{jet}Jet{n}_mass_{cut}"      , it2 ),
            ( "{jet}Jet{n}_msoftdrop_{cut}" , it3 ),
            ( "{jet}Jet{n}_tau1_{cut}"      , it3 ),
            ( "{jet}Jet{n}_tau2_{cut}"      , it3 ),
            ( "{jet}Jet{n}_tau3_{cut}"      , it3 ),
            ( "{jet}Jet{n}_tau21_{cut}"     , it3 ),
            ( "{jet}Jet{n}_tau32_{cut}"     , it3 ),

            ( "MET_pt{cut}"      , it5 ),
            ( "HT{jet}{cut}"     , it5 ),
            ( "ST{jet}{cut}"     , it5 ),
            ( "METrHT{jet}{cut}" , it51),
            ( "METrST{jet}{cut}" , it5 ),

            ( "deltaR_{jet}Jet{n1}_{jet}Jet{n2}_{cut}"   , it6 ),
            ( "deltaPhi_{jet}Jet{n1}_{jet}Jet{n2}_{cut}" , it6 ),
            ( "deltaEta_{jet}Jet{n1}_{jet}Jet{n2}_{cut}" , it6 ),

            ( "deltaPhi_{jet}Jet{n}_MET_{cut}" , it2 ),
            ( "deltaPhiMin_{jet}Jet_MET_{cut}" , it7 ),

            ( "{jet}Jet{n1}_{jet}Jet{n2}_mass_{cut}" , it6 ),
            ( "{jet}Jet{n1}_{jet}Jet{n2}_pt_{cut}"   , it6 ),


            ( "deltaR_{jet}Jet{n}_{jet}Jet{n}PFCands_{cut}" , it2 ),
            ( "sum_{jet}Jet{n}PFCandsPt_deltaR0.0To0.1_r_{jet}Jet{n}Pt_{cut}" , it2 ),
            ( "sum_{jet}Jet{n}PFCandsPt_deltaR0.1To0.2_r_{jet}Jet{n}Pt_{cut}" , it2 ),
            ( "sum_{jet}Jet{n}PFCandsPt_deltaR0.2To0.3_r_{jet}Jet{n}Pt_{cut}" , it2 ),
            ( "sum_{jet}Jet{n}PFCandsPt_deltaR0.3To0.4_r_{jet}Jet{n}Pt_{cut}" , it2 ),

        ]

        simpleVariables = [
            "MET_phi",
            "MTak8_ge2ak8",
            "MTak4_ge2ak4",
            "RTak8_ge2ak8",
            "RTak4_ge2ak4",
        ]

        self.variables = make_variables_names(variablesDescription) + simpleVariables


        ## In the process method, a genWeights dict is defined, storing the gen weights ak array to use for the different variables
        #  The matching between genWeights and variables arrays is automatic
        #  But in case of degenerate cases, the automatic matching will fail
        #  Here we define which gen weights to use for the different variables in the degenerate case
        self.genWeightsInfo = {}


        ## List of cuts necessary for computing some variables (i.e. need at least 2 jets for computing delta R between leading 2 jets!)
        #  Note that the following convention is used:
        #  if cut is in genWeights[variable] then self.sumGenWeightsInfo[variable] = cut
        self.cuts = ["noCut"] 
        for jet in self.jets:
            for njet in range(1, self.nJetMax+1):
                self.cuts.append("ge"+str(njet)+jet)
        

        ## Make dict filled with coffea histograms to be used in the accumulator
        histograms = {}

        # Loop over all variables to histogram
        for variable in self.variables:
            # x axis label
            label = variable
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


        ## Define dict storing the different arrays
        varArrays = {}
        jaggedVarArrays = {}

        ## Define dict storing the different masks used
        masks = {}

        ## Compute variables to histogram / defined in constructor

        # Basic jet variables
        jetVariables = ["pt", "eta", "mass", "msoftdrop", "tau1", "tau2", "tau3"]


        # Event-level variables not requiring jet info
        jaggedVarArrays["MET_4vector"] = make_PtEtaPhiMLorentzVector(
            get_from_events(events, "MET_pt"),
            ak.zeros_like(get_from_events(events, "MET_pt")),
            get_from_events(events, "MET_phi"),
            ak.zeros_like(get_from_events(events, "MET_pt")),
        )
        varArrays["MET_pt"] = jaggedVarArrays["MET_4vector"].pt
        varArrays["MET_phi"] = jaggedVarArrays["MET_4vector"].phi


        # Looping over all jet types
        for jet in self.jets:
            # This could be refined fer Delphes etc...
            if jet == "ak8":
                PFnanoJet = "FatJet"
            else:
                PFnanoJet = "Jet"


            # Making jet 4-vectors
            jaggedVarArrays[jet+"Jet_4vector"] = make_PtEtaPhiMLorentzVector(
                get_from_events(events, PFnanoJet+"_pt"),
                get_from_events(events, PFnanoJet+"_eta"),
                get_from_events(events, PFnanoJet+"_phi"),
                get_from_events(events, PFnanoJet+"_mass"),
            )

            # Making jet constituents 4-vectors
            if self.fileType == "PFnano102X":
                if jet == "ak8": prefix = "Fat"
                else: prefix = ""
            elif self.fileType == "PFnano106X":
                prefix = ""
            # the else case cannot happen, it has already been tackled
            jaggedVarArrays[jet+"JetPFCands4vector"] = make_PtEtaPhiMLorentzVector(
                get_from_events(events, prefix+"JetPFCands_pt"),
                get_from_events(events, prefix+"JetPFCands_eta"),
                get_from_events(events, prefix+"JetPFCands_phi"),
                get_from_events(events, prefix+"JetPFCands_mass"),
            )


            # Reading jet "basic" variables for all jets in each event (flatten the jagged array)
            varArrays[jet+"Jet_n"] = get_from_events(events, "n"+PFnanoJet)
            for var in jetVariables:
                jaggedVarArrays[jet+"Jet_"+var] = get_from_events(events, PFnanoJet+"_"+var)
                if jaggedVarArrays[jet+"Jet_"+var] is not None:
                    varArrays[jet+"Jet_"+var] = ak.flatten(jaggedVarArrays[jet+"Jet_"+var])
                else:
                    varArrays[jet+"Jet_"+var] = None
            
            # Computing simple composed variables (e.g. tau21 = tau2/tau1)
            if get_from_events(events, PFnanoJet+"_tau21") is not None:
                jaggedVarArrays[jet+"Jet_tau21"] = get_from_events(events, PFnanoJet+"_tau21")
                varArrays[PFnanoJet+"_tau21"] = ak.flatten(jaggedVarArrays[jet+"Jet_tau21"])
                jetVariables.append("tau21")
            elif (varArrays[jet+"Jet_tau1"] is not None) and (varArrays[jet+"Jet_tau2"] is not None):
                jaggedVarArrays[jet+"Jet_tau21"] = divide_ak_arrays(jaggedVarArrays[jet+"Jet_tau2"], jaggedVarArrays[jet+"Jet_tau1"])
                varArrays[jet+"Jet_tau21"] = ak.flatten(jaggedVarArrays[jet+"Jet_tau21"])
                jetVariables.append("tau21")
            else:
                jaggedVarArrays[jet+"Jet_tau21"] = None
                varArrays[jet+"Jet_tau21"] = None

            if get_from_events(events, PFnanoJet+"_tau32") is not None:
                jaggedVarArrays[jet+"Jet_tau32"] = get_from_events(events, PFnanoJet+"_tau32")
                varArrays[PFnanoJet+"_tau32"] = ak.flatten(jaggedVarArrays[jet+"Jet_tau32"])
                jetVariables.append("tau32")
            elif (varArrays[jet+"Jet_tau2"] is not None) and (varArrays[jet+"Jet_tau3"] is not None):
                jaggedVarArrays[jet+"Jet_tau32"] = divide_ak_arrays(jaggedVarArrays[jet+"Jet_tau3"], jaggedVarArrays[jet+"Jet_tau2"])
                varArrays[jet+"Jet_tau32"] = ak.flatten(jaggedVarArrays[jet+"Jet_tau32"])
                jetVariables.append("tau32")
            else:
                jaggedVarArrays[jet+"Jet_tau32"] = None
                varArrays[jet+"Jet_tau32"] = None


            # Making array of the above quantities for leading, subleading ... jets for event with more than 1, 2 ... jets
            for njet in range(1, self.nJetMax+1):
                geNJet  =  "ge" + str(njet) + jet   # shorthand
                masks[geNJet] = (get_from_events(events, "n"+PFnanoJet) > njet-1)
                masks[geNJet+"_bc"] = ak.flatten(ak.broadcast_arrays(masks[geNJet], get_from_events(events, PFnanoJet+"_pt"))[0])

                for ijet in range(1, njet+1):
                    sijet = str(ijet)
                    for var in jetVariables:
                        if jaggedVarArrays[jet+"Jet_"+var] is not None:
                            varArrays[jet+"Jet"+sijet+"_"+var+"_"+geNJet] = jaggedVarArrays[jet+"Jet_"+var][masks[geNJet]][:, ijet-1]
                        else:
                            varArrays[jet+"Jet"+sijet+"_"+var+"_"+geNJet] = None


            # Making some quantites involving jet and jet constituents 4-vectors
            for njet in range(1, self.nJetMax+1):

                geNJet = "ge" + str(njet) + jet   # shorthand
                maskGeNJet = masks[geNJet]        # shorthand

                # MET 4-vector for events with at least njet jets
                jaggedVarArrays["MET_4vector_"+geNJet] = jaggedVarArrays["MET_4vector"][maskGeNJet]

                # jet 4-vector for events with at least njet jets
                jaggedVarArrays[jet+"Jet_4vector_"+geNJet] = jaggedVarArrays[jet+"Jet_4vector"][maskGeNJet]

                for ijet in range(1, njet+1):
                    sijet = str(ijet)
                    jaggedVarArrays[jet+"Jet"+sijet+"_4vector_"+geNJet] = jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet-1]

                    if self.fileType == "PFnano102X":
                        if jet == "ak8": prefix = "Fat"
                        else: prefix = ""
                        masks[jet+"JetIdx"+sijet+"_"+geNJet] = (get_from_events(events, prefix+"JetPFCands_jetIdx")[maskGeNJet] == ijet-1)
                        mask = masks[jet+"JetIdx"+sijet+"_"+geNJet]  # shorthand
                        jaggedVarArrays[jet+"Jet"+sijet+"PFCands4vector_"+geNJet] = jaggedVarArrays[jet+"JetPFCands4vector"][maskGeNJet][mask]

                    elif self.fileType == "PFnano106X":
                        masks[jet+"JetIdx"+sijet+"_"+geNJet] = (get_from_events(events, "JetPFCands"+jet.upper()+"_jetIdx")[maskGeNJet] == ijet-1)
                        masks[jet+"CandIdx"+sijet+"_"+geNJet] = events["JetPFCands"+jet.upper()+"_candIdx"][maskGeNJet][masks[jet+"JetIdx"+sijet+"_"+geNJet]]
                        mask = masks[jet+"CandIdx"+sijet+"_"+geNJet]  # shorthand
                        jaggedVarArrays[jet+"Jet"+sijet+"PFCands4vector_"+geNJet] = jaggedVarArrays[jet+"JetPFCands4vector"][maskGeNJet][mask]



                    # deltaR between constituents and jet axis
                    jaggedVarArrays[jet+"Jet"+sijet+"_4vector_"+geNJet+"_bc"] = ak.broadcast_arrays(jaggedVarArrays[jet+"Jet"+sijet+"_4vector_"+geNJet], jaggedVarArrays[jet+"Jet"+sijet+"PFCands4vector_"+geNJet])[0]
                    jaggedVarArrays["deltaR_"+jet+"Jet"+sijet+"_"+jet+"Jet"+sijet+"PFCands_"+geNJet] = jaggedVarArrays[jet+"Jet"+sijet+"PFCands4vector_"+geNJet].delta_r(jaggedVarArrays[jet+"Jet"+sijet+"_4vector_"+geNJet+"_bc"])
                    varArrays["deltaR_"+jet+"Jet"+sijet+"_"+jet+"Jet"+sijet+"PFCands_"+geNJet] = ak.flatten(jaggedVarArrays["deltaR_"+jet+"Jet"+sijet+"_"+jet+"Jet"+sijet+"PFCands_"+geNJet])


                    # Sum jet constituents pT / jet pT   for constituents having x.x <= deltaR < y.y wrt jet axis
                    for dRCut in ((0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4)):
                        deltaRCut = str(dRCut[0]) + "To" + str(dRCut[1])
                        masks["deltaR_"+deltaRCut] = (jaggedVarArrays["deltaR_"+jet+"Jet"+sijet+"_"+jet+"Jet"+sijet+"PFCands_"+geNJet] >= dRCut[0]) & (jaggedVarArrays["deltaR_"+jet+"Jet"+sijet+"_"+jet+"Jet"+sijet+"PFCands_"+geNJet] < dRCut[1])
                        varArrays["sum_"+jet+"Jet"+sijet+"PFCandsPt_deltaR"+deltaRCut+"_r_"+jet+"Jet"+sijet+"Pt_"+geNJet] = \
                            ak.sum(jaggedVarArrays[jet+"Jet"+sijet+"PFCands4vector_"+geNJet].pt[masks["deltaR_"+deltaRCut]], axis=1) / jaggedVarArrays[jet+"Jet"+sijet+"_4vector_"+geNJet].pt

                # delta R, phi eta between any pair of jets
                # mass and pt of any pair of jets
                if njet >= 2:
                    for ijet1, ijet2 in make_pairs(njet):
                        sijet1 = str(ijet1+1)
                        sijet2 = str(ijet2+1)
                        varArrays["deltaR_"+jet+"Jet"+sijet1+"_"+jet+"Jet"+sijet2+"_"+geNJet] = jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet1].delta_r(jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet2])
                        varArrays["deltaPhi_"+jet+"Jet"+sijet1+"_"+jet+"Jet"+sijet2+"_"+geNJet] = abs(jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet1].delta_phi(jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet2]))
                        varArrays["deltaEta_"+jet+"Jet"+sijet1+"_"+jet+"Jet"+sijet2+"_"+geNJet] = abs(jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet1].eta - jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet2].eta)

                        varArrays[jet+"Jet"+sijet1+"_"+jet+"Jet"+sijet2+"_mass_"+geNJet] = (jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet1] + jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet2]).mass
                        varArrays[jet+"Jet"+sijet1+"_"+jet+"Jet"+sijet2+"_pt_"+geNJet] = (jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet1] + jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet2]).pt


                # delta phi between any jet and MET
                for ijet in range(njet):
                    sijet = str(ijet+1)
                    varArrays["deltaPhi_"+jet+"Jet"+sijet+"_MET_"+geNJet] = abs(jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, ijet].delta_phi(jaggedVarArrays["MET_4vector_"+geNJet]))
                
                # delta phi min between any jet and MET
                if njet >=2:
                    listOfAkArrays = [ varArrays["deltaPhi_"+jet+"Jet"+str(ijet+1)+"_MET_"+geNJet] for ijet in range(njet) ]
                    varArrays["deltaPhiMin_"+jet+"Jet_MET_"+geNJet] = ak.min(ak.Array(listOfAkArrays), axis=0)



            # Some other event-level variables requiring jet info
            # MET for different cuts
            for njet in range(1, self.nJetMax+1):
                geNJet = "ge" + str(njet) + jet   # shorthand
                maskGeNJet = masks[geNJet]        # shorthand
                varArrays["MET_pt_"+geNJet] = varArrays["MET_pt"][maskGeNJet]
                varArrays["MET_phi_"+geNJet] = varArrays["MET_phi"][maskGeNJet]

            # HT, ST
            varArrays["HT"+jet] = ak.sum(jaggedVarArrays[jet+"Jet_pt"], axis=1)
            varArrays["ST"+jet] = varArrays["HT"+jet] + varArrays["MET_pt"]

            # M2
            if self.nJetMax >= 2:
                geNJet = "ge2" + jet           # shorthand
                j1j2 = jet+"Jet1_"+jet+"Jet2"  # shorthand
                j1j2_4vector = jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, 0] + jaggedVarArrays[jet+"Jet_4vector_"+geNJet][:, 1]
                MET_phi = ak.zip({ "phi": varArrays["MET_phi_"+geNJet] })

                deltaPhi_j1j2_MET = abs(j1j2_4vector.delta_phi(MET_phi))
                j1j2_m2 = varArrays[j1j2+"_mass_"+geNJet]**2
                j1j2_pt = varArrays[j1j2+"_pt_"+geNJet]**2
                j1j2_pt2 = j1j2_pt**2
                MET_pt = varArrays["MET_pt_"+geNJet]
                varArrays["MT"+jet+"_"+geNJet] = np.sqrt( j1j2_m2  +  2 * ( np.sqrt(j1j2_m2 + j1j2_pt2) * MET_pt - MET_pt * j1j2_pt * np.cos(deltaPhi_j1j2_MET) ) )

            # RT
            varArrays["RT"+jet+"_"+geNJet] = varArrays["MET_pt_"+geNJet] / varArrays["MT"+jet+"_"+geNJet]
         
            # Ratios of MET with HT, ST
            varArrays["METrST"+jet] = varArrays["MET_pt"] / varArrays["ST"+jet]
            for njet in range(1, self.nJetMax+1):
                geNJet = "ge" + str(njet) + jet   # shorthand
                maskGeNJet = masks[geNJet]        # shorthand
                varArrays["HT"+jet+"_"+geNJet] = varArrays["HT"+jet][maskGeNJet]
                varArrays["ST"+jet+"_"+geNJet] = varArrays["ST"+jet][maskGeNJet]
                varArrays["METrHT"+jet+"_"+geNJet] = varArrays["MET_pt"][maskGeNJet] / varArrays["HT"+jet][maskGeNJet]
                varArrays["METrST"+jet+"_"+geNJet] = varArrays["MET_pt"][maskGeNJet] / varArrays["ST"+jet][maskGeNJet]


        return varArrays, masks


    def initialize_gen_weights(self, events, masks):
        """Make events gen weights for the different variables."""

        genWeights = {}
        genWeights["noCut"] = events["genWeight"]
        genWeights["noCut_ak8bc"] = ak.flatten(ak.broadcast_arrays(genWeights["noCut"], events["FatJet_pt"])[0])
        genWeights["noCut_ak4bc"] = ak.flatten(ak.broadcast_arrays(genWeights["noCut"], events["Jet_pt"])[0])

        ## To be automatised
        for jet in self.jets:
            for njet in range(1, self.nJetMax+1):
                geNJet = "ge" + str(njet) + jet   # shorthand
                maskGeNJet = masks[geNJet]        # shorthand
                genWeights[geNJet] = genWeights["noCut"][maskGeNJet]

                for ijet in range(njet):
                    sijet = str(ijet+1)   # shorthand

                    if self.fileType == "PFnano102X":
                        if jet == "ak8": prefix = "Fat"
                        else: prefix = ""
                        genWeights[geNJet+"_"+jet+"Jet"+sijet+"PFCandsbc"] = ak.flatten(ak.broadcast_arrays(genWeights[geNJet], events[prefix+"JetPFCands_pt"][maskGeNJet][masks[jet+"JetIdx"+sijet+"_"+geNJet]])[0])
                    elif self.fileType == "PFnano106X":
                        genWeights[geNJet+"_"+jet+"Jet"+sijet+"PFCandsbc"] = ak.flatten(ak.broadcast_arrays(genWeights[geNJet], events["JetPFCands_pt"][maskGeNJet][masks[jet+"JetIdx"+sijet+"_"+geNJet]])[0])

        return genWeights


    def fill_histograms(self, varArrays, genWeights, output):
        """Fill histograms."""

        size2gen = {}
        for key in genWeights.keys():
            size = str(ak.size(genWeights[key]))
            if size not in size2gen.keys():
                size2gen[size] = [key]
            else:
                size2gen[size].append(key)

        if len(size2gen.keys()) != len(genWeights.keys()):
            print("\nWARNING: Some gen weights have the same length, automatic matching between variable and genWeight arrays may be incorrect.")
            print("         Check the cutflow to see whether these variables have 100% efficiency.")
            print("         If not, the gen weights for the relevant variables can be defined in genWeightsInfo.")
            print("         The following genWeights have the same lengths:")
            for key, value in size2gen.items():
                if len(value) > 1:
                    print("%s with length %s" %(value, key))

        for variable in self.variables:
            # Find corresponding gen weights
            size = str(ak.size(varArrays[variable]))
            if (size not in size2gen.keys()) and (variable not in self.genWeightsInfo.keys()):
                print("\nWARNING: No gen weight for variable %s. Skipping." %variable)
            else:
                if variable in self.genWeightsInfo.keys():
                    genWeight = genWeights[self.genWeightsInfo[variable]]
                else:
                    genWeight = genWeights[size2gen[size][0]]
                fill_histogram(output, variable, varArrays[variable], genWeight)



    def process(self, events):
        """Make arrays and genWeights for all interesting variable, and fill histograms."""

        ## Define accumulator
        output = self.accumulator.identity()

        ## Make varArrays that will be used to fill histograms
        #  Also return masks used, will be useful for making gen weights
        varArrays, masks = self.initialize_ak_arrays(events)

        ## Running some sanity checks
        checkPassed = check_dict_keys(self.variables, varArrays, "variable", "ak array", errorKeysDict1NotInKeysDict2=False)
        if not checkPassed:
            print("\nWARNING: Only variable defined in constructor will be histogrammed.")

        ## Make gen weights
        genWeights = self.initialize_gen_weights(events, masks)

        ## Make sum of gen weights
        # Make the sum of gen weights for the cuts defined in contructors
        sumGenWeights = {}
        for cut in self.cuts:
            sumGenWeights[cut] = ak.sum(genWeights[cut])
        # Save cutflow information
        for cut in self.cuts:
            output["cutflow"][cut] += sumGenWeights[cut]

        ## Fill histograms
        self.fill_histograms(varArrays, genWeights, output)

        return output


    def postprocess(self, accumulator):
        return accumulator

