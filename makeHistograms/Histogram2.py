from coffea import processor
import awkward as ak
import sys

sys.path.append("../utilities/")
import utilities as utl
import awkwardArrayUtilities as akutl
import coffeaUtilities as cfutl
import PtEtaPhiMLorentzVectorUtilities as vecutl
import initializeAkArrayHelper as init_helper

from HistogramDefault import HistogramDefault


class Histogram2(HistogramDefault):
    """Make histograms of additional variables in input ROOT file."""


    def __init__(self, binning_info, file_type):
        """Define the variables to histogram, their binning and their coffea histogram object."""

        self.binning_info = binning_info
        self.file_type = file_type

        self.njet_max = 4
        self.jets = ["ak4", "ak8"]

        ## Iterable over jets
        it1 = [ {"jet": "ak4"}, {"jet": "ak8"} ]

        ## Iterable over jets, cuts and jet number
        it2 = []
        for njet in range(1, self.njet_max+1):
            for ijet in range(1, njet+1):
                for jet in self.jets:
                    it2.append({"jet": jet, "cut": "ge"+str(njet)+jet, "n": str(ijet)})

        
        variables_description = [
            ( "{jet}Jet_ptD"          , it1 ),
            ( "{jet}Jet_girth"        , it1 ),

            ( "{jet}Jet{n}_ptD_{cut}"      , it2 ),
            ( "{jet}Jet{n}_girth_{cut}"    , it2 ),

        ]

        self.variables = self.make_variables_names(variables_description)


        ## In the process method, a gen_weights dict is defined, storing the gen weights ak array to use for the different variables
        #  The matching between gen_weights and variables arrays is automatic
        #  But in case of degenerate cases, the automatic matching will fail
        #  Here we define which gen weights to use for the different variables in the degenerate case
        self.gen_weights_info = {}


        ## List of cuts necessary for computing some variables (i.e. need at least 2 jets for computing delta R between leading 2 jets!)
        self.cuts = ["noCut"] 
        for jet in self.jets:
            for njet in range(1, self.njet_max+1):
                self.cuts.append("ge"+str(njet)+jet)
        
        ## Get binning of for all variables to histogram
        self.get_binning()

        ## Make dict filled with coffea histograms to be used in the accumulator
        self.define_histograms()

        ## Define accumulator
        self._accumulator = processor.dict_accumulator({
            **self.histograms,
            **{
                # Cutflow in case cuts are applied (e.g. for delta phi between leding jets, at least 2 jets are needed!)
                'cutflow': processor.defaultdict_accumulator(float)
            }
            })


    def initialize_ak_arrays(self, events):
        """Fill a dict with the ak arrays that will be used to fill the histograms."""

 
        ## Define dict storing the different arrays
        var_arrays = {}
        jagged_var_arrays = {}

        ## Define dict storing the different masks used
        masks = {}

        ## Compute variables to histogram / defined in constructor

        # Basic jet variables
        jet_variables = ["ptD", "girth"]


        # Looping over all jet types
        for jet in self.jets:
            # This could be refined fer Delphes etc...
            jet_collection = "FatJet" if jet == "ak8" else "Jet"

            # Making jet constituents 4-vectors
            if self.file_type == "PFnano102X":
                if jet == "ak8": prefix = "Fat"
                else: prefix = ""
            elif self.file_type == "PFnano106X":
                prefix = ""
            # the else case cannot happen, it has already been tackled


            # Reading jet "basic" variables for all jets in each event (flatten the jagged array)
            init_helper.read_basic_variables(events, jet, jet_collection, jet_variables, jagged_var_arrays, var_arrays)
           
            init_helper.make_njet_masks(events, jet, jet_collection, self.njet_max, masks, jet_variables[0])

            # Making array of the above quantities for leading, subleading ... jets for event with more than 1, 2 ... jets
            for njet in range(1, self.njet_max+1):
                init_helper.compute_variables_per_jet(jet_variables, jet, njet, jagged_var_arrays, var_arrays, masks)

        return var_arrays, masks


    def initialize_gen_weights(self, events, masks):
        """Make events gen weights for the different variables."""

        gen_weights = {}
        gen_weights["noCut"] = events["genWeight"]
        gen_weights["noCut_ak8bc"] = ak.flatten(ak.broadcast_arrays(gen_weights["noCut"], events["FatJet_ptD"])[0])
        gen_weights["noCut_ak4bc"] = ak.flatten(ak.broadcast_arrays(gen_weights["noCut"], events["Jet_ptD"])[0])

        ## To be automatised
        for jet in self.jets:
            for njet in range(1, self.njet_max+1):
                ge_njet = "ge" + str(njet) + jet   # shorthand
                mask_ge_njet = masks[ge_njet]        # shorthand
                gen_weights[ge_njet] = gen_weights["noCut"][mask_ge_njet]

        return gen_weights


