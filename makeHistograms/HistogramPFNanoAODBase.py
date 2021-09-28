from coffea import processor
import awkward as ak
import sys

sys.path.append("../utilities/")
import utilities as utl
import coffeaUtilities as cfutl
import PtEtaPhiMLorentzVectorUtilities as vecutl
import initializeAkArrayHelper as init_helper

from HistogramDefault import HistogramDefault

      

class HistogramPFNanoAODBase(HistogramDefault):
    """Coffea processor for accumulating histograms of selected variables.

    Some variables are directly read from the ROOT file while other are
    computed on the fly.

    This must implement the following methods:
        * __init__
        * initialize_ak_arrays
        * initialize_gen_weights

    __init__ defines the following attributes:
        * self.njet_max
        * self.jet_types
        * self.variables
        * self.gen_weights_info
        * self.cuts
        * self._accumulator

    initialize_ak_arrays returns the following objects:
        * var_arrays (dict): arrays for all variables
        * masks (dict): masks used for the necessary cuts to define some variables

    initialize_gen_weights returns the following objects:
        * gen_weights(dict): generator weights for all variables
    """


    def __init__(self, binning_info):
        """Define the variables to histogram, their binning and their coffea histogram object.

        Args:
            binning_info (dict): Form of the dictionary:
		{
		    "noregex": {
			"ak8Jet1_ak8Jet2_mass"  :  [500  ,  0   , 5000 ],
		    },
		    "regex": {
			"ak[48]Jet_n"           :  [20   ,  0   , 20   ],
		    }
		}

        Returns:
            None
        """

        self.binning_info = binning_info

        self.njet_max = 4
        self.jet_types = ["ak4", "ak8"]

        ## Iterable over jets
        it1 = [ {"jet_type": "ak4"}, {"jet_type": "ak8"} ]
        it11 = [ {"jet_type": "ak8"} ]

        ## Iterable over jets, cuts and jet number
        it2 = []
        it3 = []
        for njet in range(1, self.njet_max+1):
            for ijet in range(1, njet+1):
                for jet_type in self.jet_types:
                    it2.append({"jet_type": jet_type, "cut": "ge"+str(njet)+jet_type, "n": str(ijet)})
                for jet_type in ["ak8"]:
                    it3.append({"jet_type": jet_type, "cut": "ge"+str(njet)+jet_type, "n": str(ijet)})

        ## Iterable over cuts
        it4 = [ {"cut": ""} ]
        for njet in range(1, self.njet_max+1):
            it4.append({"cut": "_ge"+str(njet)+jet_type})

        ## Iterable over cuts and jets
        it5 = []
        for njet in range(0, self.njet_max+1):
            for jet_type in self.jet_types:
                if njet == 0:
                    it5.append({"jet_type": jet_type, "cut": ""})
                else:
                    it5.append({"jet_type": jet_type, "cut": "_ge"+str(njet)+jet_type})

        ## Iterable over cuts and jets
        it51 = []
        for njet in range(1, self.njet_max+1):
            for jet_type in self.jet_types:
                    it51.append({"jet_type": jet_type, "cut": "_ge"+str(njet)+jet_type})

 
        ## Iterable over jet_type, cuts and jet pair numbers
        it6 = []
        for jet_type in self.jet_types:
           for njet in range(2, self.njet_max+1):
               for ijet1, ijet2 in init_helper.make_pairs(njet):
                   it6.append({"jet_type": jet_type, "cut": "ge"+str(njet)+jet_type, "n1": str(ijet1+1), "n2": str(ijet2+1)})

        ## Iterable over jet and cuts 
        it7 = []
        for jet_type in self.jet_types:
           for njet in range(2, self.njet_max+1):
               it7.append({"jet_type": jet_type, "cut": "ge"+str(njet)+jet_type})

        
        variables_description = [
            ( "{jet_type}Jet_n"           , it1 ),
            ( "{jet_type}Jet_pt"          , it1 ),
            ( "{jet_type}Jet_eta"         , it1 ),
            ( "{jet_type}Jet_mass"        , it1 ),
            ( "{jet_type}Jet_msoftdrop"   , it11 ),
            ( "{jet_type}Jet_tau1"        , it11 ),
            ( "{jet_type}Jet_tau2"        , it11 ),
            ( "{jet_type}Jet_tau3"        , it11 ),
            ( "{jet_type}Jet_tau4"        , it11 ),
            ( "{jet_type}Jet_tau21"       , it11 ),
            ( "{jet_type}Jet_tau32"       , it11 ),
            ( "{jet_type}Jet_tau43"       , it11 ),
            ( "{jet_type}Jet_n2b1"        , it11 ),
            ( "{jet_type}Jet_n3b1"        , it11 ),

            ( "{jet_type}Jet{n}_pt_{cut}"        , it2 ),
            ( "{jet_type}Jet{n}_mass_{cut}"      , it2 ),
            ( "{jet_type}Jet{n}_msoftdrop_{cut}" , it3 ),
            ( "{jet_type}Jet{n}_tau1_{cut}"      , it3 ),
            ( "{jet_type}Jet{n}_tau2_{cut}"      , it3 ),
            ( "{jet_type}Jet{n}_tau3_{cut}"      , it3 ),
            ( "{jet_type}Jet{n}_tau4_{cut}"      , it3 ),
            ( "{jet_type}Jet{n}_tau21_{cut}"     , it3 ),
            ( "{jet_type}Jet{n}_tau32_{cut}"     , it3 ),
            ( "{jet_type}Jet{n}_tau43_{cut}"     , it3 ),
            ( "{jet_type}Jet{n}_n2b1_{cut}"      , it3 ),
            ( "{jet_type}Jet{n}_n3b1_{cut}"      , it3 ),

            ( "MET_pt{cut}"      , it5 ),
            ( "HT{jet_type}{cut}"     , it5 ),
            ( "ST{jet_type}{cut}"     , it5 ),
            ( "METrHT{jet_type}{cut}" , it51),
            ( "METrST{jet_type}{cut}" , it5 ),

            ( "deltaR_{jet_type}Jet{n1}_{jet_type}Jet{n2}_{cut}"   , it6 ),
            ( "deltaPhi_{jet_type}Jet{n1}_{jet_type}Jet{n2}_{cut}" , it6 ),
            ( "deltaEta_{jet_type}Jet{n1}_{jet_type}Jet{n2}_{cut}" , it6 ),
            ( "deltaPhi_{jet_type}Jet{n1}{jet_type}Jet{n2}_MET_{cut}" , it6 ),

            ( "deltaPhi_{jet_type}Jet{n}_MET_{cut}" , it2 ),
            ( "deltaPhiMin_{jet_type}Jet_MET_{cut}" , it7 ),

            ( "{jet_type}Jet{n1}_{jet_type}Jet{n2}_mass_{cut}" , it6 ),
            ( "{jet_type}Jet{n1}_{jet_type}Jet{n2}_pt_{cut}"   , it6 ),
        ]

        simple_variables = [
            "MET_phi",
            "MTak8_ge2ak8",
            "MTak4_ge2ak4",
            "MTwrongak8_ge2ak8",
            "MTwrongak4_ge2ak4",
            "RTak8_ge2ak8",
            "RTak4_ge2ak4",
        ]

        self.variables = self.make_variables_names(variables_description) + simple_variables


        ## In the process method, a gen_weights dict is defined, storing the gen weights ak array to use for the different variables
        #  The matching between gen_weights and variables arrays is automatic
        #  But in case of degenerate cases, the automatic matching will fail
        #  Here we define which gen weights to use for the different variables in the degenerate case
        self.gen_weights_info = {}


        ## List of cuts necessary for computing some variables (i.e. need at least 2 jets for computing delta R between leading 2 jets!)
        self.cuts = ["noCut"] 
        for jet_type in self.jet_types:
            for njet in range(1, self.njet_max+1):
                self.cuts.append("ge"+str(njet)+jet_type)
        

        ## Get binning of for all variables to histogram
        self.make_binning()

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
        """Read or compute all variables to histogram.

        Read or compute all variables to histogram.
        Make masks for the cuts necessary to define some variables.

        Args:
            events (ak.Array): the Events TTree opened with uproot.
        
        Returns:
            (tuple) tuple containing:
                var_arrays (dict[str, ak.Array])
                masks (dict[str, ak.Array])
        """

        ## Define dict storing the different arrays
        var_arrays = {}
        jagged_var_arrays = {}

        ## Define dict storing the different masks used
        masks = {}

        ## Compute variables to histogram / defined in constructor

        # Basic jet variables
        jet_variables = [
            "pt", "eta", "mass", "msoftdrop",
            "tau1", "tau2", "tau3", "tau4", "n2b1", "n3b1"]


        # Event-level variables not requiring jet info
        jagged_var_arrays["MET_4vector"] = vecutl.make_PtEtaPhiMLorentzVector(
            cfutl.get_from_events(events, "MET_pt"),
            ak.zeros_like(cfutl.get_from_events(events, "MET_pt")),
            cfutl.get_from_events(events, "MET_phi"),
            ak.zeros_like(cfutl.get_from_events(events, "MET_pt")),
        )
        var_arrays["MET_pt"] = jagged_var_arrays["MET_4vector"].pt
        var_arrays["MET_phi"] = jagged_var_arrays["MET_4vector"].phi


        # Looping over all jet types
        for jet_type in self.jet_types:
            # This could be refined fer Delphes etc...
            jet_collection = "FatJet" if jet_type == "ak8" else "Jet"


            # Making jet 4-vectors
            jagged_var_arrays[jet_type+"Jet_4vector"] = vecutl.make_PtEtaPhiMLorentzVector(
                cfutl.get_from_events(events, jet_collection+"_pt"),
                cfutl.get_from_events(events, jet_collection+"_eta"),
                cfutl.get_from_events(events, jet_collection+"_phi"),
                cfutl.get_from_events(events, jet_collection+"_mass"),
            )

            # Making jet constituents 4-vectors
            jagged_var_arrays[jet_type+"JetPFCands4vector"] = vecutl.make_PtEtaPhiMLorentzVector(
                cfutl.get_from_events(events, "PFCands_pt"),
                cfutl.get_from_events(events, "PFCands_eta"),
                cfutl.get_from_events(events, "PFCands_phi"),
                cfutl.get_from_events(events, "PFCands_mass"),
            )


            # Reading jet "basic" variables for all jets in each event (flatten the jagged array)
            init_helper.read_basic_variables(events, jet_type, jet_collection, jet_variables, jagged_var_arrays, var_arrays)
           
            # Computing ratios of nsubjettiness
            for tauA, tauB, tauRatio in (("tau2", "tau1", "tau21"), ("tau3", "tau2", "tau32"), ("tau4", "tau3", "tau43")):
                tauRatioBranchName = jet_collection + "_" + tauRatio
                tauAVariableName = jet_type + "Jet_" + tauA
                tauBVariableName = jet_type + "Jet_" + tauB
                tauRatioVariableName = jet_type + "Jet_" + tauRatio
                init_helper.compute_nsubjettiness_ratio(events, tauRatioBranchName, tauAVariableName, tauBVariableName, tauRatioVariableName, jagged_var_arrays, var_arrays)
                if var_arrays[tauRatioVariableName] is not None: jet_variables.append(tauRatio)


            init_helper.make_njet_masks(events, jet_type, jet_collection, self.njet_max, masks, jet_variables[0])
            for njet in range(1, self.njet_max+1):
                init_helper.make_masked_jagged_array_shorthands(jet_type, njet, masks, jagged_var_arrays)

            # Making array of the above quantities for leading, subleading ... jets for event with more than 1, 2 ... jets
            for njet in range(1, self.njet_max+1):
                init_helper.compute_variables_per_jet(jet_variables, jet_type, njet, jagged_var_arrays, var_arrays, masks)

            # Making some quantites involving jet and jet constituents 4-vectors
            for njet in range(1, self.njet_max+1):

                ge_njet = "ge" + str(njet) + jet_type   # shorthand
                mask_ge_njet = masks[ge_njet]        # shorthand

                # delta R, phi eta between any pair of jets
                # mass and pt of any pair of jets
                if njet >= 2:
                    for ijet1, ijet2 in init_helper.make_pairs(njet):
                        # shorthands
                        sijet1 = str(ijet1+1)
                        sijet2 = str(ijet2+1)
                        j1 = jagged_var_arrays[jet_type+"Jet_4vector_"+ge_njet][:, ijet1]
                        j2 = jagged_var_arrays[jet_type+"Jet_4vector_"+ge_njet][:, ijet2]
                        met = jagged_var_arrays["MET_4vector_"+ge_njet]
                        suffix = jet_type+"Jet"+sijet1+"_"+jet_type+"Jet"+sijet2+"_"+ge_njet
                        
                        var_arrays["deltaR_"+suffix] = vecutl.delta_r(j1, j2)
                        var_arrays["deltaPhi_"+suffix] = vecutl.abs_delta_phi(j1, j2)
                        var_arrays["deltaEta_"+suffix] = vecutl.abs_delta_eta(j1, j2)

                        suffix = jet_type+"Jet"+sijet1+jet_type+"Jet"+sijet2+"_MET_"+ge_njet
                        var_arrays["deltaPhi_"+suffix] = vecutl.abs_delta_phi(j1+j2, met)

                        prefix = jet_type+"Jet"+sijet1+"_"+jet_type+"Jet"+sijet2
                        var_arrays[prefix + "_mass_"+ge_njet] = vecutl.mass(j1, j2)
                        var_arrays[prefix + "_pt_"+ge_njet] = vecutl.pt(j1, j2)


                # delta phi between any jet and MET
                for ijet in range(njet):
                    sijet = str(ijet+1)
                    j = jagged_var_arrays[jet_type+"Jet_4vector_"+ge_njet][:, ijet]
                    met = jagged_var_arrays["MET_4vector_"+ge_njet]
                    var_arrays["deltaPhi_"+jet_type+"Jet"+sijet+"_MET_"+ge_njet] = vecutl.abs_delta_phi(j, met)
                
                # delta phi min between any jet and MET
                if njet >=2:
                    listOfAkArrays = [ var_arrays["deltaPhi_"+jet_type+"Jet"+str(ijet+1)+"_MET_"+ge_njet] for ijet in range(njet) ]
                    var_arrays["deltaPhiMin_"+jet_type+"Jet_MET_"+ge_njet] = ak.min(ak.Array(listOfAkArrays), axis=0)



            # Some other event-level variables requiring jet info
            # MET for different cuts
            for njet in range(1, self.njet_max+1):
                ge_njet = "ge" + str(njet) + jet_type   # shorthand
                mask_ge_njet = masks[ge_njet]           # shorthand
                var_arrays["MET_pt_"+ge_njet] = var_arrays["MET_pt"][mask_ge_njet]
                var_arrays["MET_phi_"+ge_njet] = var_arrays["MET_phi"][mask_ge_njet]

            # HT, ST
            var_arrays["HT"+jet_type] = ak.sum(jagged_var_arrays[jet_type+"Jet_pt"], axis=1)
            var_arrays["ST"+jet_type] = var_arrays["HT"+jet_type] + var_arrays["MET_pt"]

            # MT
            if self.njet_max >= 2:
                ge_njet = "ge2" + jet_type           # shorthand
                jj_4vector = jagged_var_arrays[jet_type+"Jet_4vector_"+ge_njet][:, 0] + jagged_var_arrays[jet_type+"Jet_4vector_"+ge_njet][:, 1]
                met_4vector = jagged_var_arrays["MET_4vector_"+ge_njet]
                var_arrays["MT"+jet_type+"_"+ge_njet] = vecutl.mt(jj_4vector, met_4vector)
                var_arrays["MTwrong"+jet_type+"_"+ge_njet] = vecutl.mt_wrong(jj_4vector, met_4vector)


            # RT
            var_arrays["RT"+jet_type+"_"+ge_njet] = var_arrays["MET_pt_"+ge_njet] / var_arrays["MT"+jet_type+"_"+ge_njet]
         
            # Ratios of MET with HT, ST
            var_arrays["METrST"+jet_type] = var_arrays["MET_pt"] / var_arrays["ST"+jet_type]
            for njet in range(1, self.njet_max+1):
                ge_njet = "ge" + str(njet) + jet_type   # shorthand
                mask_ge_njet = masks[ge_njet]           # shorthand
                var_arrays["HT"+jet_type+"_"+ge_njet] = var_arrays["HT"+jet_type][mask_ge_njet]
                var_arrays["ST"+jet_type+"_"+ge_njet] = var_arrays["ST"+jet_type][mask_ge_njet]
                var_arrays["METrHT"+jet_type+"_"+ge_njet] = var_arrays["MET_pt"][mask_ge_njet] / var_arrays["HT"+jet_type][mask_ge_njet]
                var_arrays["METrST"+jet_type+"_"+ge_njet] = var_arrays["MET_pt"][mask_ge_njet] / var_arrays["ST"+jet_type][mask_ge_njet]


        return var_arrays, masks


    def initialize_gen_weights(self, events, masks):
        """Make events gen weights for the different variables.

        Args:
            events (ak.Array): the Events TTree open with uproot.
            masks (dict[str, ak.Array])

        Returns:
            dict[str, ak.Array]
        """

        gen_weights = {}
        gen_weights["noCut"] = events["genWeight"]
        gen_weights["noCut_ak8bc"] = ak.flatten(ak.broadcast_arrays(gen_weights["noCut"], events["FatJet_pt"])[0])
        gen_weights["noCut_ak4bc"] = ak.flatten(ak.broadcast_arrays(gen_weights["noCut"], events["Jet_pt"])[0])

        for jet_type in self.jet_types:
            for njet in range(1, self.njet_max+1):
                ge_njet = "ge" + str(njet) + jet_type   # shorthand
                mask_ge_njet = masks[ge_njet]           # shorthand
                gen_weights[ge_njet] = gen_weights["noCut"][mask_ge_njet]

        return gen_weights


