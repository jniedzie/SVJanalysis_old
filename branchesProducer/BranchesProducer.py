from coffea import processor
import awkward as ak
import numpy as np
import sys

sys.path.append("../utilities/")
import coffeaUtilities as cfutl
import awkwardArrayUtilities as akutl
import PtEtaPhiMLorentzVectorUtilities as vecutl
import physicsUtilities as phutl
import variablesComputation.awkwardArray.jetVariables as jetvar


## Define some helper functions - Could be put in separate file later

def jet_algo_2_jet_collection(jet_algo):
    """Helper function converting jet algorithm name into jet collection name.

    Args:
        jet_algo (str)

    Returns:
        str
    """

    table = {
        "ak4": "Jet",
        "ak8": "FatJet",
    }

    return table[jet_algo.lower()]


def get_jets(events, jet_type):
    """Get jet collection for a given jet type.

    Args:
        events (ak.Array): the Events TTree opened with uproot.
        jet_type (str): AK4 or AK8

    Returns:
        awkward.Array: ak array with field pt, eta, phi, mass,
            that can be used as a coffea PtEtaPhiMLorentzVector.
    """

    jet_collection = jet_algo_2_jet_collection(jet_type)
    jet_branch = eval("events." + jet_collection)

    jets = ak.zip(
        {
            "pt"    : jet_branch.pt,
            "mass"  : jet_branch.mass,
            "eta"   : jet_branch.eta,
            "phi"   : jet_branch.phi,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    return jets


def get_pf_cands(events, jet_type):
    """Get jet pf candidates collection for a given jet type.

    Args:
        events (ak.Array): the Events TTree opened with uproot.
        jet_type (str): AK4 or AK8

    Returns:
        awkward.Array: ak array with field jetIdx, pt, eta, phi, mass,
            charge and pdgId that can be used as a coffea
            PtEtaPhiMLorentzVector.
    """

    # For now a quick and dirty hack to switch between 102X and 106X
    # Will be obsolete when using only 106X

    ## For 102X
    if jet_type == "AK8":
        jet_pf_cands = events.FatJetPFCands
    else:
        jet_pf_cands = events.JetPFCands
    pf_cands = ak.zip(
        {
            "jetIdx": jet_pf_cands.jetIdx,
            "pt"    : jet_pf_cands.pt,
            "mass"  : jet_pf_cands.mass,
            "eta"   : jet_pf_cands.eta,
            "phi"   : jet_pf_cands.phi,
            "charge": jet_pf_cands.charge,
            "pdgId" : jet_pf_cands.pdgId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    
    ## For 106X
    #jet_pf_cands_idx_info = eval("events.JetPFCands"+jet_type)
    #jet_pf_cands = events.JetPFCands[jet_pf_cands_idx_info.candIdx]
    #pf_cands = ak.zip(
    #    {
    #        "jetIdx": jet_pf_cands_idx_info.jetIdx,
    #        "pt"    : jet_pf_cands.pt,
    #        "mass"  : jet_pf_cands.mass,
    #        "eta"   : jet_pf_cands.eta,
    #        "phi"   : jet_pf_cands.phi,
    #        "charge": jet_pf_cands.charge,
    #        "pdgId" : jet_pf_cands.pdgId,
    #    },
    #    with_name="PtEtaPhiMLorentzVector",
    #)

    return pf_cands


def get_collection_size(events, collection):
    return ak.count(getattr(events, collection).pt, axis=1)


class BranchesProducer(processor.ProcessorABC):
    """Processor for computing new branches to be added to a new ROOT file.

    The class must implement the following methods:
        * __init__: Define the branches to be computed
        * process: Calculation of the new quantities and accumulation into
                   the branches.
    """


    def __init__(self):
        """Define the branches to be computed.

        New branches to be computed are defined in the class attribute
        _accumulator (coffea.processor.dict_accumulator).

        Args:
            None
        
        Returns:
            None
        """

        self.jet_types = ["AK4", "AK8"]  # capitals used because capitals are used in PFNanoAOD

        
        branches = {}

        branches["genWeight"]: cfutl.column_accumulator(np.float64)

        for jet_type in self.jet_types:
            jet_collection = jet_algo_2_jet_collection(jet_type)
            branches["n" + jet_collection]: cfutl.column_accumulator(np.int64)
            branches[jet_collection + "_ptD"]: cfutl.column_accumulator(np.float64)
            branches[jet_collection + "_girth"]: cfutl.column_accumulator(np.float64)
            branches[jet_collection + "_axisMajor"]: cfutl.column_accumulator(np.float64)
            branches[jet_collection + "_axisMinor"]: cfutl.column_accumulator(np.float64)
            branches[jet_collection + "_axisAvg"]: cfutl.column_accumulator(np.float64)
            if jet_type == "AK8":
                branches[jet_collection + "_chHEF"]: cfutl.column_accumulator(np.float64)
                branches[jet_collection + "_neHEF"]: cfutl.column_accumulator(np.float64)
 
        ## Define accumulator
        self._accumulator = processor.dict_accumulator({
            **branches
            })


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        """Compute new quantities to be added to the ROOT NTuple.

        Args:
            events (ak.Array): the Events TTree opened with uproot.

        Returns:
            dict[str, coffea.processor.dict_accumulator]
        """

        ## Define accumulator
        output = self.accumulator.identity()

        output["genWeight"] = cfutl.accumulate(events.genWeight)

        for jet_type in self.jet_types:

            jet_collection = jet_algo_2_jet_collection(jet_type)

            jets = get_jets(events, jet_type)
            jet_pf_cands = get_pf_cands(events, jet_type)
            njets = get_collection_size(events, jet_collection)

            output["n" + jet_collection] = cfutl.accumulate(njets)
            output[jet_collection + "_ptD"] = cfutl.accumulate(jetvar.calculate_ptD(jet_pf_cands)[0])
            output[jet_collection + "_girth"] = cfutl.accumulate(jetvar.calculate_girth(jet_pf_cands, jets, njets)[0])
            axis_major, axis_minor, axis_avg = jetvar.calculate_axes(jet_pf_cands, jets, njets)
            output[jet_collection + "_axisMajor"] = cfutl.accumulate(axis_major)
            output[jet_collection + "_axisMinor"] = cfutl.accumulate(axis_minor)
            output[jet_collection + "_axisAvg"] = cfutl.accumulate(axis_avg)

            if jet_type == "AK8":
                output[jet_collection + "_chHEF"] = cfutl.accumulate(jetvar.calculate_chHEF(jet_pf_cands, jets, njets)[0])
                output[jet_collection + "_neHEF"] = cfutl.accumulate(jetvar.calculate_neHEF(jet_pf_cands, jets, njets)[0])

        return output


    def postprocess(self, accumulator):
        return accumulator

