from coffea import processor
import awkward as ak
import sys

sys.path.append("../utilities/")
import nameUtilities as nameutl
import coffeaUtilities as cfutl

import objects as obj


MET_FILTERS = (
    "goodVertices",
    "globalSuperTightHalo2016Filter",
    "HBHENoiseFilter",
    "HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "BadPFMuonFilter",
    "eeBadScFilter",
    )



def apply_tchannel_cuts(events, accumulator):
    """t channel pre-selection cuts.

    Args:
        events (awkward.highlevel.Array): the Events TTree opened with uproot
        accumulator (coffea.processor.accumulator.dict_accumulator)

    Returns:
        None
    """

    ## Highest efficiency HLT trigger
    events = events[events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60]
    accumulator["cutflow"]["trigger"] = ak.sum(events.genWeight)

    ## MET filters
    for METFilter in MET_FILTERS:
        events = events[getattr(events.Flag, METFilter)]
        accumulator["cutflow"][METFilter] = ak.sum(events.genWeight)

    return events


def fill_collection_variables(accumulator, objects, collection_name):
    """Add collection variables to accumulator.

    Args:
        accumulator (coffea.processor.accumulator.dict_accumulator)
        objects (Jets, PfCands, JetPfCandsMatchingTable)
        collection_name (str)

    Returns:
        None
    """

    accumulator["n"+collection_name] = cfutl.accumulate(objects.n)

    for variable in objects.variables:
        accumulator[collection_name + "_" + variable] = cfutl.accumulate(getattr(objects, variable))
 
    return


def fill_jets_variables(accumulator, good_jets, jet_algo_name):
    """Add jet variables to accumulator.

    Args:
        accumulator (coffea.processor.accumulator.dict_accumulator)
        good_jets (Jets)
        jet_algo_name (str)

    Returns:
        None
    """

    jet_collection_name = nameutl.jet_algo_name_to_jet_collection_name(jet_algo_name)
    fill_collection_variables(accumulator, good_jets, jet_collection_name)

    return


def fill_pf_cands_variables(accumulator, good_pf_cands):
    """Add PF candidates variables to accumulator.

    Args:
        accumulator (coffea.processor.accumulator.dict_accumulator)
        good_pf_cands (PfCands)

    Returns:
        None
    """

    fill_collection_variables(accumulator, good_pf_cands, "PFCands")

    return
 

def fill_jet_pf_cands_matching_table(accumulator, good_jet_pf_cands_matching_table, jet_algo_name):
    """Add jet PF candidates matching table to accumulator.

    Args:
        accumulator (coffea.processor.accumulator.dict_accumulator)
        good_jet_pf_cands_matching_table (JetPfCandsMatchingTable)
        jet_algo_name (str)

    Returns:
        None
    """

    jet_collection_name = nameutl.jet_algo_name_to_jet_collection_name(jet_algo_name)
    fill_collection_variables(accumulator, good_jet_pf_cands_matching_table, jet_collection_name + "PFCands")

    return


def fill_event_variables(accumulator, events):
    """Add event variables to accumulator.

    Args:
        accumulator (coffea.processor.accumulator.dict_accumulator)
        events (awkward.highlevel.Array): the Events TTree opened with uproot

    Returns:
        None
    """

    accumulator["MET_phi"]          = cfutl.accumulate(events.MET.phi)
    accumulator["MET_pt"]           = cfutl.accumulate(events.MET.pt)
    accumulator["MET_significance"] = cfutl.accumulate(events.MET.significance)
    accumulator["MET_sumEt"]        = cfutl.accumulate(events.MET.sumEt)
    accumulator["PuppiMET_phi"]     = cfutl.accumulate(events.MET.phi)
    accumulator["PuppiMET_pt"]      = cfutl.accumulate(events.MET.pt)
    accumulator["PuppiMET_sumEt"]   = cfutl.accumulate(events.MET.sumEt)
    accumulator["RawMET_phi"]       = cfutl.accumulate(events.MET.phi)
    accumulator["RawMET_pt"]        = cfutl.accumulate(events.MET.pt)
    accumulator["RawMET_sumEt"]     = cfutl.accumulate(events.MET.sumEt)

    accumulator["genWeight"]        = cfutl.accumulate(events.genWeight)

    return



class Preselection_tchannel(processor.ProcessorABC):
    """Make coffea accumulator containing events branches and cutflow.

    Selects events:
       * passing higest efficiency trigger (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, evaluated on baseline model)
       * passing MET filters: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    Selects objects:
       * AK8 jets satisfying pt > 200 GeV and |eta| < 2.4
       * AK4 jets satisfying pt > 30  GeV and |eta| < 2.4
    """

    def __init__(self, input_file_type="PFNanoAOD_106X_v02"):
        """
        Args:
            input_file_type (str)

        Returns:
            None
        """

        self.input_file_type = input_file_type

        # By using dict_accumulator, the branches type do not have to be defined
        self._accumulator = processor.dict_accumulator()


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        """Apply cuts and fill accumulator.

        Args:
	    events (awkward.highlevel.Array): the Events TTree opened with uproot

        Returns:
	    coffea.processor.accumulator.dict_accumulator
        """

        accumulator = self.accumulator.identity()

        accumulator["cutflow"] = processor.defaultdict_accumulator(float)
        accumulator["cutflow"]["all"] = ak.sum(events.genWeight)

        ## Event selection
        events = apply_tchannel_cuts(events, accumulator)
        good_ak4_jets = obj.Jets(events, "AK4", self.input_file_type, cut="(jet.pt > 30 ) & (np.abs(jet.eta) < 2.4)")
        good_ak8_jets = obj.Jets(events, "AK8", self.input_file_type, cut="(jet.pt > 200) & (np.abs(jet.eta) < 2.4)")

        # Good PF candidates
        good_pf_cands = obj.PfCands(events, self.input_file_type)
        good_ak4_jet_pf_cands_matching_table = obj.JetPfCandsMatchingTable(events, good_ak4_jets, "AK4", self.input_file_type)
        good_ak8_jet_pf_cands_matching_table = obj.JetPfCandsMatchingTable(events, good_ak8_jets, "AK8", self.input_file_type)

        if ak.count(events.genWeight) > 0:  # if there are events left after pre-selection cuts for that chunk
            fill_event_variables(accumulator, events)
            fill_jets_variables(accumulator, good_ak4_jets, "AK4")
            fill_jets_variables(accumulator, good_ak8_jets, "AK8")
            fill_pf_cands_variables(accumulator, good_pf_cands)
            fill_jet_pf_cands_matching_table(accumulator, good_ak4_jet_pf_cands_matching_table, "AK4")
            fill_jet_pf_cands_matching_table(accumulator, good_ak8_jet_pf_cands_matching_table, "AK8")

        return accumulator


    def postprocess(self, accumulator):
        return accumulator

