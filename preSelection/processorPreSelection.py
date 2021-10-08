import sys

import awkward as ak
from coffea import processor

sys.path.append("../utilities/")
import nameUtilities as nameutl
import coffeaUtilities as cfutl
import objects as obj
from metFilters import met_filters
from triggers import selected_triggers


def capitalize(word):
    return word[0].upper() + word[1:]


def get_number_of_events(events):
    return ak.sum(events.genWeight)


def apply_trigger_cuts(events, accumulator, trigger_list):
    """Filter events using an or of all triggers.

    Args:
        events (awkward.highlevel.Array): the Events TTree opened with uproot
        accumulator (coffea.processor.accumulator.dict_accumulator)

    Returns:
        awkward.highlevel.Array
    """

    trigger_selection = None
    for trigger in trigger_list:
        trigger_branch = getattr(events.HLT, trigger)
        if trigger_selection is None:
            trigger_selection = trigger_branch
        else:
            trigger_selection = (trigger_selection | trigger_branch)
    events = events[trigger_selection]
    accumulator["cutflow"]["Trigger"] += get_number_of_events(events)

    return events


def apply_met_filters_cuts(events, accumulator, met_filters):
    """MET filters cuts.

    Args:
        events (awkward.highlevel.Array): the Events TTree opened with uproot
        accumulator (coffea.processor.accumulator.dict_accumulator)

    Returns:
        awkward.highlevel.Array
    """

    for met_filter in met_filters:
        events = events[getattr(events.Flag, met_filter)]
        accumulator["cutflow"][capitalize(met_filter)] = get_number_of_events(events)

    return events


def fill_collection_variables(accumulator, objects, collection_name):
    """Add collection variables to accumulator.

    Args:
        accumulator (coffea.processor.accumulator.dict_accumulator)
        objects (Jets, PfCands, JetPfCandsMatchingTable, Met)
        collection_name (str)

    Returns:
        None
    """

    if hasattr(objects, "n"):
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


def fill_single_event_level_variables(accumulator, single_event_level_collections):
    """Add single event level variables to accumulator.

    Args:
        accumulator (coffea.processor.accumulator.dict_accumulator)
        single_event_level_collections (SingleEventLevel)

    Returns:
        None
    """

    for variable in single_event_level_collections.variables:
        accumulator[variable] = cfutl.accumulate(getattr(single_event_level_collections, variable))

    return



class Preselection_tchannel(processor.ProcessorABC):
    """Make coffea accumulator containing events branches and cutflow.

    Selects objects:
       * AK8 jets satisfying pt > 200 GeV and |eta| < 2.4
       * AK4 jets satisfying pt > 30  GeV and |eta| < 2.4
    Selects events:
       * passing all MET and JetHT triggers
       * passing MET filters: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
       * having 2 good AK4 jets
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
        accumulator["cutflow"]["NoCut"] = get_number_of_events(events)

        ## Event selection
        events = apply_trigger_cuts(events, accumulator, selected_triggers)
        events = apply_met_filters_cuts(events, accumulator, met_filters)

        # Good jets
        good_ak4_jets = obj.Jets(events, "AK4", self.input_file_type, cut="(jet.pt > 30 ) & (np.abs(jet.eta) < 2.4)")
        good_ak8_jets = obj.Jets(events, "AK8", self.input_file_type, cut="(jet.pt > 200) & (np.abs(jet.eta) < 2.4)")

        # Keeping only events with more than 2 good AK4 jets
        filter_ = (good_ak4_jets.n >= 2)
        events = events[filter_]
        good_ak4_jets.apply_cut(filter_)
        good_ak8_jets.apply_cut(filter_)
        accumulator["cutflow"]["2GoodAk4Jets"] = get_number_of_events(events)

        # Good PF candidates
        good_pf_cands = obj.PfCands(events, self.input_file_type)
        good_ak4_jet_pf_cands_matching_table = obj.JetPfCandsMatchingTable(events, good_ak4_jets, "AK4", self.input_file_type)
        good_ak8_jet_pf_cands_matching_table = obj.JetPfCandsMatchingTable(events, good_ak8_jets, "AK8", self.input_file_type)

        # MET
        met_flavors = ["MET", "PuppiMET", "RawMET"]
        mets = { met_flavor: obj.Met(events, met_flavor, self.input_file_type) for met_flavor in met_flavors }

        # Single event level quantities
        single_event_level = obj.SingleEventLevel(events, self.input_file_type)

        if ak.count(events.genWeight) > 0:  # if there are events left after pre-selection cuts for that chunk
            fill_single_event_level_variables(accumulator, single_event_level)
            for met_flavor, met in mets.items():
                fill_collection_variables(accumulator, met, met_flavor)
            fill_jets_variables(accumulator, good_ak4_jets, "AK4")
            fill_jets_variables(accumulator, good_ak8_jets, "AK8")
            fill_pf_cands_variables(accumulator, good_pf_cands)
            fill_jet_pf_cands_matching_table(accumulator, good_ak4_jet_pf_cands_matching_table, "AK4")
            fill_jet_pf_cands_matching_table(accumulator, good_ak8_jet_pf_cands_matching_table, "AK8")

        # Explicitly freeing the memory to avoid memory leaks
        del events
        del good_ak4_jets
        del good_ak8_jets
        del good_pf_cands
        del good_ak4_jet_pf_cands_matching_table
        del good_ak8_jet_pf_cands_matching_table
        del mets
        del single_event_level

        return accumulator


    def postprocess(self, accumulator):
        return accumulator

