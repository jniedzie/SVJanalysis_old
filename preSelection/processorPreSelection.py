from coffea import processor
from coffea.nanoevents.methods import vector
import awkward as ak
import numpy as np
from collections import OrderedDict
import sys

sys.path.append("../utilities/")
import nameUtilities as nameutl
import coffeaUtilities as cfutl
import PtEtaPhiMLorentzVectorUtilities as vecutl

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



def apply_tchannel_cuts(events, output):
    """t channel pre-selection cuts."""

    ## Highest efficiency HLT trigger
    events = events[events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60]
    output["cutflow"]["trigger"] += ak.sum(events.genWeight)

    ## MET filters
    for METFilter in MET_FILTERS:
        events = events[getattr(events.Flag, METFilter)]
        output["cutflow"][METFilter] += ak.sum(events.genWeight)

    ## Lepton veto
    # to be written...

    return events


def init_pf_cands_dict():

    pf_cands_dict = {
        "nPFCands"       : cfutl.column_accumulator(np.int64),
        "PFCands_eta"    : cfutl.column_accumulator(object),
        "PFCands_mass"   : cfutl.column_accumulator(object),
        "PFCands_phi"    : cfutl.column_accumulator(object),
        "PFCands_pt"     : cfutl.column_accumulator(object),
        "PFCands_trkChi2": cfutl.column_accumulator(object),
        "PFCands_vtxChi2": cfutl.column_accumulator(object),
        "PFCands_charge" : cfutl.column_accumulator(object),
        "PFCands_pdgId"  : cfutl.column_accumulator(object),
        }

    return pf_cands_dict


def init_ak4_jets_dicts():

    pf_cands_matching_dict = {
        "nJetPFCands"          : cfutl.column_accumulator(np.int64),
        "JetPFCands_jetIdx"    : cfutl.column_accumulator(object),
        "JetPFCands_pFCandsIdx": cfutl.column_accumulator(object),
        }

    jet_dict = {
        "nJet"         : cfutl.column_accumulator(np.int64),
        "Jet_pt"       : cfutl.column_accumulator(object),
        "Jet_eta"      : cfutl.column_accumulator(object),
        "Jet_phi"      : cfutl.column_accumulator(object),
        "Jet_mass"     : cfutl.column_accumulator(object),
        "Jet_chHEF"    : cfutl.column_accumulator(object),
        "Jet_neHEF"    : cfutl.column_accumulator(object),
    }

    return jet_dict, pf_cands_matching_dict


def init_ak8_jets_dicts():

    pf_cands_matching_dict = {
        "nFatJetPFCands"          : cfutl.column_accumulator(np.int64),
        "FatJetPFCands_jetIdx"    : cfutl.column_accumulator(object),
        "FatJetPFCands_pFCandsIdx": cfutl.column_accumulator(object),
        }

    jet_dict = {
        "nFatJet"         : cfutl.column_accumulator(np.int64),
        "FatJet_pt"       : cfutl.column_accumulator(object),
        "FatJet_eta"      : cfutl.column_accumulator(object),
        "FatJet_phi"      : cfutl.column_accumulator(object),
        "FatJet_mass"     : cfutl.column_accumulator(object),
        "FatJet_msoftdrop": cfutl.column_accumulator(object),
        "FatJet_n2b1"     : cfutl.column_accumulator(object),
        "FatJet_n3b1"     : cfutl.column_accumulator(object),
        "FatJet_tau1"     : cfutl.column_accumulator(object),
        "FatJet_tau2"     : cfutl.column_accumulator(object),
        "FatJet_tau3"     : cfutl.column_accumulator(object),
        "FatJet_tau4"     : cfutl.column_accumulator(object),
    }

    return jet_dict, pf_cands_matching_dict


def init_event_variables_dict():

    event_variables_dict = {
        # MET
        "MET_phi"         : cfutl.column_accumulator(np.float64),
        "MET_pt"          : cfutl.column_accumulator(np.float64),
        "MET_significance": cfutl.column_accumulator(np.float64),
        "MET_sumEt"       : cfutl.column_accumulator(np.float64),
        "PuppiMET_phi"    : cfutl.column_accumulator(np.float64),
        "PuppiMET_pt"     : cfutl.column_accumulator(np.float64),
        "PuppiMET_sumEt"  : cfutl.column_accumulator(np.float64),
        "RawMET_phi"      : cfutl.column_accumulator(np.float64),
        "RawMET_pt"       : cfutl.column_accumulator(np.float64),
        "RawMET_sumEt"    : cfutl.column_accumulator(np.float64),
        # Event weight
        "genWeight"       : cfutl.column_accumulator(np.float64),
    }

    return event_variables_dict



def init_accumulator(dicts):

    accumulator = processor.dict_accumulator({
        ## Events and jet collections variables
        **dicts,

        ## Cut efficiencies
        **{ "cutflow": processor.defaultdict_accumulator(int) },

        })

    return accumulator


def fill_collection_variables(output, objects, collection_name):

    output["n"+collection_name] += cfutl.accumulate(objects.n)

    for variable in objects.variables:
        output[collection_name + "_" + variable] += cfutl.accumulate(getattr(objects, variable))
 
    return

def fill_jets_variables(output, good_jets, jet_algo_name):

    jet_collection_name = nameutl.jet_algo_name_to_jet_collection_name(jet_algo_name)
    fill_collection_variables(output, good_jets, jet_collection_name)

    return


def fill_pf_cands_variables(output, good_pf_cands):

    fill_collection_variables(output, good_pf_cands, "PFCands")

    return
 

def fill_jet_pf_cands_matching_table(output, good_jet_pf_cands_matching_table, jet_algo_name):

    jet_collection_name = nameutl.jet_algo_name_to_jet_collection_name(jet_algo_name)
    fill_collection_variables(output, good_jet_pf_cands_matching_table, jet_collection_name + "PFCands")

    return


def fill_event_variables(output, events):

    #output["MET_phi"]          += cfutl.accumulate(ak.to_numpy(ak.flatten(events.MET_phi         , axis=None)))
    output["MET_phi"]          += cfutl.accumulate(events.MET.phi)
    output["MET_pt"]           += cfutl.accumulate(events.MET.pt)
    output["MET_significance"] += cfutl.accumulate(events.MET.significance)
    output["MET_sumEt"]        += cfutl.accumulate(events.MET.sumEt)
    output["PuppiMET_phi"]     += cfutl.accumulate(events.MET.phi)
    output["PuppiMET_pt"]      += cfutl.accumulate(events.MET.pt)
    output["PuppiMET_sumEt"]   += cfutl.accumulate(events.MET.sumEt)
    output["RawMET_phi"]       += cfutl.accumulate(events.MET.phi)
    output["RawMET_pt"]        += cfutl.accumulate(events.MET.pt)
    output["RawMET_sumEt"]     += cfutl.accumulate(events.MET.sumEt)

    output["genWeight"]        += cfutl.accumulate(events.genWeight)

    return



class Preselection_tchannel(processor.ProcessorABC):
    """
    Make jagged array with selected events and objects.

    Selects events:
       * passing higest efficiency trigger (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, evaluated on baseline model)
       * passing MET filters: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    Selects objects:
       * AK8 jets satisfying pt > 200 GeV and |eta| < 2.4
       * AK4 jets satisfying pt > 30  GeV and |eta| < 2.4
       * PF candidates for good AK4 and AK8 jets (102X only)

    """

    def __init__(self, input_file_type="master"):
        """Define all variables to be stored. Transform data sturcture to match PF nano AOD master."""

        self.input_file_type = input_file_type

        ## Define accumulator
        pf_cands_dict = init_pf_cands_dict()
        ak4_jet_dict, ak4_jet_pf_cands_matching_dict = init_ak4_jets_dicts()
        ak8_jet_dict, ak8_jet_pf_cands_matching_dict = init_ak8_jets_dicts()
        event_variables_dict = init_event_variables_dict()
        self._accumulator = init_accumulator(
            {**event_variables_dict,
             **ak4_jet_dict,
             **ak8_jet_dict,
             **ak4_jet_pf_cands_matching_dict,
             **ak8_jet_pf_cands_matching_dict,
             **pf_cands_dict,
             }
        )


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        """Apply cuts."""

        output = self.accumulator.identity()

        output["cutflow"]["all"] += ak.sum(events.genWeight)

        ## Event selection
        events = apply_tchannel_cuts(events, output)
        good_ak4_jets = obj.GoodJets(events, "AK4", self.input_file_type, cut="(jet.pt > 30 ) & (np.abs(jet.eta) < 2.4)")
        good_ak8_jets = obj.GoodJets(events, "AK8", self.input_file_type, cut="(jet.pt > 200) & (np.abs(jet.eta) < 2.4)")

        # Good PF candidates
        good_pf_cands = obj.GoodPfCands(events, self.input_file_type)
        good_ak4_jet_pf_cands_matching_table = obj.GoodJetPfCandsMatchingTable(events, good_ak4_jets, "AK4", self.input_file_type)
        good_ak8_jet_pf_cands_matching_table = obj.GoodJetPfCandsMatchingTable(events, good_ak8_jets, "AK8", self.input_file_type)

        if ak.count(events.genWeight) > 0:  # if there are events left after pre-selection cuts for that chunk
            fill_event_variables(output, events)
            fill_jets_variables(output, good_ak4_jets, "AK4")
            fill_jets_variables(output, good_ak8_jets, "AK8")
            fill_pf_cands_variables(output, good_pf_cands)
            fill_jet_pf_cands_matching_table(output, good_ak4_jet_pf_cands_matching_table, "AK4")
            fill_jet_pf_cands_matching_table(output, good_ak8_jet_pf_cands_matching_table, "AK8")

        return output


    def postprocess(self, accumulator):
        return accumulator

