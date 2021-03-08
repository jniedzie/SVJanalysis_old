from coffea import processor
import awkward as ak
import numpy as np
from collections import OrderedDict

import objects as obj


MET_FILTERS = (
    "goodVertices",
    "globalSuperTightHalo2016Filter",
    "HBHENoiseFilter",
    "HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "BadPFMuonFilter",
    "eeBadScFilter",
    "ecalBadCalibFilterV2",
    )


def np_acc_int():
    return processor.column_accumulator(np.array([], dtype=np.int64))

def np_acc_float():
    return processor.column_accumulator(np.array([], dtype=np.float64))

def np_acc_object():
    return processor.column_accumulator(np.array([], dtype=object))



class Preselection(processor.ProcessorABC):
    """
    Make jagged array with selected events and objects.

    Selects events:
       * passing higest efficiency trigger
       * passing MET filters: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    Selects objects:
       * AK8 jets satisfying pt > 200 GeV and |eta| < 2.4
       * AK4 jets satisfying pt > 30  GeV and |eta| < 2.4
       * PF candidates for good AK4 and AK8 jets

    """

    def __init__(self):
        """Define all variables to be stored. Keep same data structure as PFnanoAOD"""

        self._accumulator = processor.dict_accumulator({
            "nFatJet"         : np_acc_int(),
            "FatJet_pt"       : np_acc_object(),
            "FatJet_eta"      : np_acc_object(),
            "FatJet_phi"      : np_acc_object(),
            "FatJet_mass"     : np_acc_object(),
            "FatJet_msoftdrop": np_acc_object(),
            "FatJet_n2b1"     : np_acc_object(),
            "FatJet_n3b1"     : np_acc_object(),
            "FatJet_tau1"     : np_acc_object(),
            "FatJet_tau2"     : np_acc_object(),
            "FatJet_tau3"     : np_acc_object(),
            "FatJet_tau4"     : np_acc_object(),

            "nJet"     : np_acc_int(),
            "Jet_mass" : np_acc_object(),
            "Jet_pt"   : np_acc_object(),
            "Jet_eta"  : np_acc_object(),
            "Jet_phi"  : np_acc_object(),
            "Jet_chHEF": np_acc_object(),
            "Jet_neHEF": np_acc_object(),

            "nJetPFCands"       : np_acc_int(),
            "JetPFCands_jetIdx" : np_acc_object(),
            "JetPFCands_eta"    : np_acc_object(),
            "JetPFCands_mass"   : np_acc_object(),
            "JetPFCands_phi"    : np_acc_object(),
            "JetPFCands_pt"     : np_acc_object(),
            "JetPFCands_trkChi2": np_acc_object(),
            "JetPFCands_vtxChi2": np_acc_object(),
            "JetPFCands_charge" : np_acc_object(),
            "JetPFCands_pdgId"  : np_acc_object(),
        
            "nFatJetPFCands"       : np_acc_int(),
            "FatJetPFCands_jetIdx" : np_acc_object(),
            "FatJetPFCands_eta"    : np_acc_object(),
            "FatJetPFCands_mass"   : np_acc_object(),
            "FatJetPFCands_phi"    : np_acc_object(),
            "FatJetPFCands_pt"     : np_acc_object(),
            "FatJetPFCands_trkChi2": np_acc_object(),
            "FatJetPFCands_vtxChi2": np_acc_object(),
            "FatJetPFCands_charge" : np_acc_object(),
            "FatJetPFCands_pdgId"  : np_acc_object(),

            "MET_phi"         : np_acc_float(),
            "MET_pt"          : np_acc_float(),
            "MET_significance": np_acc_float(),
            "MET_sumEt"       : np_acc_float(),
            "PuppiMET_phi"    : np_acc_float(),
            "PuppiMET_pt"     : np_acc_float(),
            "PuppiMET_sumEt"  : np_acc_float(),
            "RawMET_phi"      : np_acc_float(),
            "RawMET_pt"       : np_acc_float(),
            "RawMET_sumEt"    : np_acc_float(),

            "cutflow": processor.dict_accumulator(
                OrderedDict(
                    **{
                        "all": processor.value_accumulator(int, initial=0),
                        "trigger": processor.value_accumulator(int, initial=0),
                    },
                    **{
                        METFilter: processor.value_accumulator(int, initial=0) for METFilter in MET_FILTERS
                    }
                    )
                ),

            })



    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        """Apply cuts."""

        output = self.accumulator.identity()

        ## Event selection
        output["cutflow"]["all"] += len(events)

        # Highest efficiency HLT trigger
        events = events[events.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60]
        output["cutflow"]["trigger"] += len(events)

        # MET filters
        for METFilter in MET_FILTERS:
            events = events[getattr(events, "Flag_"+METFilter)]
            output["cutflow"][METFilter] += len(events)

        # Lepton veto


        ## Good object selection
        # Good Jet
        goodFatJets = obj.GoodFatJets(events)
        goodFatJetPFCands = obj.GoodFatJetPFCands(events, goodFatJets)

        # Good Jet
        goodJets = obj.GoodJets(events)
        goodJetPFCands = obj.GoodJetPFCands(events, goodJets)
        

#        ## Example of difference between the two schemas
#        ## "schema": BaseSchema
#        output["muon_pt"] += processor.column_accumulator(ak.to_numpy(ak.flatten(events.Muon_pt)))
#        ## "schema": NanoAODSchema
#        output["muon_pt"] += processor.column_accumulator(ak.to_numpy(ak.flatten(events.Muon.pt)))
#
#        ## Without object pre-selection
#        output["nFatJet"] += processor.column_accumulator(ak.to_numpy(ak.flatten(events.nFatJet, axis=None)))
#
#        output["FatJet_pt"]        += processor.column_accumulator(np.array(ak.to_list(events.FatJet_pt)       , dtype=object))
#        output["FatJet_eta"]       += processor.column_accumulator(np.array(ak.to_list(events.FatJet_eta)      , dtype=object))
#        output["FatJet_phi"]       += processor.column_accumulator(np.array(ak.to_list(events.FatJet_phi)      , dtype=object))
#        output["FatJet_mass"]      += processor.column_accumulator(np.array(ak.to_list(events.FatJet_mass)     , dtype=object))
#        output["FatJet_msoftdrop"] += processor.column_accumulator(np.array(ak.to_list(events.FatJet_msoftdrop), dtype=object))
#        output["FatJet_n2b1"]      += processor.column_accumulator(np.array(ak.to_list(events.FatJet_n2b1)     , dtype=object))
#        output["FatJet_n3b1"]      += processor.column_accumulator(np.array(ak.to_list(events.FatJet_n3b1)     , dtype=object))
#        output["FatJet_tau1"]      += processor.column_accumulator(np.array(ak.to_list(events.FatJet_tau1)     , dtype=object))
#        output["FatJet_tau2"]      += processor.column_accumulator(np.array(ak.to_list(events.FatJet_tau2)     , dtype=object))
#        output["FatJet_tau3"]      += processor.column_accumulator(np.array(ak.to_list(events.FatJet_tau3)     , dtype=object))
#        output["FatJet_tau4"]      += processor.column_accumulator(np.array(ak.to_list(events.FatJet_tau4)     , dtype=object))
#
#        output["nFatJetPFCands"]       += processor.column_accumulator(ak.to_numpy(ak.flatten(events.nFatJetPFCands, axis=None)))
#        output["FatJetPFCands_jetIdx"] += processor.column_accumulator(np.array(ak.to_list(events.FatJetPFCands_jetIdx) , dtype=object))
#        output["FatJetPFCands_eta"]    += processor.column_accumulator(np.array(ak.to_list(events.FatJetPFCands_eta)    , dtype=object))
#        output["FatJetPFCands_mass"]   += processor.column_accumulator(np.array(ak.to_list(events.FatJetPFCands_mass)   , dtype=object))
#        output["FatJetPFCands_phi"]    += processor.column_accumulator(np.array(ak.to_list(events.FatJetPFCands_phi)    , dtype=object))
#        output["FatJetPFCands_pt"]     += processor.column_accumulator(np.array(ak.to_list(events.FatJetPFCands_pt)     , dtype=object))
#        output["FatJetPFCands_trkChi2"]+= processor.column_accumulator(np.array(ak.to_list(events.FatJetPFCands_trkChi2), dtype=object))
#        output["FatJetPFCands_vtxChi2"]+= processor.column_accumulator(np.array(ak.to_list(events.FatJetPFCands_vtxChi2), dtype=object))
#        output["FatJetPFCands_charge"] += processor.column_accumulator(np.array(ak.to_list(events.FatJetPFCands_charge) , dtype=object))
#        output["FatJetPFCands_pdgId"]  += processor.column_accumulator(np.array(ak.to_list(events.FatJetPFCands_pdgId)  , dtype=object))


        # Fat Jets
        output["nFatJet"] += processor.column_accumulator(ak.to_numpy(ak.flatten(goodFatJets.n, axis=None)))

        output["FatJet_pt"]        += processor.column_accumulator(np.array(ak.to_list(goodFatJets.pt)       , dtype=object))
        output["FatJet_eta"]       += processor.column_accumulator(np.array(ak.to_list(goodFatJets.eta)      , dtype=object))
        output["FatJet_phi"]       += processor.column_accumulator(np.array(ak.to_list(goodFatJets.phi)      , dtype=object))
        output["FatJet_mass"]      += processor.column_accumulator(np.array(ak.to_list(goodFatJets.mass)     , dtype=object))
        output["FatJet_msoftdrop"] += processor.column_accumulator(np.array(ak.to_list(goodFatJets.msoftdrop), dtype=object))
        output["FatJet_n2b1"]      += processor.column_accumulator(np.array(ak.to_list(goodFatJets.n2b1)     , dtype=object))
        output["FatJet_n3b1"]      += processor.column_accumulator(np.array(ak.to_list(goodFatJets.n3b1)     , dtype=object))
        output["FatJet_tau1"]      += processor.column_accumulator(np.array(ak.to_list(goodFatJets.tau1)     , dtype=object))
        output["FatJet_tau2"]      += processor.column_accumulator(np.array(ak.to_list(goodFatJets.tau2)     , dtype=object))
        output["FatJet_tau3"]      += processor.column_accumulator(np.array(ak.to_list(goodFatJets.tau3)     , dtype=object))
        output["FatJet_tau4"]      += processor.column_accumulator(np.array(ak.to_list(goodFatJets.tau4)     , dtype=object))

        # Fat Jet PF candidates
        output["nFatJetPFCands"]       += processor.column_accumulator(ak.to_numpy(ak.flatten(goodFatJetPFCands.n, axis=None)))
        output["FatJetPFCands_jetIdx"] += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.jetIdx) , dtype=object))
        output["FatJetPFCands_eta"]    += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.eta)    , dtype=object))
        output["FatJetPFCands_mass"]   += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.mass)   , dtype=object))
        output["FatJetPFCands_phi"]    += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.phi)    , dtype=object))
        output["FatJetPFCands_pt"]     += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.pt)     , dtype=object))
        output["FatJetPFCands_trkChi2"]+= processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.trkChi2), dtype=object))
        output["FatJetPFCands_vtxChi2"]+= processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.vtxChi2), dtype=object))
        output["FatJetPFCands_charge"] += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.charge) , dtype=object))
        output["FatJetPFCands_pdgId"]  += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.pdgId)  , dtype=object))


        # AK4 Jets
        output["nJet"] += processor.column_accumulator(ak.to_numpy(ak.flatten(goodJets.n, axis=None)))

        output["Jet_mass"]  += processor.column_accumulator(np.array(ak.to_list(goodJets.mass) , dtype=object))
        output["Jet_pt"]    += processor.column_accumulator(np.array(ak.to_list(goodJets.pt)   , dtype=object))
        output["Jet_eta"]   += processor.column_accumulator(np.array(ak.to_list(goodJets.eta)  , dtype=object))
        output["Jet_phi"]   += processor.column_accumulator(np.array(ak.to_list(goodJets.phi)  , dtype=object))
        output["Jet_chHEF"] += processor.column_accumulator(np.array(ak.to_list(goodJets.chHEF), dtype=object))
        output["Jet_neHEF"] += processor.column_accumulator(np.array(ak.to_list(goodJets.neHEF), dtype=object))

        # AK4 Jet PF candidates
        output["nJetPFCands"]       += processor.column_accumulator(ak.to_numpy(ak.flatten(goodJetPFCands.n, axis=None)))
        output["JetPFCands_jetIdx"] += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.jetIdx) , dtype=object))
        output["JetPFCands_eta"]    += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.eta)    , dtype=object))
        output["JetPFCands_mass"]   += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.mass)   , dtype=object))
        output["JetPFCands_phi"]    += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.phi)    , dtype=object))
        output["JetPFCands_pt"]     += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.pt)     , dtype=object))
        output["JetPFCands_trkChi2"]+= processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.trkChi2), dtype=object))
        output["JetPFCands_vtxChi2"]+= processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.vtxChi2), dtype=object))
        output["JetPFCands_charge"] += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.charge) , dtype=object))
        output["JetPFCands_pdgId"]  += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.pdgId)  , dtype=object))


        # MET
        output["MET_phi"]          += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_phi         , axis=None)))
        output["MET_pt"]           += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_pt          , axis=None)))
        output["MET_significance"] += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_significance, axis=None)))
        output["MET_sumEt"]        += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_sumEt       , axis=None)))
        output["PuppiMET_phi"]     += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_phi         , axis=None)))
        output["PuppiMET_pt"]      += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_pt          , axis=None)))
        output["PuppiMET_sumEt"]   += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_sumEt       , axis=None)))
        output["RawMET_phi"]       += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_phi         , axis=None)))
        output["RawMET_pt"]        += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_pt          , axis=None)))
        output["RawMET_sumEt"]     += processor.column_accumulator(ak.to_numpy(ak.flatten(events.MET_sumEt       , axis=None)))


        return output


    def postprocess(self, accumulator):
        return accumulator

