from coffea import processor
import awkward as ak
import numpy as np
from collections import OrderedDict
import sys

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


## Define some short-hands for column and value accumulator
def column_accumulator(type_):
    return processor.column_accumulator(np.array([], dtype=type_))

def value_accumulator(type_, initial=0):
    return processor.value_accumulator(type_, initial=initial)


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

    def __init__(self, inputFileType="PFnano102X"):
        """Define all variables to be stored. Keep same data structure as PFnanoAOD"""

        self.inputFileType = inputFileType

        ## Jets, jetPFCands collections and matching of the two vary between PFnano102X and PFnano106X
        if self.inputFileType == "PFnano102X":
            jetPFCandsDict = {
                "nJetPFCands"       : column_accumulator(np.int64),
                "JetPFCands_eta"    : column_accumulator(object),
                "JetPFCands_mass"   : column_accumulator(object),
                "JetPFCands_phi"    : column_accumulator(object),
                "JetPFCands_pt"     : column_accumulator(object),
                "JetPFCands_trkChi2": column_accumulator(object),
                "JetPFCands_vtxChi2": column_accumulator(object),
                "JetPFCands_charge" : column_accumulator(object),
                "JetPFCands_pdgId"  : column_accumulator(object),
            
                "nFatJetPFCands"       : column_accumulator(np.int64),
                "FatJetPFCands_eta"    : column_accumulator(object),
                "FatJetPFCands_mass"   : column_accumulator(object),
                "FatJetPFCands_phi"    : column_accumulator(object),
                "FatJetPFCands_pt"     : column_accumulator(object),
                "FatJetPFCands_trkChi2": column_accumulator(object),
                "FatJetPFCands_vtxChi2": column_accumulator(object),
                "FatJetPFCands_charge" : column_accumulator(object),
                "FatJetPFCands_pdgId"  : column_accumulator(object),
                }
            jetPFCandsMatchingDict = {
                "JetPFCands_jetIdx"   : column_accumulator(object),
                "FatJetPFCands_jetIdx": column_accumulator(object),
                }

        elif self.inputFileType == "PFnano106X":
            jetPFCandsDict = {
                "nJetPFCands"       : column_accumulator(np.int64),
                "JetPFCands_eta"    : column_accumulator(object),
                "JetPFCands_mass"   : column_accumulator(object),
                "JetPFCands_phi"    : column_accumulator(object),
                "JetPFCands_pt"     : column_accumulator(object),
                "JetPFCands_trkChi2": column_accumulator(object),
                "JetPFCands_vtxChi2": column_accumulator(object),
                "JetPFCands_charge" : column_accumulator(object),
                "JetPFCands_pdgId"  : column_accumulator(object),
                }
            jetPFCandsMatchingDict = {
                "nJetPFCandsAK4"       : column_accumulator(np.int64),
                "JetPFCandsAK4_jetIdx" : column_accumulator(object),
                "JetPFCandsAK4_candIdx": column_accumulator(object),
                "nJetPFCandsAK8"       : column_accumulator(np.int64),
                "JetPFCandsAK8_jetIdx" : column_accumulator(object),
                "JetPFCandsAK8_candIdx": column_accumulator(object),
                }

        else:
            print("ERROR: Unknown inputFileType %s" %inputFileType)
            sys.exit()


        ## Dictionary to store cut efficiencies
        cutflowDict = { "cutflow": processor.dict_accumulator(OrderedDict(
                          **{ "all": value_accumulator(int),
                              "trigger": value_accumulator(int),
                          },
                          **{
                              METFilter: value_accumulator(int) for METFilter in MET_FILTERS
                          }
                      ))}


        ## Define accumulator
        self._accumulator = processor.dict_accumulator({
            ## Event level variables
            **{
            # MET
            "MET_phi"         : column_accumulator(np.float64),
            "MET_pt"          : column_accumulator(np.float64),
            "MET_significance": column_accumulator(np.float64),
            "MET_sumEt"       : column_accumulator(np.float64),
            "PuppiMET_phi"    : column_accumulator(np.float64),
            "PuppiMET_pt"     : column_accumulator(np.float64),
            "PuppiMET_sumEt"  : column_accumulator(np.float64),
            "RawMET_phi"      : column_accumulator(np.float64),
            "RawMET_pt"       : column_accumulator(np.float64),
            "RawMET_sumEt"    : column_accumulator(np.float64),
            # Event weight
            "genWeight"       : column_accumulator(np.float64),
             },

            ## Jet collections variables
            **{
            "nFatJet"         : column_accumulator(np.int64),
            "FatJet_pt"       : column_accumulator(object),
            "FatJet_eta"      : column_accumulator(object),
            "FatJet_phi"      : column_accumulator(object),
            "FatJet_mass"     : column_accumulator(object),
            "FatJet_msoftdrop": column_accumulator(object),
            "FatJet_n2b1"     : column_accumulator(object),
            "FatJet_n3b1"     : column_accumulator(object),
            "FatJet_tau1"     : column_accumulator(object),
            "FatJet_tau2"     : column_accumulator(object),
            "FatJet_tau3"     : column_accumulator(object),
            "FatJet_tau4"     : column_accumulator(object),

            "nJet"     : column_accumulator(np.int64),
            "Jet_mass" : column_accumulator(object),
            "Jet_pt"   : column_accumulator(object),
            "Jet_eta"  : column_accumulator(object),
            "Jet_phi"  : column_accumulator(object),
            "Jet_chHEF": column_accumulator(object),
            "Jet_neHEF": column_accumulator(object),
            },

            ## JetPFCands variables
            **jetPFCandsDict,

            ## Jet and JetPFCands matching variables
            **jetPFCandsMatchingDict,

            ## Cutflow informations
            **cutflowDict,

            })



    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        """Apply cuts."""

        output = self.accumulator.identity()

        ## Event selection
        output["cutflow"]["all"] += ak.sum(events.genWeight)

        # Highest efficiency HLT trigger
        events = events[events.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60]
        output["cutflow"]["trigger"] += ak.sum(events.genWeight)

        # MET filters
        for METFilter in MET_FILTERS:
            events = events[getattr(events, "Flag_"+METFilter)]
            output["cutflow"][METFilter] += ak.sum(events.genWeight)

        # Lepton veto


        ## Good object selection
        # Good Jet
        goodFatJets = obj.GoodFatJets(events)
        goodFatJetPFCands = obj.GoodFatJetPFCands(events, goodFatJets, self.inputFileType)

        # Good Jet
        goodJets = obj.GoodJets(events)
        goodJetPFCands = obj.GoodJetPFCands(events, goodJets, self.inputFileType)
        

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

        # Fat Jet PF candidates - only in PFnano102X
        if self.inputFileType == "PFnano102X":
            output["nFatJetPFCands"]       += processor.column_accumulator(ak.to_numpy(ak.flatten(goodFatJetPFCands.n, axis=None)))
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
        output["JetPFCands_eta"]    += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.eta)    , dtype=object))
        output["JetPFCands_mass"]   += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.mass)   , dtype=object))
        output["JetPFCands_phi"]    += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.phi)    , dtype=object))
        output["JetPFCands_pt"]     += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.pt)     , dtype=object))
        output["JetPFCands_trkChi2"]+= processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.trkChi2), dtype=object))
        output["JetPFCands_vtxChi2"]+= processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.vtxChi2), dtype=object))
        output["JetPFCands_charge"] += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.charge) , dtype=object))
        output["JetPFCands_pdgId"]  += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.pdgId)  , dtype=object))

        # Jet PF candidates variables for matching with AK4/AK8 jet
        if self.inputFileType == "PFnano102X":
            output["JetPFCands_jetIdx"]    += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.jetIdx)    , dtype=object))
            output["FatJetPFCands_jetIdx"] += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.jetIdx) , dtype=object))
        elif self.inputFileType == "PFnano106X":
            output["nJetPFCandsAK4"] += processor.column_accumulator(np.array(goodJetPFCands.nAK4))
            output["JetPFCandsAK4_jetIdx"]  += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.jetIdx)    , dtype=object))
            output["JetPFCandsAK4_candIdx"] += processor.column_accumulator(np.array(ak.to_list(goodJetPFCands.candIdx)   , dtype=object))
            output["nJetPFCandsAK8"] += processor.column_accumulator(np.array(goodFatJetPFCands.nAK8))
            output["JetPFCandsAK8_jetIdx"]  += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.jetIdx) , dtype=object))
            output["JetPFCandsAK8_candIdx"] += processor.column_accumulator(np.array(ak.to_list(goodFatJetPFCands.candIdx), dtype=object))

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

        output["genWeight"]        += processor.column_accumulator(ak.to_numpy(ak.flatten(events.genWeight       , axis=None)))

        return output


    def postprocess(self, accumulator):
        return accumulator

