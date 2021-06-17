from coffea import processor
import awkward as ak
import numpy as np
import processorUtils as utl


class BranchesMaker(processor.ProcessorABC):
    """Example for calculating ptD."""


    def __init__(self, schema, pf_nano_aod_version):

        self.schema = schema
        self.pf_nano_aod_version = pf_nano_aod_version

        branches = { 
            "FatJet_ptD"      : utl.column_accumulator(np.float64),
            "nFatJet"         : utl.column_accumulator(np.float64),
        }
 
        ## Define accumulator
        self._accumulator = processor.dict_accumulator({
            **branches
            })


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        """Compute new quantities to be added to the NTuple."""

        ## Define accumulator
        output = self.accumulator.identity()

        pf_cands = utl.get_pf_cands(events, self.schema, self.pf_nano_aod_version)
        output["FatJet_ptD"] = processor.column_accumulator(np.array(ak.to_list(utl.calculate_ptD(pf_cands)), dtype=object))
        output["nFatJet"] = processor.column_accumulator(ak.to_numpy(utl.get_n_fat_jet(events, self.schema)))

        return output


    def postprocess(self, accumulator):
        return accumulator

