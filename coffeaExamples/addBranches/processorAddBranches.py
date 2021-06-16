from coffea import processor
import awkward as ak
import numpy as np
import sys


sys.path.append("../../pythonUtils/")
import awkwardArrayUtils as akutl


## Define some short-hands for column and value accumulator
def column_accumulator(type_):
    return processor.column_accumulator(np.array([], dtype=type_))


def calculate_ptD(pf_cands):
    """Calculate ptD for each jet in all events.

    Calculate ptD for each jet in all events using optimized columnar
    operations. However, cannot avoid copy of the array at the end to
    transpose / swap axes of the ak array.

    Args:
        pf_cands (awkward.Array):
            Ak array where axis 0 is events axis and axis 1 has pt and jetIdx
            fields.
    
    Returns:
         awkward.Array

    Examples:
        >>> pf_cands_idx = ak.Array([[0,0,1,1,1,2,2,2], [0,0,1,2,2,3]])
        >>> jet_cands_pt = ak.Array([[180,150,150,130,70,180,100,20], [500, 300, 400, 300, 250, 100]])
        >>> pf_cands = ak.zip({"jetIdx": jet_cands_idx, "pt": jet_cands_pt})
        >>> calculate_ptD(pf_cands)
        [[0.71, 0.601, 0.69], [0.729, 1, 0.71, 1]]
    """

    max_idx = ak.max(pf_cands.jetIdx)

    for jet_idx in range(max_idx+1):
        pf_cands_pt = pf_cands[(pf_cands.jetIdx == jet_idx)].pt

        ptD_denominator = ak.sum(pf_cands_pt, axis=1)
        ptD_numerator = np.sqrt(ak.sum(pf_cands_pt**2, axis=1))
        ptD_tmp = akutl.divide_ak_arrays(ptD_numerator, ptD_denominator, division_by_zero_value=-1)

        if jet_idx == 0:
            ptD = [ptD_tmp]
        else:
            ptD = ak.concatenate((ptD, [ptD_tmp]), axis=0)

    # Need to transpose / swap axes of the ak array to have events in the first axis
    # instead of jets. Maybe there is something smarter than that.
    np_array = ak.to_numpy(ptD).T
    ak_array = ak.Array([np_array])     # Can't avoid that copy in memory!!!
    ptD = ak_array[ak_array != -1][0]

    return ptD


def get_pf_cands(events, schema, pf_nano_aod_version):

    if schema == "Base":
        if pf_nano_aod_version == "102X":
            pf_cands = ak.zip({
                "jetIdx": events.FatJetPFCands_jetIdx,
                "pt": events.FatJetPFCands_pt,
            })
        if pf_nano_aod_version == "106X":
            pf_cands = ak.zip({
                "jetIdx": events.JetPFCandsAK8_jetIdx,
                "pt": events.JetPFCands_pt[events.JetPFCandsAK8_candIdx]
            })
    elif schema == "PFNanoAOD":
        if pf_nano_aod_version == "102X":
            pf_cands = ak.zip({
                "jetIdx": events.FatJetPFCands.jetIdx,
                "pt": events.FatJetPFCands.pt,
            })
        if pf_nano_aod_version == "106X":
            pf_cands = ak.zip({
                "jetIdx": events.JetPFCandsAK8.jetIdx,
                "pt": events.JetPFCands[events.JetPFCandsAK8.candIdx].pt,
            })

    return pf_cands


def get_n_fat_jet(events, schema):
    if schema == "Base":
        return events.nFatJet
    elif schema == "PFNanoAOD":
        return ak.count(events.FatJet.pt, axis=1)


class BranchesMaker(processor.ProcessorABC):
    """Example for calculating ptD."""


    def __init__(self, schema, pf_nano_aod_version):

        self.schema = schema
        self.pf_nano_aod_version = pf_nano_aod_version

        branches = { 
            "FatJet_ptD"      : column_accumulator(np.float64),
            "nFatJet"         : column_accumulator(np.float64),
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

        pf_cands = get_pf_cands(events, self.schema, self.pf_nano_aod_version)
        output["FatJet_ptD"] = processor.column_accumulator(np.array(ak.to_list(calculate_ptD(pf_cands)), dtype=object))
        output["nFatJet"] = processor.column_accumulator(ak.to_numpy(get_n_fat_jet(events, self.schema)))

        return output


    def postprocess(self, accumulator):
        return accumulator

