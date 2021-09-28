import sys

import awkward as ak

sys.path.append("../utilities/")
import nameUtilities as nameutl
import PtEtaPhiMLorentzVectorUtilities as vecutl


def get_met(events, met_collection_name="MET"):
    """Get met collection.

    Args:
        events (awkward.Array): the Events TTree opened with uproot.
        jet_type (str): AK4 or AK8

    Returns:
        awkward.Array: ak array with field pt, eta, phi, mass,
            that can be used as a coffea PtEtaPhiMLorentzVector.
    """

    met_branch = eval("events." + met_collection_name)

    met = ak.zip(
        {
            "pt"    : met_branch.pt,
            "mass"  : ak.zeros_like(met_branch.pt),
            "eta"   : ak.zeros_like(met_branch.pt),
            "phi"   : met_branch.phi,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    return met


def get_jets(events, jet_type):
    """Get jet collection for a given jet type.

    Args:
        events (awkward.Array): the Events TTree opened with uproot.
        jet_type (str): AK4 or AK8

    Returns:
        awkward.Array: ak array with field pt, eta, phi, mass,
            that can be used as a coffea PtEtaPhiMLorentzVector.
    """

    jet_collection_name = nameutl.jet_algo_name_to_jet_collection_name(jet_type)
    jet_branch = eval("events." + jet_collection_name)

    jets = ak.zip(
        {
            "pt"    : jet_branch.pt,
            "mass"  : jet_branch.mass,
            "eta"   : jet_branch.eta,
            "phi"   : jet_branch.phi,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    jets = ak.with_field(jets, vecutl.momentum(jets), where="momentum")

    return jets


def get_pf_cands(events, jet_type):
    """Get jet pf candidates collection for a given jet type.

    Args:
        events (awkward.Array): the Events TTree opened with uproot.
        jet_type (str): AK4 or AK8

    Returns:
        awkward.Array: ak array with field jetIdx, pt, eta, phi, mass,
            charge and pdgId that can be used as a coffea
            PtEtaPhiMLorentzVector.
    """

    jet_collection_name = nameutl.jet_algo_name_to_jet_collection_name(jet_type)
    jet_pf_cands_idx_info = eval("events." + jet_collection_name + "PFCands")
    jet_pf_cands = events.PFCands[jet_pf_cands_idx_info.pFCandsIdx]
    pf_cands = ak.zip(
        {
            "jetIdx": jet_pf_cands_idx_info.jetIdx,
            "pt"    : jet_pf_cands.pt,
            "mass"  : jet_pf_cands.mass,
            "eta"   : jet_pf_cands.eta,
            "phi"   : jet_pf_cands.phi,
            "charge": jet_pf_cands.charge,
            "pdgId" : jet_pf_cands.pdgId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    pf_cands = ak.with_field(pf_cands, vecutl.rapidity(pf_cands), where="rapidity")

    return pf_cands


def get_collection_size(events, collection):
    return ak.count(getattr(events, collection).pt, axis=1)


