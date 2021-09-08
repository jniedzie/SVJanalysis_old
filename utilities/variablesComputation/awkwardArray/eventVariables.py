import awkward as ak
import numpy as np
import sys

sys.path.append("../../../")
import awkwardArrayUtilities as akutl
import PtEtaPhiMLorentzVectorUtilities as vecutl


def calculate_number_of_jets(jets):
    """Calculate number of jets in all events.

    Args:
        jets (awkward.Array): Jagged ak array where axis 0 is the event axis
            and axis 1 is the jet axis.

    Returns:
        awkward.Array
    """

    return ak.num(jets, axis=1)


def calculate_razor_MR(jets, njets=None, nan_value=0):
    """Calculate M_R Razor variable for all events using 2 leading jets.

    Args:
        jets (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with fields pt, rapidity, phi and mass and with name
            PtEtaPhiMLorentzVector.
        njets (awkward.Array, optional, default=None):
            Ak array with one axis with number of jets in each event.
            If None, will be computed from jets.
        nan_value (float, optional, default=0):
            Value to use when the event has less than 2 jets.

    Returns:
        awkward.Array
    """

    if njets is None:
        njets = calculate_number_of_jets(jets)

    jets = ak.mask(jets, njets>=2)

    jets0 = jets[:, 0]
    jets1 = jets[:, 1]
    mr = np.sqrt( (jets0.momentum + jets1.momentum) ** 2 - (jets0.pz + jets1.pz) ** 2 )
    mr = ak.fill_none(mr, nan_value)

    return mr


def calculate_razor_MRT(jets, met, njets=None, nan_value=0):
    """Calculate M^R_T Razor variable for all events using 2 leading jets.

    Args:
        jets (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with fields pt, eta, phi and mass and with name
            PtEtaPhiMLorentzVector.
        met (awkward.Array):
            Missing transverse energy. 1D ak array with fields pt, eta, phi,
            mass and with name PtEtaPhiMLorentzVector.
        njets (awkward.Array, optional, default=None):
            Ak array with one axis with number of jets in each event.
            If None, will be computed from jets.
        nan_value (float, optional, default=0):
            Value to use when the event has less than 2 jets.

    Returns:
        awkward.Array
    """

    if njets is None:
        njets = calculate_number_of_jets(jets)

    jets = ak.mask(jets, njets>=2)
    met = ak.mask(met, njets>=2)

    jet0 = jets[:, 0]
    jet1 = jets[:, 1]
    mrt = np.sqrt( met.pt * (jet0.pt + jet1.pt) - met.dot(jet0 + jet1) ) / np.sqrt(2)
    mrt = ak.fill_none(mrt, nan_value)

    return mrt


def calculate_razor_R(mr, mrt, nan_value=0):
    """Calculate R Razor variable for all events from M_R and M^R_T variables.

    Args:
        mr (awkward.Array)
        mrt (awkward.Array)
    
    Returns:
        awkward.Array
    """

    return akutl.divide_ak_arrays(mrt, mr, division_by_zero_value=nan_value)


def calculate_MT(jets, met, jet_indices=[0, 1], njets=None, nan_value=0):
    """Calculate transverse mass variable for all events using 2 leading jets.

    Args:
        jets (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with fields pt, rapidity, phi and mass and with name
            PtEtaPhiMLorentzVector.
        njets (awkward.Array, optional, default=None):
            Ak array with one axis with number of jets in each event.
            If None, will be computed from jets.
        nan_value (float, optional, default=0):
            Value to use when the event has less than 2 jets.

    Returns:
        awkward.Array
    """

    if njets is None:
        njets = calculate_number_of_jets(jets)

    max_index = max(jet_indices)

    jets = ak.mask(jets, njets>=max_index+1)
    met = ak.mask(met, njets>=max_index+1)

    jet0 = jets[:, jet_indices[0]]
    jet1 = jets[:, jet_indices[1]]
    dijet = jet0 + jet1
    mt = np.sqrt( dijet.mass**2 + 2 * ( np.sqrt(dijet.mass**2 + dijet.pt**2) * met.pt - met.dot(dijet) ) )
    mt = ak.fill_none(mt, nan_value)

    return mt
