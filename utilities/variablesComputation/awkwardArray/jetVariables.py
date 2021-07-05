import awkward as ak
import numpy as np
import sys

sys.path.append("../../../")
import awkwardArrayUtilities as akutl
import PtEtaPhiMLorentzVectorUtilities as vecutl
import physicsUtilities as phutl



def get_max_jet_idx(pf_cands):
    """Returns the maximum jet index from PF candidates.

    The PF candidates ak.Array passed in argument must have a jetIdx field.

    Args:
        pf_cands (ak.Array)

    Returns:
        int
    """

    return ak.max(pf_cands.jetIdx)


def calculate_ptD_1_jet(jet_idx, pf_cands, nan_value=-1):
    """Calculate ptD for all jets of a given index for all events.

    Args:
        jet_idx (int)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt and jetIdx.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).

    Returns:
        [awkward.Array]

    Examples:
	>>> pf_cands_idx = ak.Array([[0, 0, 1, 1, 1], [0, 0, 1, 1, 2, 2]])
	>>> pf_cands_pt = ak.Array([[280, 150, 150, 130, 70], [500, 300, 400, 300, 150, 100]])
	>>> pf_cands = ak.zip({"jetIdx": pf_cands_idx, "pt": pf_cands_pt})
	>>> jet_idx = 2
	>>> calculate_ptD_1_jet(jet_idx, pf_cands)
        [<Array [0.601, 0.714] type='2 * float64'>]
	>>> jet_idx = 2
	>>> calculate_ptD_1_jet(jet_idx, pf_cands)
        [<Array [-1, 0.721] type='2 * float64'>]
    """

    pf_cands_pt = pf_cands[(pf_cands.jetIdx == jet_idx)].pt

    numerator = np.sqrt(ak.sum(pf_cands_pt**2, axis=1))
    denominator = ak.sum(pf_cands_pt, axis=1)
    ptD = akutl.divide_ak_arrays(numerator, denominator, division_by_zero_value=nan_value)

    return [ptD]


def calculate_girth_1_jet(jet_idx, pf_cands, jets, njets, nan_value=-1):
    """Calculate girth for all jets of a given index for all events.

    Args:
        jet_idx (int)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt, eta, phi, mass and jetIdx,
            and with name PtEtaPhiMLorentzVector.
        jets (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with fields pt, eta, phi and mass and with name
            PtEtaPhiMLorentzVector.
        njets (awkward.Array):
            Ak array with one axis with number of jets in each event.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).

    Returns:
        [awkward.Array]

    Examples:
	>>> pf_cands_idx = ak.Array([[0, 0, 1, 1, 1], [0, 0, 1, 1, 2, 2]])
	>>> pf_cands_pt = ak.Array([[280, 150, 150, 130, 70], [500, 300, 400, 300, 150, 100]])
	>>> pf_cands_eta = ak.Array([[0.8, 0.87, -0.3, -0.2, -0.25], [0.5, 0.45, -1.4, -1.5, 0.1, 0.2]])
	>>> pf_cands_phi = ak.Array([[0.2, 0.28, 3.1, -3.1, 2.8], [-1.5, -1.2, 2.0, 2.5, -2.6, -2.9]])
	>>> pf_cands_mass = ak.Array([[50, 40, 20, 25, 23], [100, 50, 20, 45, 30, 15]])
	>>> pf_cands = ak.zip({"jetIdx": pf_cands_idx, "pt": pf_cands_pt, "eta": pf_cands_eta, "phi": pf_cands_phi, "mass": pf_cands_mass})
	>>> jets_pt = ak.Array([[430, 350], [800, 700, 250]])
	>>> jets_eta = ak.Array([[0.8, -0.28], [0.48, -1.45, 0.15]])
	>>> jets_phi = ak.Array([[0.25, 3.0], [-1.4, 2.2, -2.55]])
	>>> jets_mass = ak.Array([[20, 40], [15, 80, 50]])
	>>> jets = ak.zip({"pt": jets_pt, "eta": jets_eta, "phi": jets_phi, "mass": jets_mass})
	>>> njets = ak.Array([2, 3])
	>>> jet_idx = 1
	>>> calculate_girth_1_jet(jet_idx, pf_cands, jets, njets)
        [<Array [0.158, 0.248] type='2 * float64'>]
    """


    jet_pf_cands = pf_cands[(pf_cands.jetIdx == jet_idx)]
    jet = ak.mask(jets, njets > jet_idx)[:, jet_idx]

    # The following line introduces None in place of events with no jet with index jet_idx
    jet_broadcasted = akutl.broadcast(jet, jet_pf_cands)[0]
    delta_r = vecutl.delta_r(jet_pf_cands, jet_broadcasted)

    numerator = ak.sum(jet_pf_cands.pt * delta_r, axis=1)
    denominator = ak.sum(jet_pf_cands.pt, axis=1)
    girth = akutl.divide_ak_arrays(numerator, denominator, division_by_zero_value=nan_value)
    girth = ak.fill_none(girth, nan_value)

    return [girth]


def calculate_axes_1_jet(jet_idx, pf_cands, jets, njets, nan_value=-1):
    """Calculate axis major, minor and avg for all jets of a given index for all events.

    Args:
        jet_idx (int)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt, eta, phi, mass and jetIdx,
            and with name PtEtaPhiMLorentzVector.
        jets (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with fields pt, eta, phi and mass and with name
            PtEtaPhiMLorentzVector.
        njets (awkward.Array):
            Ak array with one axis with number of jets in each event.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).

    Returns:
        [awkward.Array, awkward.Array, awkward.Array]
    """

    jet_pf_cands = pf_cands[(pf_cands.jetIdx == jet_idx)]
    jet = ak.mask(jets, njets > jet_idx)[:, jet_idx]

    # The following line introduces None in place of events with no jet with index jet_idx
    jet_broadcasted = akutl.broadcast(jet, jet_pf_cands)[0]
    jet_pf_cands_delta_phi = vecutl.delta_phi(jet_pf_cands, jet_broadcasted)
    jet_pf_cands_delta_eta = vecutl.delta_eta(jet_pf_cands, jet_broadcasted)
    jet_pf_cands_pt2 = jet_pf_cands.pt**2

    # lambda1 = Tr(M)/2 + sqrt( (Tr(M)/2)**2 - Det(M) ) = u + v
    # lambda2 =               -                         = u - v
    # Tr(M)/2 = (M11 + M22)/2 
    # Det(M) = M11*M22 - M12**2
    sum_jet_pf_cands_pt2 = ak.sum(jet_pf_cands_pt2, axis=1)
    m11 = ak.sum(jet_pf_cands_pt2 * jet_pf_cands_delta_eta**2, axis=1)
    m22 = ak.sum(jet_pf_cands_pt2 * jet_pf_cands_delta_phi**2, axis=1)
    m12 = -ak.sum(jet_pf_cands_pt2 * jet_pf_cands_delta_eta * jet_pf_cands_delta_phi, axis=1)
    tr_m = m11 + m22
    det_m = m11*m22 - m12**2
    u = tr_m/2
    v = np.sqrt( u**2 - det_m )
    lambda1 = u + v
    lambda2 = u - v

    axis_major = np.sqrt( lambda1 / sum_jet_pf_cands_pt2 )
    axis_minor = np.sqrt( lambda2 / sum_jet_pf_cands_pt2 )
    axis_avg   = np.sqrt( axis_major**2 + axis_minor**2 )

    axis_major = ak.fill_none(ak.nan_to_num(axis_major, nan=nan_value), nan_value)
    axis_minor = ak.fill_none(ak.nan_to_num(axis_minor, nan=nan_value), nan_value)
    axis_avg   = ak.fill_none(ak.nan_to_num(axis_avg  , nan=nan_value), nan_value)
    
    return [axis_major, axis_minor, axis_avg]


def calculate_HEF_1_jet(charged, jet_idx, pf_cands, nan_value=-1):
    """Calculate energy fraction for charged or neutral hadrons.

    Args:
        charged (bool)
        jet_idx (int)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields charge, energy and jetIdx.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).

    Returns:
        [awkward.Array]
    """

    jet_pf_cands = pf_cands[(pf_cands.jetIdx == jet_idx)]
    hadrons_jet_pf_cands = jet_pf_cands[akutl.is_in(jet_pf_cands.pdgId, phutl.pdg_id("hadrons"))]

    if charged:
        selected_hadrons_jet_pf_cands = hadrons_jet_pf_cands[(hadrons_jet_pf_cands.charge != 0)]
    else:
        selected_hadrons_jet_pf_cands = hadrons_jet_pf_cands[(hadrons_jet_pf_cands.charge == 0)]

    numerator = ak.sum(selected_hadrons_jet_pf_cands.energy, axis=1)
    denominator = ak.sum(jet_pf_cands.energy, axis=1)
    energy_fraction = akutl.divide_ak_arrays(numerator, denominator, division_by_zero_value=nan_value)

    return [energy_fraction]


def calculate_chHEF_1_jet(jet_idx, pf_cands, nan_value=-1):
    """Return energy fraction for charged hadrons.

    See doctring of calculate_HEF_1_jet.
    """

    return calculate_HEF_1_jet(True, jet_idx, pf_cands, nan_value=-1)


def calculate_neHEF_1_jet(jet_idx, pf_cands, nan_value=-1):
    """Return energy fraction for neutral hadrons.

    See doctring of calculate_HEF_1_jet.
    """

    return calculate_HEF_1_jet(False, jet_idx, pf_cands, nan_value=-1)


def calculate_jet_variable(variable, pf_cands, jets, njets, nan_value=-1):
    """Return the desired calculated jet variable for all jets and all events.

    Args:
        variable (str)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with required fields.
        jets (awkward.Array or None):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with required field.
        njets (awkward.Array or None):
            Ak array with one axis with number of jets in each event.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).

    Returns:
        [awkward.Array]
    """

    max_idx = get_max_jet_idx(pf_cands)

    for jet_idx in range(max_idx+1):
        if variable == "ptD":
            ak_arrays = calculate_ptD_1_jet(jet_idx, pf_cands, nan_value)
        elif variable == "girth":
            ak_arrays = calculate_girth_1_jet(jet_idx, pf_cands, jets, njets, nan_value)
        elif variable == "axes":
            ak_arrays = calculate_axes_1_jet(jet_idx, pf_cands, jets, njets, nan_value)
        elif variable == "chHEF":
            ak_arrays = calculate_chHEF_1_jet(jet_idx, pf_cands, nan_value)
        elif variable == "neHEF":
            ak_arrays = calculate_neHEF_1_jet(jet_idx, pf_cands, nan_value)
        
        if jet_idx == 0:
            ak_array_list = len(ak_arrays) * [ ak.Array([]) ]
        for idx, ak_array in enumerate(ak_arrays):
            ak_array_list[idx] = ak.concatenate((ak_array_list[idx], [ak_array]), axis=0)

    # Need to transpose / swap axes of the ak array to have events in the first axis
    # instead of jets. Maybe there is something smarter than that.
    for idx, ak_array in enumerate(ak_array_list):
        ak_array = akutl.swap_axes(ak_array)
        ak_array = ak_array[ak_array != nan_value][0]
        ak_array_list[idx] = ak_array

    return ak_array_list


def calculate_ptD(pf_cands):
    """Return ptD for each jet for all events.

    See doctring of calculate_ptD_1_jet.

    Args:
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt and jetIdx.
    
    Returns:
         [awkward.Array]

    Examples:
	>>> pf_cands_idx = ak.Array([[0, 0, 1, 1, 1], [0, 0, 1, 1, 2, 2]])
	>>> pf_cands_pt = ak.Array([[280, 150, 150, 130, 70], [500, 300, 400, 300, 150, 100]])
        >>> pf_cands = ak.zip({"jetIdx": pf_cands_idx, "pt": pf_cands_pt})
        >>> ak.to_list(calculate_ptD(pf_cands)[0])
        [[0.739, 0.601], [0.729, 0.714, 0.721]]
    """

    return calculate_jet_variable("ptD", pf_cands, None, None)


def calculate_girth(pf_cands, jets, njets):
    """Return girth for each jet for all events.
    
    See doctring of calculate_girth_1_jet.
    """

    return calculate_jet_variable("girth", pf_cands, jets, njets)


def calculate_axes(pf_cands, jets, njets):
    """Return axis major, axis minor and axis avg for each jet for all events.

    See doctring of calculate_axes_1_jet.
    """

    return calculate_jet_variable("axes", pf_cands, jets, njets)


def calculate_chHEF(pf_cands, jets, njets):
    """Return charged hadron energy fraction for each jet for all events.
    
    See doctring of calculate_HEF_1_jet.
    """

    return calculate_jet_variable("chHEF", pf_cands, jets, njets)


def calculate_neHEF(pf_cands, jets, njets):
    """Return neutral hadron energy fraction for each jet for all events.

    See doctring of calculate_HEF_1_jet.
    """

    return calculate_jet_variable("neHEF", pf_cands, jets, njets)


