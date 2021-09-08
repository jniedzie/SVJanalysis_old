import energyflow as ef
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


def calculate_sum_pfcands_pt_1_jet(jet_idx, pf_cands):
    """Returns the sum of PF candidates pT for for all jets of a given index for all events.

    Returns 0 for events having no jet with jet index jet_idx.

    Args:
        jet_idx (int)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt and jetIdx.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).
        
    Returns:
        awkward.Array
    """

    jet_pf_cands = pf_cands[(pf_cands.jetIdx == jet_idx)]

    return ak.sum(jet_pf_cands.pt, axis=1)


def calculate_generalized_angularity_1_jet(jet_idx, pf_cands, jets, njets, jet_radius, beta, kappa, sum_pfcands_pt=None, nan_value=-1):
    """Calculate generalized angularity for all jets of a given index for all events.

    Args:
        jet_idx (int)
        pf_cands (awkward.Array or None):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt, eta, phi, mass and jetIdx,
            and with name PtEtaPhiMLorentzVector.
            Can be None if kappa=0.
        jets (awkward.Array or None):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with fields pt, eta, phi and mass and with name
            PtEtaPhiMLorentzVector.
            Can be None if beta=0.
        njets (awkward.Array):
            Ak array with one axis with number of jets in each event.
        beta (float): Angular distance exponent.
        kappa (float): Transverse momentum fraction exponent.
        sum_pfcands_pt (awkward.Array, optional, default=None):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis.
            Sum of jet constituents pT for all jets in all events.
            If None, will be computed from pf_cands.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).

    Returns:
        awkward.Array

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
        >>> # Girth with jet radius R=1
        >>> calculate_generalized_angularity_1_jet(jet_idx, pf_cands, jets, njets, 1, 1, 1)
        [<Array [0.158, 0.248] type='2 * float64'>]
        >>> # ptD squared with more jets than jet_idx+1
        >>> calculate_generalized_anglarity_1_jet(jet_idx, pf_cands, None, None, 1, 0, 2)
        [<Array [0.601, 0.714] type='2 * float64'>]
        >>> # ptD squared with fewer jets than jet_idx+1, note the -1 
        >>> jet_idx = 2
        >>> calculate_generalized_anglarity_1_jet(jet_idx, pf_cands)
        [<Array [-1, 0.721] type='2 * float64'>]
    """

    jet_pf_cands = pf_cands[(pf_cands.jetIdx == jet_idx)]

    if beta == 0 and kappa == 0:
        generalized_angularity = ak.num(jet_pf_cands.pt, axis=1)

    else:
        if beta == 0:
            angular_term = 1.
        else:
            # The following line introduces None in place of events with no jet with index jet_idx
            jet = ak.mask(jets, njets > jet_idx)[:, jet_idx]
            jet_broadcasted = akutl.broadcast(jet, jet_pf_cands)[0]
            delta_r = vecutl.delta_r(jet_pf_cands, jet_broadcasted)
            angular_term = delta_r**beta
    
        if kappa == 0:
            sum_pfcands_pt = 1.
            momentum_term = 1.
        else:
            if sum_pfcands_pt is not None:
                sum_pfcands_pt = sum_pfcands_pt[:, jet_idx]
            else:
                sum_pfcands_pt = calculate_sum_pfcands_pt_1_jet(jet_idx, pf_cands)
            momentum_term = jet_pf_cands.pt**kappa

        numerator = ak.sum(momentum_term * angular_term, axis=1)
        denominator = sum_pfcands_pt**kappa * jet_radius**beta
        generalized_angularity = akutl.divide_ak_arrays(numerator, denominator, division_by_zero_value=nan_value)
        generalized_angularity = ak.fill_none(generalized_angularity, nan_value)
    
    return generalized_angularity


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
        tuple: (awkward.Array, awkward.Array, awkward.Array)
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
    
    return axis_major, axis_minor, axis_avg


def calculate_efps_1_jet(jet_idx, pf_cands, degree):
    """Calculate EFPs for all jets of a given index for all events.

    Calculates Energy Flow Polynomials for graphs with degree (= number of
    edges) at least "degree" for all jets of a given index for all events.
    For jets with less constituents than required, the corresponding EFP is
    0. This is used to return 0 for events without jet with index "jet_idx".

    Args:
        jet_idx (int)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt, rapidity, phi, mass and jetIdx,
            and with name PtEtaPhiMLorentzVector.

    Returns:
        numpy.ndarray
    """

    jet_pf_cands_ptyphim = akutl.ak_to_ptyphim_four_vectors(pf_cands, jet_idx)
    
    efpset = ef.EFPSet('d<='+str(degree), measure='hadr', beta=1, normed=True, verbose=False, coords="ptyphim")
    efps = efpset.batch_compute(jet_pf_cands_ptyphim)

    return tuple(efps.T)
    

def calculate_ecf_e2_1_jet(jet_idx, pf_cands, beta, sum_pfcands_pt=None, nan_value=-1):
    """Calculate 2-point ECF for all jets of a given index for all events.

    Args:
        jet_idx (int)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt, rapidity, phi, mass and jetIdx,
            and with name PtEtaPhiMLorentzVector.
        beta (float):
            angular exponent
        sum_pfcands_pt (awkward.Array, optional, default=None):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis.
            Sum of jet constituents pT for all jets in all events.
            If None, will be computed from pf_cands.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).

    Returns:
        awkward.Array
    """

    jet_pf_cands = pf_cands[(pf_cands.jetIdx == jet_idx)]
    combinations = ak.combinations(jet_pf_cands, 2, axis=1, fields=["i", "j"])
    delta_r_ij = vecutl.delta_r(combinations.i, combinations.j, use_rapidity=True)
    pti = combinations.i.pt
    ptj = combinations.j.pt
    ecf = ak.sum(pti*ptj*delta_r_ij**beta, axis=1)

    if sum_pfcands_pt is not None:
        ecf = akutl.divide_ak_arrays(ecf, sum_pfcands_pt[:, jet_idx]**2, division_by_zero_value=nan_value)
    else:
        ecf = akutl.divide_ak_arrays(ecf, calculate_sum_pfcands_pt_1_jet(jet_idx, pf_cands)**2, division_by_zero_value=nan_value)

    return ecf


def calculate_ecfs_e3_1_jet(jet_idx, pf_cands, beta, calculate_ecfgs=True, sum_pfcands_pt=None, nan_value=-1):
    """Calculate 3-point generalized ECFs for all jets of a given index for all events.

    Args:
        jet_idx (int)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt, rapidity, phi, mass and jetIdx,
            and with name PtEtaPhiMLorentzVector.
        beta (float):
            angular exponent
        calculate_ecfgs (bool, optional, default=True):
            False to calculate only 3e3, True for all generalized ECFs.
        sum_pfcands_pt (awkward.Array, optional, default=None):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis.
            Sum of jet constituents pT for all jets in all events.
            If None, will be computed from pf_cands.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).

    Returns:
        awkward.Array
    """

    jet_pf_cands = pf_cands[(pf_cands.jetIdx == jet_idx)]
    combinations = ak.combinations(jet_pf_cands, 3, axis=1, fields=["i", "j", "k"])
    delta_r_ij = vecutl.delta_r(combinations.i, combinations.j, use_rapidity=True)
    delta_r_ik = vecutl.delta_r(combinations.i, combinations.k, use_rapidity=True)
    delta_r_jk = vecutl.delta_r(combinations.j, combinations.k, use_rapidity=True)
    pti = combinations.i.pt
    ptj = combinations.j.pt
    ptk = combinations.k.pt
    pt_product = pti * ptj * ptk

    if calculate_ecfgs:
        sorted_array = ak.sort([delta_r_ij, delta_r_ik, delta_r_jk], axis=0)
        ecfgs_base = [pt_product] + 3 * [None]
        ecfgs = 3 * [None]
        for idx in range(1, 4):
            ecfgs_base[idx] = ecfgs_base[idx-1] * sorted_array[idx-1]**beta
            ecfgs[idx-1] = ak.sum(ecfgs_base[idx], axis=1)
    else:
        ecf = ak.sum(pt_product * (delta_r_ij*delta_r_ik*delta_r_jk)**beta, axis=1)

    ## Normalize by the sum of constituents pT
    if sum_pfcands_pt is not None:
        norm = sum_pfcands_pt[:, jet_idx]**3
    else:
        norm = calculate_sum_pfcands_pt_1_jet(jet_idx, pf_cands)**3

    if calculate_ecfgs:
        for idx in range(3):
            ecfgs[idx] = akutl.divide_ak_arrays(ecfgs[idx], norm, division_by_zero_value=nan_value)
    else:
        ecf = akutl.divide_ak_arrays(ecf, norm, division_by_zero_value=nan_value)

    if calculate_ecfgs:
        return tuple(ecfgs)
    else:
        return ecf


def calculate_ecfs_e4_1_jet(jet_idx, pf_cands, beta, calculate_ecfgs=True, sum_pfcands_pt=None, nan_value=-1):
    """Calculate 4-point generalized ECFs for all jets of a given index for all events.

    Args:
        jet_idx (int)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt, rapidity, phi, mass and jetIdx,
            and with name PtEtaPhiMLorentzVector.
        beta (float):
            angular exponent
        calculate_ecfgs (bool, optional, default=True):
            False to calculate only 6e4, True for all generalized ECFs.
        sum_pfcands_pt (awkward.Array, optional, default=None):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis.
            Sum of jet constituents pT for all jets in all events.
            If None, will be computed from pf_cands.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).

    Returns:
        awkward.Array
    """

    jet_pf_cands = pf_cands[(pf_cands.jetIdx == jet_idx)]
    combinations = ak.combinations(jet_pf_cands, 4, axis=1, fields=["i", "j", "k", "l"])
    delta_r_ij = vecutl.delta_r(combinations.i, combinations.j, use_rapidity=True)
    delta_r_ik = vecutl.delta_r(combinations.i, combinations.k, use_rapidity=True)
    delta_r_il = vecutl.delta_r(combinations.i, combinations.l, use_rapidity=True)
    delta_r_jk = vecutl.delta_r(combinations.j, combinations.k, use_rapidity=True)
    delta_r_jl = vecutl.delta_r(combinations.j, combinations.l, use_rapidity=True)
    delta_r_kl = vecutl.delta_r(combinations.k, combinations.l, use_rapidity=True)
    pti = combinations.i.pt
    ptj = combinations.j.pt
    ptk = combinations.k.pt
    ptl = combinations.l.pt
    pt_product = pti * ptj * ptk * ptl

    if calculate_ecfgs:
        sorted_array = ak.sort([delta_r_ij, delta_r_ik, delta_r_il, delta_r_jk, delta_r_jl, delta_r_kl], axis=0)

        ecfgs_base = [pt_product] + 6 * [None]
        ecfgs = 6 * [None]
        for idx in range(1, 7):
            ecfgs_base[idx] = ecfgs_base[idx-1] * sorted_array[idx-1]**beta
            ecfgs[idx-1] = ak.sum(ecfgs_base[idx], axis=1)
    else:
        ecf = ak.sum(pt_product * (delta_r_ij*delta_r_ik*delta_r_il*delta_r_jk*delta_r_jl*delta_r_kl)**beta, axis=1)

    ## Normalize by the sum of constituents pT
    if sum_pfcands_pt is not None:
        norm = sum_pfcands_pt[:, jet_idx]**4
    else:
        norm = calculate_sum_pfcands_pt_1_jet(jet_idx, pf_cands)**4

    if calculate_ecfgs:
        for idx in range(6):
            ecfgs[idx] = akutl.divide_ak_arrays(ecfgs[idx], norm, division_by_zero_value=nan_value)
    else:
        ecf = akutl.divide_ak_arrays(ecf, norm, division_by_zero_value=nan_value)

    if calculate_ecfgs:
        return tuple(ecfgs)
    else:
        return ecf


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
        awkward.Array
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

    return energy_fraction


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


def calculate_jet_variable(variable, pf_cands, jets, njets=None, sum_pfcands_pt=None, nan_value=-1, jagged=True, params={}):
    """Return the desired calculated jet variable for all jets and all events.

    Args:
        variable (str)
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with required fields.
        jets (awkward.Array or None):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with required field.
        njets (awkward.Array, optional, default=None):
            Ak array with one axis with number of jets in each event.
            If None, will be computed from jets.
        sum_pfcands_pt (awkward.Array, optional, default=None):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis.
            Sum of jet constituents pT for all jets in all events.
            If None, will be computed from jets in *_1_jet functions.
        nan_value (float, optional, default=-1):
            Value to use when the event has less than jet_idx+1 jet(s).
        jagged (bool, optional, default=True):
            Return a jagged (True) or rectangular (False) array.
        params (dict, optional):
            Parameters on which depend the function to compute the variable.

    Returns:
        awkward.Array or tuple[awkward.Array]
    """

    if njets is None and jets is not None:
        njets = ak.num(jets, axis=1)

    max_idx = get_max_jet_idx(pf_cands)

    for jet_idx in range(max_idx+1):
        if variable == "sum_pfcands_pt":
            ak_arrays = calculate_sum_pfcands_pt_1_jet(jet_idx, pf_cands)
        elif variable == "generalized_angularity":
            ak_arrays = calculate_generalized_angularity_1_jet(jet_idx, pf_cands, jets, njets, params["jet_radius"], params["beta"], params["kappa"], sum_pfcands_pt=sum_pfcands_pt, nan_value=nan_value)
        elif variable == "ptD":
            ak_arrays = calculate_ptD_1_jet(jet_idx, pf_cands, nan_value=nan_value)
        elif variable == "girth":
            ak_arrays = calculate_girth_1_jet(jet_idx, pf_cands, jets, njets, nan_value=nan_value)
        elif variable == "axes":
            ak_arrays = calculate_axes_1_jet(jet_idx, pf_cands, jets, njets, nan_value=nan_value)
        elif variable == "efps":
            ak_arrays = calculate_efps_1_jet(jet_idx, pf_cands, params["degree"])
        elif variable == "ecf_e2":
            ak_arrays = calculate_ecf_e2_1_jet(jet_idx, pf_cands, params["beta"], sum_pfcands_pt=sum_pfcands_pt, nan_value=nan_value)
        elif variable == "ecfs_e3":
            ak_arrays = calculate_ecfs_e3_1_jet(jet_idx, pf_cands, params["beta"], calculate_ecfgs=params["calculate_ecfgs"], sum_pfcands_pt=sum_pfcands_pt, nan_value=nan_value)
        elif variable == "ecfs_e4":
            ak_arrays = calculate_ecfs_e4_1_jet(jet_idx, pf_cands, params["beta"], calculate_ecfgs=params["calculate_ecfgs"], sum_pfcands_pt=sum_pfcands_pt, nan_value=nan_value)
        elif variable == "chHEF":
            ak_arrays = calculate_chHEF_1_jet(jet_idx, pf_cands, nan_value=nan_value)
        elif variable == "neHEF":
            ak_arrays = calculate_neHEF_1_jet(jet_idx, pf_cands, nan_value=nan_value)
        else:
            print("ERROR: Unknown variable %s. Exit!" %variable)
            sys.exit()

        if not isinstance(ak_arrays, tuple):
            is_tuple = False
            ak_arrays = (ak_arrays, )
        else:
            is_tuple = True

        if jet_idx == 0:
            ak_array_list = len(ak_arrays) * [ ak.Array([]) ]
        for idx, ak_array in enumerate(ak_arrays):
            ak_array_list[idx] = ak.concatenate((ak_array_list[idx], [ak_array]), axis=0)

    # Need to transpose / swap axes of the ak array to have events in the first axis
    # instead of jets. Maybe there is something smarter than that.
    for idx, ak_array in enumerate(ak_array_list):
        ak_array = akutl.swap_axes(ak_array)
        if jagged:
            ak_array = ak_array[ak_array != nan_value][0]
        else:
            ak_array = ak_array[0]
        ak_array_list[idx] = ak_array

    if not is_tuple:
        return ak_array_list[0]
    else:
        return ak_array_list


def calculate_sum_pfcands_pt(pf_cands, jagged=True):
    """Return ECF1 for each jet for all events.

    See doctring of calculate_sum_pfcands_pt.
    """

    return calculate_jet_variable("sum_pfcands_pt", pf_cands, None, None, nan_value=0., jagged=jagged)


def calculate_generalized_angularity(pf_cands, jets, jet_radius, beta, kappa, njets=None, sum_pfcands_pt=None, nan_value=-1, jagged=True):
    """Return generalized angularity with parameters beta and kappa for each jet for all events.

    See docstring of calculate_generalized_angularity_1_jet.
    """

    return calculate_jet_variable("generalized_angularity", pf_cands, jets, njets=njets, \
                                 sum_pfcands_pt=sum_pfcands_pt, jagged=jagged, nan_value=nan_value, \
                                 params={"jet_radius": jet_radius, "beta": beta, "kappa": kappa})


def calculate_ptD(pf_cands, sum_pfcands_pt=None, jagged=True):
    """Return ptD for each jet for all events.

    See docstring of calculate_generalized_angularity_1_jet.

    Args:
        pf_cands (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the PF
            candidates axis with fields pt and jetIdx.
    
    Returns:
         awkward.Array

    Examples:
        >>> pf_cands_idx = ak.Array([[0, 0, 1, 1, 1], [0, 0, 1, 1, 2, 2]])
        >>> pf_cands_pt = ak.Array([[280, 150, 150, 130, 70], [500, 300, 400, 300, 150, 100]])
        >>> pf_cands = ak.zip({"jetIdx": pf_cands_idx, "pt": pf_cands_pt})
        >>> ak.to_list(calculate_ptD(pf_cands)[0])
        [[0.739, 0.601], [0.729, 0.714, 0.721]]
    """

    ptD2 = calculate_generalized_angularity(pf_cands, None, 1, 0, 2, sum_pfcands_pt=sum_pfcands_pt, jagged=jagged)
    ptD = np.sqrt(ptD2)
    return ptD


def calculate_lha(pf_cands, jets, jet_radius, njets=None, sum_pfcands_pt=None, jagged=True):
    """Return Les Houches Angularity for each jet for all events.
    
    See docstring of calculate_generalized_angularity_1_jet.
    """

    return calculate_generalized_angularity(pf_cands, jets, jet_radius, 0.5, 1, njets=njets, sum_pfcands_pt=sum_pfcands_pt, jagged=jagged)


def calculate_girth(pf_cands, jets, jet_radius, njets=None, sum_pfcands_pt=None, jagged=True):
    """Return girth for each jet for all events.
    
    See docstring of calculate_generalized_angularity_1_jet.
    """

    return calculate_generalized_angularity(pf_cands, jets, jet_radius, 1, 1, njets=njets, sum_pfcands_pt=sum_pfcands_pt, jagged=jagged)


def calculate_thrust(pf_cands, jets, jet_radius, njets=None, sum_pfcands_pt=None, jagged=True):
    """Return thrust for each jet for all events.
    
    See docstring of calculate_generalized_angularity_1_jet.
    """

    return calculate_generalized_angularity(pf_cands, jets, jet_radius, 2, 1, njets=njets, sum_pfcands_pt=sum_pfcands_pt, jagged=jagged)


def calculate_multiplicity(pf_cands, nan_value=0, jagged=True):
    """Return multiplicity for each jet for all events.
    
    See docstring of calculate_generalized_angularity_1_jet.
    """

    return calculate_generalized_angularity(pf_cands, None, 1., 0, 0, nan_value=0, jagged=jagged)


def calculate_axes(pf_cands, jets, njets=None, jagged=True):
    """Return axis major, axis minor and axis avg for each jet for all events.

    See doctring of calculate_axes_1_jet.
    """

    return calculate_jet_variable("axes", pf_cands, jets, njets=njets, jagged=jagged)


def calculate_efps(pf_cands, degree, jagged=True):
    """Return EFPs up to a certain degree for each jet for all events.

    See doctring of calculate_efps_1_jet.
    """

    return calculate_jet_variable("efps", pf_cands, None, nan_value=0, jagged=jagged, params={"degree": degree})


def calculate_ecf_e2(pf_cands, beta, sum_pfcands_pt=None, jagged=True):
    """Return ECF2 for each jet for all events.

    See doctring of calculate_ecf_e2_1_jet.
    """

    return calculate_jet_variable("ecf_e2", pf_cands, None, sum_pfcands_pt=sum_pfcands_pt, nan_value=-1, jagged=jagged, params={"beta": beta})


def calculate_ecfs_e3(pf_cands, beta, sum_pfcands_pt=None, jagged=True, calculate_ecfgs=True):
    """Return ECF3 for each jet for all events.

    See doctring of calculate_ecf_e3_1_jet.
    """

    return calculate_jet_variable("ecfs_e3", pf_cands, None, sum_pfcands_pt=sum_pfcands_pt, nan_value=-1, jagged=jagged, params={"beta": beta, "calculate_ecfgs": calculate_ecfgs})


def calculate_ecfs_e4(pf_cands, beta, sum_pfcands_pt=None, jagged=True, calculate_ecfgs=True):
    """Return ECF4 for each jet for all events.

    See doctring of calculate_ecf_e4_1_jet.
    """

    return calculate_jet_variable("ecfs_e4", pf_cands, None, sum_pfcands_pt=sum_pfcands_pt, nan_value=-1, jagged=jagged, params={"beta": beta, "calculate_ecfgs": calculate_ecfgs})


def calculate_chHEF(pf_cands, jets, njets=None, jagged=True):
    """Return charged hadron energy fraction for each jet for all events.
    
    See doctring of calculate_HEF_1_jet.
    """

    return calculate_jet_variable("chHEF", pf_cands, jets, njets=njets, jagged=jagged)


def calculate_neHEF(pf_cands, jets, njets=None, jagged=True):
    """Return neutral hadron energy fraction for each jet for all events.

    See doctring of calculate_HEF_1_jet.
    """

    return calculate_jet_variable("neHEF", pf_cands, jets, njets=njets, jagged=jagged)



def calculate_ecf_c(ei, ej, ek):
    """Calculate C-series ECFs.
    
    ei (resp. ej, ek) is the i-point (resp. j-point, k-point) normalized ECF.
    To calculate C-series ECFs: j=i+1 and k=i+2.

    Args:
        ei (awkward.Array): i-point ECF
        ej (awkward.Array): j-point ECF
        ek (awkward.Array): k-point ECF

    Returns:
        awkward.Array
    """

    return ei*ek/ej**2

def calculate_ecf_d(ei, ej, ek):
    """Calculate D-series ECFs.
    
    ei (resp. ej, ek) is the i-point (resp. j-point, k-point) normalized ECF.
    To calculate D-series ECFs: j=i+1 and k=i+2.

    Args:
        ei (awkward.Array): i-point ECF
        ej (awkward.Array): j-point ECF
        ek (awkward.Array): k-point ECF

    Returns:
        awkward.Array
    """

    return ei**3*ek/ej**3

def calculate_ecf_m(v1ei, v1ej):
    """Calculate M-series ECFs.
    
    v1ei (resp. v1ej) is the first i-point (resp. j-point) normalized
    generalized ECFs.
    To calculate M-series ECFs: j=i+1.

    Args:
        v1ei (awkward.Array): 1st i-point ECFG
        v1ej (awkward.Array): 1st j-point ECFG

    Returns:
        awkward.Array
    """

    return v1ej/v1ei

def calculate_ecf_n(v1ei, v2ej):
    """Calculate N-series ECFs.
    
    v1ei (resp. v2ej) is the first i-point (resp.i second j-point) normalized
    generalized ECFs.
    To calculate N-series ECFs: j=i+1.

    Args:
        v1ei (awkward.Array): 1st i-point ECFG
        v2ej (awkward.Array): 2nd j-point ECFG

    Returns:
        awkward.Array
    """

    return v2ej/v1ei**2


def calculate_delta_phi(jets, met, jagged=True):
    """Calculate azimuthal angle between jets and MET.

    Args:
        jets (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with name PtEtaPhiMLorentzVector.
        met (Awkward.Array):
            Ak array with 1 axis with field phi.

    Returns:
        awkward.Array
    """

    met = akutl.broadcast(met, jets)[0]
    return vecutl.delta_phi(jets, met)
