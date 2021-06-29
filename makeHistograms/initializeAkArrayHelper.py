import awkward as ak
import sys

sys.path.append("../pythonUtils/")
import awkwardArrayUtilities as akutl
import coffeaUtilities as cfutl

 
def read_basic_variables(events, jet, jet_collection, jet_variables, jagged_var_arrays, var_arrays):
    """Read basic jet variables.

    Args:
        events (ak.Array)
        jet (str)
        jet_collection (str)
        jet_variables (list[str])
        jagged_var_arrays (dict[str, ak.Array])
        var_arrays (dict[str, ak.Array])

    Returns:
        None
    """

    var_arrays[jet+"Jet_n"] = cfutl.get_from_events(events, "n"+jet_collection)
    for var in jet_variables:
        jagged_var_arrays[jet+"Jet_"+var] = cfutl.get_from_events(events, jet_collection+"_"+var)
        if jagged_var_arrays[jet+"Jet_"+var] is not None:
            var_arrays[jet+"Jet_"+var] = ak.flatten(jagged_var_arrays[jet+"Jet_"+var])
        else:
            var_arrays[jet+"Jet_"+var] = None
    return


def make_njet_masks(events, jet, jet_collection, njet_max, masks, variable):
    """Make masks to filter events depending on their number of jets.

    Args:
        events (ak.Array)
        jet (str)
        jet_collection (str)
        njet_max (int)
        masks (dict[str, ak.Array])
        variables (str)

    Returns:
        None
    """

    njet_branch = cfutl.get_from_events(events, "n"+jet_collection)
    jet_variable_branch = cfutl.get_from_events(events, jet_collection+"_"+variable)
    for njet in range(1, njet_max+1):
        ge_njet  =  "ge" + str(njet) + jet   # shorthand
        if njet_branch is None:
            masks[ge_njet] = (njet_branch > njet-1)
        else:
            masks[ge_njet] = (ak.count(jet_variable_branch, axis=1) > njet-1)
        masks[ge_njet+"_bc"] = ak.flatten(ak.broadcast_arrays(masks[ge_njet], jet_variable_branch)[0])
    return


def compute_variables_per_jet(jet_variables, jet, njet, jagged_var_arrays, var_arrays, masks):
    """For all basic variables, compute them per jet index.

    For instance, break down ak4Jet_pt into ak4Jet1_pt, ak4Jet2_pt, etc...

    Args:
        jet_variables (list[str])
        jet (str)
        njet (int)
        jagged_var_arrays (dict[str, ak.Array])
        var_arrays (dict[str, ak.Array])
        masks (dict[str, ak.Array])

    Returns:
        None
    """

    ge_njet  =  "ge" + str(njet) + jet   # shorthand

    for ijet in range(1, njet+1):
        sijet = str(ijet)
        for var in jet_variables:
            if jagged_var_arrays[jet+"Jet_"+var] is not None:
                var_arrays[jet+"Jet"+sijet+"_"+var+"_"+ge_njet] = jagged_var_arrays[jet+"Jet_"+var][masks[ge_njet]][:, ijet-1]
            else:
                var_arrays[jet+"Jet"+sijet+"_"+var+"_"+ge_njet] = None
    return


def make_masked_jagged_array_shorthands(jet, njet, masks, jagged_var_arrays):
    """Define some shorthands.

    Args:
        jet (str)
        njet (int)
        masks (dict[str, ak.Array])
        jagged_var_arrays (dict[str, ak.Array])

    Returns:
        None
    """

    ge_njet = "ge" + str(njet) + jet   # shorthand
    mask_ge_njet = masks[ge_njet]        # shorthand

    # MET 4-vector for events with at least njet jets
    jagged_var_arrays["MET_4vector_"+ge_njet] = jagged_var_arrays["MET_4vector"][mask_ge_njet]

    # jet 4-vector for events with at least njet jets
    jagged_var_arrays[jet+"Jet_4vector_"+ge_njet] = jagged_var_arrays[jet+"Jet_4vector"][mask_ge_njet]

    return


def compute_nsubjettiness_ratio(events, tau_ratio_branch_name, tauA_variable_name, tauB_variable_name, tau_ratio_variable_name, jagged_var_arrays, var_arrays):
    """Read tauRatio=tauA/tauB if available in tree or compute it on the fly if not.
    
    Args:
        events (ak.Array)
        tau_ratio_branch_name (str)
        tauA_variable_name (str)
        tauB_variable_name (str)
        tau_ratio_variable_name (str)
        jagged_var_arrays (ak.Array)
        var_arrays (ak.Array)

    Returns:
        None
    """

    tau_ratio_branch = cfutl.get_from_events(events, tau_ratio_branch_name)
    if tau_ratio_branch is not None:
        jagged_var_arrays[tau_ratio_branch_name] = tau_ratio_branch
        var_arrays[tau_ratio_branch_name] = ak.flatten(jagged_var_arrays[tau_ratio_branch_name])
    elif (var_arrays[tauA_variable_name] is not None) and (var_arrays[tauB_variable_name] is not None):
        jagged_var_arrays[tau_ratio_variable_name] = akutl.divide_ak_arrays(jagged_var_arrays[tauA_variable_name], jagged_var_arrays[tauB_variable_name])
        var_arrays[tau_ratio_variable_name] = ak.flatten(jagged_var_arrays[tau_ratio_variable_name])
    else:
        jagged_var_arrays[tau_ratio_variable_name] = None
        var_arrays[tau_ratio_variable_name] = None

    return



def make_pairs(n):
    """Return pairs of distinct integers from {0, 1, ...n}.
    
    Args:
        n (int)

    Returns:
        list[[int, int]]
    """

    pairs = []
    for i1 in range(n-1):
        for i2 in range(i1+1, n):
            pairs.append([i1, i2])
    return pairs

