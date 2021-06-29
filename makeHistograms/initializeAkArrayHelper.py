import awkward as ak
import sys

sys.path.append("../pythonUtils/")
import awkwardArrayUtilities as akutl
import coffeaUtilities as cfutl

 
def read_basic_variables(events, jet, jetCollection, jetVariables, jaggedVarArrays, varArrays):
    varArrays[jet+"Jet_n"] = cfutl.get_from_events(events, "n"+jetCollection)
    for var in jetVariables:
        jaggedVarArrays[jet+"Jet_"+var] = cfutl.get_from_events(events, jetCollection+"_"+var)
        if jaggedVarArrays[jet+"Jet_"+var] is not None:
            varArrays[jet+"Jet_"+var] = ak.flatten(jaggedVarArrays[jet+"Jet_"+var])
        else:
            varArrays[jet+"Jet_"+var] = None
    return


def make_njet_masks(events, jet, jetCollection, njet_max, masks, variable):
    njet_branch = cfutl.get_from_events(events, "n"+jetCollection)
    jet_variable_branch = cfutl.get_from_events(events, jetCollection+"_"+variable)
    for njet in range(1, njet_max+1):
        geNJet  =  "ge" + str(njet) + jet   # shorthand
        if njet_branch is None:
            masks[geNJet] = (njet_branch > njet-1)
        else:
            masks[geNJet] = (ak.count(jet_variable_branch, axis=1) > njet-1)
        masks[geNJet+"_bc"] = ak.flatten(ak.broadcast_arrays(masks[geNJet], jet_variable_branch)[0])
    return


def compute_variables_per_jet(jet_variables, jet, njet, jaggedVarArrays, varArrays, masks):
    geNJet  =  "ge" + str(njet) + jet   # shorthand

    for ijet in range(1, njet+1):
        sijet = str(ijet)
        for var in jet_variables:
            if jaggedVarArrays[jet+"Jet_"+var] is not None:
                varArrays[jet+"Jet"+sijet+"_"+var+"_"+geNJet] = jaggedVarArrays[jet+"Jet_"+var][masks[geNJet]][:, ijet-1]
            else:
                varArrays[jet+"Jet"+sijet+"_"+var+"_"+geNJet] = None
    return


def make_masked_jagged_array_shorthands(jet, njet, masks, jaggedVarArrays):
    geNJet = "ge" + str(njet) + jet   # shorthand
    maskGeNJet = masks[geNJet]        # shorthand

    # MET 4-vector for events with at least njet jets
    jaggedVarArrays["MET_4vector_"+geNJet] = jaggedVarArrays["MET_4vector"][maskGeNJet]

    # jet 4-vector for events with at least njet jets
    jaggedVarArrays[jet+"Jet_4vector_"+geNJet] = jaggedVarArrays[jet+"Jet_4vector"][maskGeNJet]

    return


def compute_nsubjettiness_ratio(events, tauRatioBranchName, tauAVariableName, tauBVariableName, tauRatioVariableName, jaggedVarArrays, varArrays):
    """Read tauRatio=tauA/tauB if available in tree or compute it on the fly if not.
    
    Args:
        events (ak.Array)
        tauRatioBranchName (str)
        tauAVariableName (str)
        tauBVariableName (str)
        tauRatioVariableName (str)
        jaggedVarArrays (ak.Array)
        varArrays (ak.Array)
    """

    tauRatioBranch = cfutl.get_from_events(events, tauRatioBranchName)
    if tauRatioBranch is not None:
        jaggedVarArrays[tauRatioBranchName] = tauRatioBranch
        varArrays[tauRatioBranchName] = ak.flatten(jaggedVarArrays[tauRatioBranchName])
    elif (varArrays[tauAVariableName] is not None) and (varArrays[tauBVariableName] is not None):
        jaggedVarArrays[tauRatioVariableName] = akutl.divide_ak_arrays(jaggedVarArrays[tauAVariableName], jaggedVarArrays[tauBVariableName])
        varArrays[tauRatioVariableName] = ak.flatten(jaggedVarArrays[tauRatioVariableName])
    else:
        jaggedVarArrays[tauRatioVariableName] = None
        varArrays[tauRatioVariableName] = None

    return



