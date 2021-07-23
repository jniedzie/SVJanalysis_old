import numpy as np
import awkward as ak
import sys

sys.path.append("../utilities/")
import nameUtilities as nameutl
    
## TODO: doctrings

## Functions needed to take care of the reindixing of jets
## e.g. modification of *_jetIdx variables due to good jets selection

def true_indices(ak_array):
    """
    Input: ak array with boolean values corresponding with structure of a branch of an Events tree
           or an iterable with similar jagged structure.
           e.g. [ [True, False, True], [True], [], [False, True] ]
    Return a jagged list with indices of the jagged array that are true.
           e.g. [ [0, 3], [0], [], [1] ]
    """

    true_indices = [ [ idx for idx, x in enumerate(y) if x ] for y in ak_array ]
    return true_indices

        
def is_in(array1, array2):
    """
    Input: 2 jagged arrays
    Return a boolean jagged array
    e.g. array1 = [ [0, 1, 1, 2, 2, 3], [0, 0, 1], [0, 1, 1, 2] ]
         array2 = [ [0, 3], [], [0, 2] ]
         return [ [True, False, False, False, False, True], [False, False, False] [True, False, False, True]]
    """

    return ak.Array([ [ True if x in y2 else False for x in y1 ] for y1, y2 in zip(array1, array2) ])


def discrete_function(X, Y, x):
    """
    Input X = (x1, x2, ..., xn)
          Y = (y1, y2, ..., yn)
          x in X
    Discrete function f(xi) = yi 
    Return f(x)
    """

    X = list(X)
    Y = list(Y)
    idx = X.index(x)
    return Y[idx]


def make_new_jet_indices(filter_, jet_indices):
    """
    Input:
      * filter_   : jagged array with boolean values corresponding the selected jets
      * jetIndices: jagged array with jet indices
    Return jagged array of jet indices for selected jets
    """

    new_jet_indices = []
    for filter_row, jet_indices_row in zip(filter_, jet_indices):
        filter_row = ak.to_numpy(filter_row)
        X = [ idx for idx, x in enumerate(filter_row) ]
        Y = np.cumsum(filter_row)-1
        new_jet_indices_row = [ discrete_function(X, Y, x) for x in jet_indices_row ]
        new_jet_indices.append(new_jet_indices_row)
    return ak.Array(new_jet_indices)


def count(ak_array, axis=1):
    """Returns None if ak array is empty, else returns ak.count."""

    if ak.size(ak.flatten(ak_array)) == 0:
        return None
    else:
        return ak.count(ak_array, axis=axis)


def make_ak_array_collection(collection, variables):
    """ ... """

    return ak.zip({
               final_variable_name: getattr(collection, initial_variable_name)
               for final_variable_name, initial_variable_name in variables.items()
               if hasattr(collection, initial_variable_name)
           })



class GoodJets():
    """Make collection of good AK8 jets."""

    def __init__(self, events, jet_algo_name, input_file_type, cut=None):

        if input_file_type == "PFNanoAOD_106X_v01" or input_file_type == "PFNanoAOD_106X_v02":
            jet_collection_name = nameutl.jet_algo_name_to_jet_collection_name(jet_algo_name)
            jets = getattr(events, jet_collection_name)

            # Jets attributes
            jet_variables = {
                "pt": "pt",
                "eta": "eta",
                "phi": "phi",
                "mass": "mass",
                "msoftdrop": "msoftdrop",
                "n2b1": "n2b1",
                "n3b1": "n3b1",
                "tau1": "tau1",
                "tau2": "tau2",
                "tau3": "tau3",
                "tau4": "tau4",
                "chHEF": "chHEF",
                "neHEF": "neHEF",
            }

        else:
            # Define here jets and jet_variables
            pass

        self.jet = make_ak_array_collection(jets, jet_variables)
        self.variables = self.jet.fields


        ## Definition of good jets
        # Can define complicated cuts here
        if cut is None:
            if jet_algo_name.lower() == "ak4":
                pass
            elif jet_algo_name.lower() == "ak8":
                pass
        # Or use the cut argument
        else:
            self.filter_ = eval(cut.replace("jet", "self.jet"))
            self.jet = self.jet[self.filter_]

        ## Number of jets after cuts
        self.n = ak.count(getattr(self, self.variables[0]), axis=1)



    def apply_cut(self, cut):
        self.jet = self.jet[cut]
        self.n = self.n[cut]
        self.filter_ = self.filter_[cut]


    @property
    def pt(self):
        return self.jet.pt

    @property
    def eta(self):
        return self.jet.eta

    @property
    def phi(self):
        return self.jet.phi

    @property
    def mass(self):
        return self.jet.mass

    @property
    def msoftdrop(self):
        return self.jet.msoftdrop

    @property
    def n2b1(self):
        return self.jet.n2b1

    @property
    def n3b1(self):
        return self.jet.n3b1

    @property
    def tau1(self):
        return self.jet.tau1

    @property
    def tau2(self):
        return self.jet.tau2

    @property
    def tau3(self):
        return self.jet.tau3

    @property
    def tau4(self):
        return self.jet.tau4

    @property
    def chHEF(self):
        return self.jet.chHEF

    @property
    def neHEF(self):
        return self.jet.neHEF



class GoodPfCands():

    def __init__(self, events, input_file_type):
        """Get .

        Args:
            events (ak.Array): the Events TTree opened with uproot.
            input_file_type (str)

        Returns:
            awkward.Array: ak array with fields 
        """

        ## For PF nano AOD 106X and master branches
        if input_file_type == "PFNanoAOD_106X_v01" or input_file_type == "PFNanoAOD_106X_v02":
            if input_file_type == "PFNanoAOD_106X_v01":
                pf_cands = events.JetPFCands
            else:
                pf_cands = events.PFCands

            # PF cands attributes
            pf_cands_variables = {
                "pt": "pt",
                "eta": "eta",
                "phi": "phi",
                "mass": "mass",
                "charge": "charge",
                "pdgId": "pdgId",
                "trkChi2": "trkChi2",
                "vtxChi2": "vtxChi2",
            }

        else:
            # Define here pf_cands and pf_cands_variables
            pass

        self.pf_cands = make_ak_array_collection(pf_cands, pf_cands_variables)
        self.variables = self.pf_cands.fields

        ## Number of jets after cuts
        self.n = ak.count(getattr(self, self.variables[0]), axis=1)


    @property
    def pt(self):
        return self.pf_cands.pt

    @property
    def eta(self):
        return self.pf_cands.eta

    @property
    def phi(self):
        return self.pf_cands.phi

    @property
    def mass(self):
        return self.pf_cands.mass

    @property
    def charge(self):
        return self.pf_cands.charge

    @property
    def pdgId(self):
        return self.pf_cands.pdgId

    @property
    def trkChi2(self):
        return self.pf_cands.trkChi2

    @property
    def vtxChi2(self):
        return self.pf_cands.vtxChi2



class GoodJetPfCandsMatchingTable():
    """..."""

    def __init__(self, events, good_jets, jet_algo_name, input_file_type):

        if input_file_type == "PFNanoAOD_106X_v01" or input_file_type == "PFNanoAOD_106X_v02":

            if input_file_type == "PFNanoAOD_106X_v01":
                jet_algo_name = jet_algo_name.upper()
                table = getattr(events, "JetPFCands"+jet_algo_name)
                cand_idx_name = "candIdx"
            else:
                jet_collection_name = nameutl.jet_algo_name_to_jet_collection_name(jet_algo_name)
                table = getattr(events, jet_collection_name + "PFCands")
                cand_idx_name = "pFCandsIdx"

            # Indices of the good jets
            good_jet_indices = true_indices(good_jets.filter_)
            # Boolean ak array to select only PF candidates in good fat jets
            self.filter_ = is_in(table.jetIdx, good_jet_indices)
            # Get jet indices of the PF candidates in good fat jets
            pf_cands_jet_idx = table.jetIdx[self.filter_]
            # Some jets may have been removed, update the indices to good jets
            self.jetIdx = make_new_jet_indices(good_jets.filter_, pf_cands_jet_idx)
            # Get PF candidates indices in good fat jets
            self.pFCandsIdx = getattr(table, cand_idx_name)[self.filter_]
            # Number of PF candidates
            self.n = count(self.jetIdx, axis=1)

        self.variables = ["jetIdx", "pFCandsIdx"]

