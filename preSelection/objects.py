import numpy as np
import awkward as ak
import sys

sys.path.append("../utilities/")
import nameUtilities as nameutl
    

## Functions needed to take care of the reindixing of jets
## e.g. modification of *_jetIdx variables due to good jets selection

def true_indices(ak_array):
    """Returns True indice in jagged array.

    Args:
        ak_array (awkard.highlevel.Array[bool]): 2D ak array with boolean values
           with structure of a branch of an Events tree or an iterable with
           similar jagged structure.

    Returns:
        list[list[int]]: a jagged list with indices of the jagged array that are true.

    Examples:
        >>> ak_array = [ [True, False, True], [True], [], [False, True] ]
        >>> true_indices(ak_array)
        [ [0, 2], [0], [], [1] ]
    """

    true_indices = [ [ idx for idx, x in enumerate(y) if x ] for y in ak_array ]
    return true_indices

        
def is_in(array1, array2):
    """
    Args:
        array1 (awkard.highlevel.Array[int])
        array2 (awkard.highlevel.Array[int])

    Returns:
        awkard.highlevel.Array[bool]

    Examples:
        >>> array1 = [ [0, 1, 1, 2, 2, 3], [0, 0, 1], [0, 1, 1, 2] ]
        >>> array2 = [ [0, 3], [], [0, 2] ]
        >>> is_in(array1, array2)
        [ [True, False, False, False, False, True], [False, False, False] [True, False, False, True]]
    """

    return ak.Array([ [ True if x in y2 else False for x in y1 ] for y1, y2 in zip(array1, array2) ])


def discrete_function(X, Y, x):
    """Discrete function f(xi) = yi.

    Args:
        X (list[T], tuple[T], array[T]): (x1, x2, ..., xn)
        Y (list[U], tuple[U], array[U]): (y1, y2, ..., yn)
        x (T): element in X
    
    Returns:
        U: f(x) = y
    """

    X = list(X)
    Y = list(Y)
    idx = X.index(x)

    return Y[idx]


def make_new_jet_indices(filter_, jet_indices):
    """Update jet indices after selecting some jets.

    Args
        filter_ (awkard.highlevel.Array):
            jagged array with boolean values corresponding the selected jets
        jetIndices (awkard.highlevel.Array): 
            jagged array with jet indices

    Returns:
        awkard.highlevel.Array: jagged array of jet indices for selected jets
    """

    new_jet_indices = []
    for filter_row, jet_indices_row in zip(filter_, jet_indices):
        filter_row = ak.to_numpy(filter_row)
        X = [ idx for idx, x in enumerate(filter_row) ]
        Y = np.cumsum(filter_row)-1
        new_jet_indices_row = [ discrete_function(X, Y, x) for x in jet_indices_row ]
        new_jet_indices.append(new_jet_indices_row)
    return ak.Array(new_jet_indices)


def make_ak_array_collection(collection, variables):
    """Make jagged array with fields.

    Args:
        collection (object having awkward.highlevel.Array attributes)
        variables (dict[str]):
            Keys are attributes names of the object collection
            Values are fields names to use in the ak array

    Returns:
        awkward.highlevel.Array
    """

    return ak.zip({
               final_variable_name: getattr(collection, initial_variable_name)
               for final_variable_name, initial_variable_name in variables.items()
               if hasattr(collection, initial_variable_name)
           })



class Jets():

    def __init__(self, events, jet_algo_name, input_file_type, cut=None):
        """
        Args:
            events (awkward.highlevel.Array): the Events TTree opened with uproot
            jet_algo_name (str)
            input_file_type (str)
            cut (str or None)

        Returns:
            None
        """

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
        """
        Args:
            cut (awkward.highlevel.Array[bool])

        Returns:
            None
        """

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



class PfCands():

    def __init__(self, events, input_file_type):
        """
        Args:
            events (awkward.highlevel.Array): the Events TTree opened with uproot
            input_file_type (str)

        Returns:
            None
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



class JetPfCandsMatchingTable():

    def __init__(self, events, good_jets, jet_algo_name, input_file_type):
        """
        Args:
            events (awkward.highlevel.Array): the Events TTree opened with uproot
            good_jets (Jets)
            jet_algo_name (str)
            input_file_type (str)

        Returns:
            None
        """

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
            self.n = ak.count(self.jetIdx, axis=1)

        self.variables = ["jetIdx", "pFCandsIdx"]

