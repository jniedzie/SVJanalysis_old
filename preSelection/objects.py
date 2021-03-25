import numpy as np
import awkward as ak

    
## Functions needed to take care of the reindixing of jets
## e.g. modification of *_jetIdx variables due to good jets selection

def true_indices(akArray):
    """
    Input: ak array with boolean values corresponding with structure of a branch of an Events tree
           or an iterable with similar jagged structure.
           e.g. [ [True, False, True], [True], [], [False, True] ]
    Return a jagged list with indices of the jagged array that are true.
           e.g. [ [0, 3], [0], [], [1] ]
    """

    trueIndices = [ [ idx for idx, x in enumerate(y) if x ] for y in akArray ]
    return trueIndices

        
def is_in(array1, array2):
    """
    Input: 2 jagged arrays
    Return a boolean jagged array
    e.g. array1 = [ [0, 1, 1, 2, 2, 3], [0, 0, 1], [0, 1, 1, 2] ]
         array2 = [ [0, 3], [], [0, 2] ]
         return [ [True, False, False, False, False, True], [False, False, False] [True, False, False, True]]
    """

    akArrayIsIn = ak.Array([ [ True if x in y2 else False for x in y1 ] for y1, y2 in zip(array1, array2) ])
    return akArrayIsIn


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


def make_new_jet_indices(filter_, jetIndices):
    """
    Input:
      * filter_   : jagged array with boolean values corresponding the selected jets
      * jetIndices: jagged array with jet indices
    Return jagged array of jet indices for selected jets
    """

    newJetIndices = []
    for filterRow, jetIndicesRow in zip(filter_, jetIndices):
        filterRow = ak.to_numpy(filterRow)
        X = [ idx for idx, x in enumerate(filterRow) ]
        Y = np.cumsum(filterRow)-1
        newJetIndicesRow = [ discrete_function(X, Y, x) for x in jetIndicesRow ]
        newJetIndices.append(newJetIndicesRow)
    return ak.Array(newJetIndices)



class GoodFatJets():
    """Make collection of good AK8 jets."""

    def __init__(self, events):

        # Definition of good AK8 jets
        self.filter_ = (events.FatJet_pt > 200) & (np.abs(events.FatJet_eta) < 2.4)

        # AK8 jets attributes
        self.pt        = events.FatJet_pt       [self.filter_]
        self.eta       = events.FatJet_eta      [self.filter_]
        self.phi       = events.FatJet_phi      [self.filter_]
        self.mass      = events.FatJet_mass     [self.filter_]
        self.msoftdrop = events.FatJet_msoftdrop[self.filter_]
        self.n2b1      = events.FatJet_n2b1     [self.filter_]
        self.n3b1      = events.FatJet_n3b1     [self.filter_]
        self.tau1      = events.FatJet_tau1     [self.filter_]
        self.tau2      = events.FatJet_tau2     [self.filter_]
        self.tau3      = events.FatJet_tau3     [self.filter_]
        self.tau4      = events.FatJet_tau4     [self.filter_]

        self.n = ak.count(self.pt, axis=1)
        

class GoodFatJetPFCands():
    """Make collection of PF candidates for good AK8 jets."""

    def __init__(self, events, goodFatJets, inputFileType):

        # PF candidates for good AK8 jets only
        if inputFileType == "PFnano102X":
            # Indices of the good jets
            goodJetIndices = true_indices(goodJets.filter_)
            # Boolean ak array to select only PF candidates in good jets
            self.filter_ = is_in(events.JetPFCands_jetIdx, goodJetIndices)
            # Get jet indices of the PF candidates in good jets
            PFCands_jetIdx = events.JetPFCands_jetIdx[self.filter_]
            # Some good jets have been removed, update the indices to good jets
            self.jetIdx = make_new_jet_indices(goodJets.filter_, PFCands_jetIdx)

            self.eta     = events.FatJetPFCands_eta    [self.filter_]
            self.mass    = events.FatJetPFCands_mass   [self.filter_]
            self.phi     = events.FatJetPFCands_phi    [self.filter_]
            self.pt      = events.FatJetPFCands_pt     [self.filter_]
            self.trkChi2 = events.FatJetPFCands_trkChi2[self.filter_]
            self.vtxChi2 = events.FatJetPFCands_vtxChi2[self.filter_]
            self.charge  = events.FatJetPFCands_charge [self.filter_]
            self.pdgId   = events.FatJetPFCands_pdgId  [self.filter_]

        elif inputFileType == "PFnano106X":
            # Indices of the good fat jets
            goodJetIndices = true_indices(goodFatJets.filter_)
            # Boolean ak array to select only PF candidates in good fat jets
            self.filter_ = is_in(events.JetPFCandsAK8_jetIdx, goodJetIndices)
            # Get jet indices of the PF candidates in good fat jets
            PFCands_jetIdx = events.JetPFCandsAK8_jetIdx[self.filter_]
            # Some fat jets may have been removed, update the indices to good jets
            self.jetIdx = make_new_jet_indices(goodFatJets.filter_, PFCands_jetIdx)
            # Get PF candidates indices in good fat jets
            self.candIdx = events.JetPFCandsAK8_candIdx[self.filter_]
            # Number of AK8 PF candidates
            self.nAK8 = ak.count(self.jetIdx, axis=1)

            # Cannot keep JetPFCands with index self.filter_ because some PF Cands not in self.filter_ may
            # also be clustered in AK4 jets!
            # This could be improved by keeping only PF Cands in both AK4 and AK8 jets and updating
            # JetPFCandsAK4_candIdx and JetPFCandsAK8_candIdx collections!
            self.eta     = events.JetPFCands_eta    
            self.mass    = events.JetPFCands_mass   
            self.phi     = events.JetPFCands_phi    
            self.pt      = events.JetPFCands_pt     
            self.trkChi2 = events.JetPFCands_trkChi2
            self.vtxChi2 = events.JetPFCands_vtxChi2
            self.charge  = events.JetPFCands_charge 
            self.pdgId   = events.JetPFCands_pdgId  

        self.n = ak.count(self.pt, axis=1)
 

class GoodJets():
    """Make collection of good AK4 jets."""

    def __init__(self, events):

        # Definition of good AK4 jets
        self.filter_ = (events.Jet_pt > 30) & (np.abs(events.Jet_eta) < 2.4)

        # AK4 jets attributes
        self.mass  = events.Jet_mass [self.filter_]
        self.pt    = events.Jet_pt   [self.filter_]
        self.eta   = events.Jet_eta  [self.filter_]
        self.phi   = events.Jet_phi  [self.filter_]
        self.chHEF = events.Jet_chHEF[self.filter_]
        self.neHEF = events.Jet_neHEF[self.filter_]

        self.n = ak.count(self.pt, axis=1)
        

class GoodJetPFCands():
    """Make collection of PF candidates for good AK4 jets."""

    def __init__(self, events, goodJets, inputFileType):

        # PF candidates for good AK4 jets only
        if inputFileType == "PFnano102X":
            # Indices of the good jets
            goodJetIndices = true_indices(goodJets.filter_)
            # Boolean ak array to select only PF candidates in good jets
            self.filter_ = is_in(events.JetPFCands_jetIdx, goodJetIndices)
            # Get jet indices of the PF candidates in good jets
            PFCands_jetIdx = events.JetPFCands_jetIdx[self.filter_]
            # Some jets may have been removed, update the indices to good jets
            self.jetIdx = make_new_jet_indices(goodJets.filter_, PFCands_jetIdx)

            # Get JetPFCands variables for PF candidates in good jets
            self.eta     = events.JetPFCands_eta    [self.filter_]
            self.mass    = events.JetPFCands_mass   [self.filter_]
            self.phi     = events.JetPFCands_phi    [self.filter_]
            self.pt      = events.JetPFCands_pt     [self.filter_]
            self.trkChi2 = events.JetPFCands_trkChi2[self.filter_]
            self.vtxChi2 = events.JetPFCands_vtxChi2[self.filter_]
            self.charge  = events.JetPFCands_charge [self.filter_]
            self.pdgId   = events.JetPFCands_pdgId  [self.filter_]

        elif inputFileType == "PFnano106X":
            # Indices of the good jets
            goodJetIndices = true_indices(goodJets.filter_)
            # Boolean ak array to select only PF candidates in good jets
            self.filter_ = is_in(events.JetPFCandsAK4_jetIdx, goodJetIndices)
            # Get jet indices of the PF candidates in good jets
            PFCands_jetIdx = events.JetPFCandsAK4_jetIdx[self.filter_]
            # Some good jets have been removed, update the indices to good jets
            self.jetIdx = make_new_jet_indices(goodJets.filter_, PFCands_jetIdx)
            # Get PF candidates indices in good jets
            self.candIdx = events.JetPFCandsAK4_candIdx[self.filter_]
            # Number of AK4 PF candidates
            self.nAK4 = ak.count(self.jetIdx, axis=1)

            # Cannot keep JetPFCands with index self.filter_ because some PF Cands not in self.filter_ may
            # also be clustered in AK8 jets!
            # This could be improved by keeping only PF Cands in both AK4 and AK8 jets and updating
            # JetPFCandsAK4_candIdx and JetPFCandsAK8_candIdx collections!
            self.eta     = events.JetPFCands_eta    
            self.mass    = events.JetPFCands_mass   
            self.phi     = events.JetPFCands_phi    
            self.pt      = events.JetPFCands_pt     
            self.trkChi2 = events.JetPFCands_trkChi2
            self.vtxChi2 = events.JetPFCands_vtxChi2
            self.charge  = events.JetPFCands_charge 
            self.pdgId   = events.JetPFCands_pdgId  

        self.n = ak.count(self.pt, axis=1)
 
