import numpy as np
import awkward as ak


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

    def __init__(self, events, goodFatJets):

        # PF candidates for good AK8 jets only
        # Works only for 102X branch of PFnanoAOD
        self.filter_ = (events.FatJetPFCands_jetIdx < goodFatJets.n)

        self.jetIdx  = events.FatJetPFCands_jetIdx [self.filter_]
        self.eta     = events.FatJetPFCands_eta    [self.filter_]
        self.mass    = events.FatJetPFCands_mass   [self.filter_]
        self.phi     = events.FatJetPFCands_phi    [self.filter_]
        self.pt      = events.FatJetPFCands_pt     [self.filter_]
        self.trkChi2 = events.FatJetPFCands_trkChi2[self.filter_]
        self.vtxChi2 = events.FatJetPFCands_vtxChi2[self.filter_]
        self.charge  = events.FatJetPFCands_charge [self.filter_]
        self.pdgId   = events.FatJetPFCands_pdgId  [self.filter_]

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

    def __init__(self, events, goodJets):

        # PF candidates for good AK4 jets only
        # Works only for 102X branch of PFnanoAOD
        self.filter_ = (events.JetPFCands_jetIdx < goodJets.n)

        self.jetIdx  = events.JetPFCands_jetIdx [self.filter_]
        self.eta     = events.JetPFCands_eta    [self.filter_]
        self.mass    = events.JetPFCands_mass   [self.filter_]
        self.phi     = events.JetPFCands_phi    [self.filter_]
        self.pt      = events.JetPFCands_pt     [self.filter_]
        self.trkChi2 = events.JetPFCands_trkChi2[self.filter_]
        self.vtxChi2 = events.JetPFCands_vtxChi2[self.filter_]
        self.charge  = events.JetPFCands_charge [self.filter_]
        self.pdgId   = events.JetPFCands_pdgId  [self.filter_]

        self.n = ak.count(self.pt, axis=1)
 
