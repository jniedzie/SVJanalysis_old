import numpy as np
import awkward as ak
from coffea.nanoevents.methods import vector
# Needed so that ak.zip({"pt": [...], "eta": [...], "phi": [...], "mass": [...]},
#                         with_name="PtEtaPhiMLorentzVector")
# is understood as a PtEtaPhiMLorentzVector from coffea.nanoevents.methods.vector
ak.behavior.update(vector.behavior)


def make_PtEtaPhiMLorentzVector(pt, eta, phi, mass):
    """Take pt, eta, phi, mass awkward arrays and return the corresponding PtEtaPhiMLorentzVector."""

    vec = ak.zip(
        {
            "pt": pt,
            "eta": eta,
            "phi": phi,
            "mass": mass,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    return vec


def delta_r(obj1, obj2):
    return obj1.delta_r(obj2)

def delta_phi(obj1, obj2):
    return obj1.delta_phi(obj2)

def abs_delta_phi(obj1, obj2):
    return abs(delta_phi(obj1, obj2))

def delta_eta(obj1, obj2):
    return obj1.eta - obj2.eta

def abs_delta_eta(obj1, obj2):
    return abs(delta_eta(obj1, obj2))

def mass(obj1, obj2):
    return (obj1+obj2).mass

def pt(obj1, obj2):
    return (obj1+obj2).pt

def mt(jj, met):

    # Note that: jj.dot(met) = jj.pt * met.pt * np.cos(jj.delta_phi(met))
    return np.sqrt( jj.mass**2 + 2*( np.sqrt(jj.mass**2 + jj.pt**2) * met.pt - jj.dot(met)) )

def mt_wrong(jj, met):

    # Note that: jj.dot(met) = jj.pt * met.pt * np.cos(jj.delta_phi(met))
    return np.sqrt( jj.mass**2 + 2*( np.sqrt(jj.mass**2 + jj.pt**4) * met.pt - jj.pt**2 * met.pt * np.cos(jj.delta_phi(met)) ) )


