import numpy as np


def delta_r(obj1, obj2):
    return obj1.delta_r(obj2)

def delta_phi(obj1, obj2):
    return abs(obj1.delta_phi(obj2))

def abs_delta_phi(obj1, obj2):
    return abs(delta_phi(obj1, obj2))

def delta_eta(obj1, obj2):
    return abs(obj1.eta - obj2.eta)

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


