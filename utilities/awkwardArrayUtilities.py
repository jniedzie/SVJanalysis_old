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



def divide_ak_arrays(ak_array1, ak_array2, division_by_zero_value=1., verbose=False):
    """Makes the division of an ak array by another one.
    
    The arrays ak_array1 and ak_array2 must have the same jagged structure.
    If division by zero for some indices, a default value to use can be
    defined, see examples.
    The output array has the same jagged structure as the input ak arrays.
    

    Args:
        ak_array1 (awkward.Array[float])
        ak_array2 (awkward.Array[float]) 
        division_by_zero_value (float, optional, default=1.)
        verbose (bool, optional, default=False)

    Returns:
        awkward.Array[float]: ak_array1 / ak_array2

    Examples:
        >>> ak_array1 = ak.Array([ [0, 3], [5], [1] ])
        >>> ak_array2 = ak.Array([ [3, 3], [0], [2] ])
        >>> divide_ak_arrays(ak_array1, ak_array2)
        [ [0, 1], [1], [0.5] ]
    """

    is_not_zero = (ak_array2!=0.)
    if (not ak.all(is_not_zero)) and verbose:
        print("The following warning about true_divide can be safely ignored.")

    raw_division = ak_array1/ak_array2
    division = ak.where(is_not_zero, raw_division, division_by_zero_value*ak.ones_like(ak_array1))

    # This implementation seems slower:
    #division = ak.Array([ [ x1/x2 if x2 != 0. else division_by_zero_value for x1, x2 in zip(y1, y2) ] for y1, y2 in zip(ak_array1, ak_a0rray2) ])

    return division



def is_in(ak_array, list_):
    """Check whether the elements of an ak array are in a list of elements.
    
    Args:
        ak_array1 (awkward.Array[T])
        list_ (list[T])

    Returns:
        awkward.Array[bool]

    Examples:
        >>> ak_array = ak.Array([ [11, 22, -11], [22], [111, 211, -11, 11] ])
        >>> list_ = [11, -11]
        >>> is_in(ak_array, list_)
        [[True, False, True], [False], [False, False, True, True]]
    """

    ak_bool = False * ak.ones_like(ak_array, dtype=bool)
    for el in list_:
        ak_bool = ak_bool | (ak_array == el)

    return ak_bool

