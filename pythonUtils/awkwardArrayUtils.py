import awkward as ak


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
        >>> ak_array1 = [ [0, 3], [5], [2] ]
        >>> ak_array2 = [ [3, 3], [0], [1] ]
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


