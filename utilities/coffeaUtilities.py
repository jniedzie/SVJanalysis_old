from coffea import processor
import awkward as ak
import numpy as np


def get_from_events(events, branch_name):
    """Return the branch from the events TTree if it exists, else return None.

    Args:
        events (ak.Array): the Events TTree opened with uproot
       branch_name (str): nane of the branch to read
    """

    if branch_name in events.fields:
        return events[branch_name]
    else:
        return None


def column_accumulator(type_):
    """Coffea processor column accumulator shorthand.

    Args:
        type_ (type): Type of the elements in the accumulator

    Returns:
        coffea.processor.column_accumulator: accumulator for type in argument
    """

    return processor.column_accumulator(np.array([], dtype=type_))


def value_accumulator(type_, initial=0):
    return processor.value_accumulator(type_, initial=initial)


def accumulate(ak_array):
    """Accumulate ak array into a coffea processor column accumulator.

    The array to accumulate can be either jagged or rectangular.
    If the array is rectangular (resp. jagged), the type can (resp. cannot)
    be retrieved afterwards.
    Jagged arrays are accumulated as object type.

    Args:
        ak_array(ak.Array)

    Returns:
        coffea.processor.column_accumulator
    """

    ak_type = str(ak.type(ak_array))
    if "var" in ak_type:  # if the array is jagged
        return processor.column_accumulator(np.array(ak.to_list(ak_array), dtype=object))
    else:
        return processor.column_accumulator(ak.to_numpy(ak_array))


