def jet_algo_2_jet_collection(jet_algo):
    """Helper function converting jet algorithm name into jet collection name.

    Args:
        jet_algo (str)

    Returns:
        str
    """

    table = {
        "ak4": "Jet",
        "ak8": "FatJet",
    }

    return table[jet_algo.lower()]


def jet_collection_2_jet_algo(jet_collection):
    """Helper function converting jet collection name into jet algorithm name.

    Args:
        jet_collection (str)

    Returns:
        str
    """

    table = {
        "jet": "ak4",
        "fatjet": "ak8",
    }

    return table[jet_collection.lower()]


