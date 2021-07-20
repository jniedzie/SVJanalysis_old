def jet_algo_name_to_jet_collection_name(jet_algo):
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


def jet_collection_name_to_jet_algo_name(jet_collection):
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


