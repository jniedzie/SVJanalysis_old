import numpy as np
import awkward as ak
import uproot3


def obj_to_ak_array(obj):
    """Convert any input type to an awkward array.

    Args:
        obj (any type convertible to ak array)

    Returns:
        awkward.highlevel.Array
    """

    if isinstance(obj, ak.highlevel.Array):
        ak_array = obj
    elif isinstance(obj, np.ndarray):
        ak_array = ak.Array(obj)
    elif isinstance(obj, np.float32) or isinstance(obj, np.float64) or isinstance(obj, np.int32) or isinstance(obj, np.int64) \
         or isinstance(obj, float) or isinstance(obj, int):
        ak_array = ak.Array([obj])
    else:
        print("Unknown type %s" %(type(obj)))

    return ak_array


def get_size_branch_names(tree):
    """Find the names of branches storing the size of jagged array branches.

    By convention they are called n<CollectionName>.
    E.g. nJet is the name of the size branch for Jet_pt.

    Args:
        tree (dict or coffea.processor.dict_accumulator)
            Keys are the names of the branches in the tree

    Returns:
        list[str]
    """

    size_branch_names = []
    for branch_name, branch in tree.items():
        nbranch_name = "n"+str(branch_name.split("_")[0])
        if (not branch_name.startswith("n")) and (nbranch_name in tree.keys()) and (nbranch_name not in size_branch_names):
            size_branch_names.append(nbranch_name)

    return size_branch_names


def make_branches(tree, debug=False):
    """Make branches and branches initailiazer to write a tree to a ROOT file.

    Args:
        tree (dict[awkward.Array] or coffea.processor.dict_accumulator)
        debug (bool)

    Returns:
        tuple: (dict[awkward0.array.jagged.JaggedArray], dict[uproot3.write.objects.TTree.newbranch])
    """

    ## Running sanity checks
    tree_type = str(type(tree)).split("'")[1]
    allowed_tree_types = ["dict", "coffea.processor.accumulator.dict_accumulator"]
    if tree_type not in allowed_tree_types:
        print("ERROR: Unknown tree type: %s" %tree_type)
        print("       Allowed tree types:")
        for type_ in allowed_tree_types: print("\t%s" %type_)
        return None, None

    ## Book dict to return
    branches = {}
    branches_init = {}

    ## Find names of the branches storing the size of jagged array branches
    size_branch_names = get_size_branch_names(tree)
    if debug:
        print("Size branch names: ", size_branch_names)

    ## Making branches and branches_init
    #  Need to use ak0 and uproot3 because writing tree to ROOT file is not yet implemented in uproot4 (Feb. 2021)
    for branch_name, branch in tree.items():

        if tree_type == "coffea.processor.accumulator.dict_accumulator":
            branch = branch.value
            dtype = branch.dtype
            branch = obj_to_ak_array(branch)
        elif tree_type == "dict":
            branch = obj_to_ak_array(branch)
            dtype = ("%s" %ak.type(branch)).split("*")[-1][1:]

        if debug:
            print("%s: %s" %(branch_name, ak.type(branch)))


        # Defining branch
        branches[branch_name] = ak.to_awkward0(branch)

        # Defining branch_init
        if not branch_name.startswith("n"):
            nbranch_name = "n"+str(branch_name.split("_")[0])
            if nbranch_name in size_branch_names:
                # Case distinction for type of the jagged-array collections, treated as "object" in the processor
                # _pFCandsIdx and _jetIdx must be saved as integers ("i4") to use array at once syntax
                if branch_name.endswith("Idx"):
                    branches_init[branch_name] = uproot3.newbranch(np.dtype("i4"), size=nbranch_name)
                else:
                    branches_init[branch_name] = uproot3.newbranch(np.dtype("f8"), size=nbranch_name)
            else:
                branches_init[branch_name] = uproot3.newbranch(dtype)
        else:
            if branch_name not in size_branch_names:
                branches_init[branch_name] = uproot3.newbranch(dtype)

    return branches, branches_init


def write_tree_to_root_file(file_, tree_name, branches, branches_init):
    """Write histograms, stored in a coffea accumulator, to a ROOT file.

    Args:
        accumulator (coffea.processor.dict_accumulator)
        output_file_name (str): Name of the ROOT file to create
        debug (bool)

    Returns:
        None
    """

    file_[tree_name] = uproot3.newtree(branches_init)
    file_[tree_name].extend(branches)

    return


