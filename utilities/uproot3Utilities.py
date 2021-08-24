import re
import sys

import numpy as np
import awkward as ak
from coffea import processor
import uproot3

import utilities
import awkwardArrayUtilities as akutl


def __send_casting_warning(dtype, new_dtype, branch_name):
    print("WARNING: Casting branch %s from %s to %s because uproot does not handle %s!" %(branch_name, dtype, new_dtype, dtype))
    return


def __send_max_value_warning(max_value_allowed, branch_name):
    print("WARNING: Values in branch %s should NOT exceed %d." %(branch_name, max_value_allowed))
    return


def __get_dtype(branch, branch_name="\b"):
    """Return branch data type.

    Args:
        branch (coffea.processor.accumulator.column_accumulator, branch, np.ndarray or
                awkward.highlevel.Array)
        branch_name (str, optional): banch name for printing out precise warnings

    Returns:
        str
    """

    if isinstance(branch, processor.accumulator.column_accumulator) \
    or isinstance(branch, np.ndarray):
        dtype = branch.dtype

    elif isinstance(branch, ak.highlevel.Array):
        dtype = ("%s" %ak.type(branch)).split("*")[-1][1:]

    else:
        supported_classes = ["coffea.processor.accumulator.column_accumulator", "numpy.ndarray", "awkward.highlevel.Array"]
        print("ERROR: Branch %s is from an unsupported class: %s" %(branch_name, type(branch)))
        print("Supported classes:")
        for class_ in supported_classes: print("\t%s" %class_)
        sys.exit()

    # uint type is not supported by uproot3
    # As long as values in the array do not exceed 2147483647 then casting brutaly from uint32 to int32 should be harmless
    re_search = re.search(r'uint([0-9]+)', dtype)
    if re_search:
        nbits = int(re_search.group(1))
        # int8 cannot be read properly by ROOT yet
        if nbits == 8:
            new_dtype = "int16"
            print_max_value_warning = False
        # int32, int64, ... are fine
        else:
            new_dtype = dtype[1:]
            print_max_value_warning = True

        __send_casting_warning(dtype, new_dtype, branch_name)
        if print_max_value_warning:
            max_value_allowed = 2**(nbits-1) - 1
            __send_max_value_warning(max_value_allowed, branch_name)
        dtype = new_dtype

    return dtype


def __get_size_branch_names(tree):
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


def __make_branch(branch, branch_name=""):
    """Make tree branch in appropriate format for writing with uproot3.

    Args:
        branch (coffea.processor.accumulator.column_accumulator,
                awkward.highlevel.Array, numpy.ndarray, numpy.float, numpy.int)
        branch_name (str)

    Returns:
        tuple (awkward0.array.jagged.JaggedArray or numpy.ndarray, int, str, str):
            branch content, number of events, data type, awkward type
    """

    if isinstance(branch, processor.accumulator.column_accumulator):
        branch = branch.value
    branch = akutl.obj_to_ak_array(branch)
    type_ = str(ak.type(branch))
    dtype = __get_dtype(branch, branch_name)
    size = ak.num(branch, axis=0)
    branch = ak.to_awkward0(branch)

    return branch, size, dtype, type_


def __make_branch_init(branch_name, dtype, size_branch_names):
    """Make tree branch initializer for uproot3.
    
    Args:
        branch_name (str)
        dtype (str)
        size_branch_names (list[str])

    Returns:
        uproot3.write.objects.TTree.newbranch or None
    """

    if not branch_name.startswith("n"):
        nbranch_name = "n"+str(branch_name.split("_")[0])
        if nbranch_name in size_branch_names:
            if dtype == "bool":
                new_dtype = "int16"
                __send_casting_warning(dtype, new_dtype, branch_name)
                dtype = new_dtype
            branch_init = uproot3.newbranch(dtype, size=nbranch_name)
        else:
            branch_init = uproot3.newbranch(dtype)
    else:
        if branch_name not in size_branch_names:
            branch_init = uproot3.newbranch(dtype)
        else:
            branch_init = None

    return branch_init


def __run_branch_checks(branches, branches_init, branches_size):
    """Run some checks on the branches and remove branches failing checks.

    Args:
        branches (dict[str, awkward0.array.jagged.JaggedArray or numpy.ndarray])
        branches_init (dict[str, uproot3.write.objects.TTree.newbranch])
        branches_size (dict[str, int])

    Returns:
        None
    """

    def delete_branch(branch_name, branches, branches_init):
        branches.pop(branch_name)
        if branch_name in branches_init.keys():
            branches_init.pop(branch_name)

        return

    size = utilities.most_common(branches_size.values())

    branches_to_delete = [ branch_name for branch_name, branch_size in branches_size.items() if branch_size != size ]
    for branch_name in branches_to_delete:
        print("\nWARNING: Branch %s has size %d while most branches have size %d." %(branch_name, branches_size[branch_name], size))
        print("         Deleting this branch!")
        delete_branch(branch_name, branches, branches_init)

    return


def make_branches(tree, debug=False):
    """Make branches and branches initailiazer to write a tree to a ROOT file.

    Args:
        tree (dict[str, awkward.Array] or coffea.processor.dict_accumulator)
        debug (bool)

    Returns:
        tuple (dict[str, awkward0.array.jagged.JaggedArray], dict[str, uproot3.write.objects.TTree.newbranch])
    """

     

    # Book dict to return
    branches = {}
    branches_init = {}

    # Book dictionary storing branch length for branch checks
    branches_size = {}

    # Find names of the branches storing the size of jagged array branches
    size_branch_names = __get_size_branch_names(tree)
    if debug:
        print("\nSize branch names:")
        print(size_branch_names)

    # Make branches and branches_init
    # Need to use ak0 and uproot3 because writing tree to ROOT file is not yet implemented in uproot4 (Feb. 2021)
    if debug:
        print("\nFilling branches and branches_init:")
    for branch_name, branch in tree.items():
        branches[branch_name], branches_size[branch_name], dtype, type_ = __make_branch(branch, branch_name)
        branch_init = __make_branch_init(branch_name, dtype, size_branch_names)
        if branch_init is not None: branches_init[branch_name] = branch_init
        if debug:
            print("%s:\ttype=%s   dtype=%s" %(branch_name, ak.type(branch), dtype))

    # Run checks: e.g. delete branches with bugs, like incorrrect length
    __run_branch_checks(branches, branches_init, branches_size)

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

