import sys
import time
import argparse

import awkward as ak
import uproot3
import uproot

sys.path.append("../utilities/")
import uproot3Utilities as uproot3utl


def make_files_list(files_arg):
    """Make list of ROOT files to merge.

    Args:
        files_arg (str): files argument from the argument parser
            Comma separated ROOT file names e.g. file1.root,file2.root or
            text file name with a ROOT file name on each line.

    Returns:
        list[str]
    """

    # If list of ROOT files in txt file
    if not files_arg.endswith(".root"):
        with open (files_arg, "r") as txt_file:
            root_files = [ x.replace("\n", "") for x in txt_file.readlines() ]
    # Else we assume coma separated list of ROOT files
    else:
        root_files = files_arg.split(",")

    return root_files


def get_all_branch_names(file_name, tree_name):
    """Return all the branches name of the tree of a ROOT file.

    Args:
        file_name (str)
        tree_name (str)

    Returns:
        list[str]
    """

    return uproot.open(file_name)[tree_name].keys()


def get_unique_branch_names(file_names, tree_name):
    """Return branches names to read in a tree for each file for branches merging algo.

    For 2 input files file1.root and file2.root with branches nJet, Jet_pt and
    nJet, Jet_eta in the tree Events, the branches to read are nJet, Jet_pt for
    file1.root and Jet_eta for file2.root.

    Args:
        file_names (list[str])
        tree_name (str)

    Returns:
        dict[str, str]: Keys are file names, values are branch names
    """

    branches_to_read = {}
    branch_names_added = []
    for ifile, file_name in enumerate(file_names):
        branches_this_file = get_all_branch_names(file_name, tree_name)
        branches_to_read[file_name] = [ branch_name for branch_name in branches_this_file if branch_name not in branch_names_added ]
        branch_names_added = branch_names_added + branches_to_read[file_name]

    return branches_to_read


def branches_merging(file_names, tree_name):
    """Return tree branches for branches merging algorithm.

    Args:
        file_names (list[str])
        tree_name (str)

    Returns:
        dict[str, awkward.highlevel.Array]:
           Keys are branch names, values are branch arrays
    """

    branches_to_read = get_unique_branch_names(file_names, tree_name)
    branches = {}
    for file_name in file_names:
        tree = uproot.open(file_name + ":" + tree_name)
        for branch_name in branches_to_read[file_name]:
            branches[branch_name] = tree[branch_name].array()

    tree_branches = { k: v for k, v in branches.items() }

    return tree_branches


def events_merging(file_names, tree_name):
    """Return tree branches for events merging algorithm.

    Args:
        file_names (list[str])
        tree_name (str)

    Returns:
        dict[str, awkward.highlevel.Array]:
           Keys are branch names, values are branch arrays
    """

    branches_to_read = get_all_branch_names(file_names[0], tree_name)
    branches = { branch_name: [] for branch_name in branches_to_read }
    for file_name in file_names:
        tree = uproot.open(file_name + ":" + tree_name)
        for branch_name in branches_to_read:
            branches[branch_name].append(tree[branch_name].array())

    tree_branches = { k: ak.concatenate(v, axis=0) for k, v in branches.items() }

    return tree_branches


def sum_merging(file_names, tree_name):
    """Return tree branches for sum merging algorithm.

    Args:
        file_names (list[str])
        tree_name (str)

    Returns:
        dict[str, awkward.highlevel.Array]:
           Keys are branch names, values are branch arrays
    """

    branches_to_read = get_all_branch_names(file_names[0], tree_name)
    branches = { branch_name: 0. for branch_name in branches_to_read }
    for file_name in file_names:
        tree = uproot.open(file_name + ":" + tree_name)
        for branch_name in branches_to_read:
            branches[branch_name] += tree[branch_name].array()[0]

    tree_branches = { k: ak.Array([v]) for k, v in branches.items() }

    return tree_branches


def metadata_merging(file_names, tree_name, efficiency_branch_name="TotalCutEfficiency",
    original_format_branch_name="OriginalFormat", cross_section_branch_name="GenCrossSection",
    gen_weight_branch_name="Cutflow/AllCuts"):
    """Return tree branches for metadata merging algorithm.

    Args:
        file_names (list[str])
        tree_name (str)

    Returns:
        dict[str, awkward.highlevel.Array]:
           Keys are branch names, values are branch arrays
    """

    tree_branches = {}

    # Get cross section and original format from first file
    file_name = file_names[0]
    with uproot.open(file_name + ":" + tree_name) as tree:
        tree_branches[cross_section_branch_name] = tree[cross_section_branch_name].array()
        tree_branches[original_format_branch_name] = tree[original_format_branch_name].array()

    # Compute total cut efficiency for the ensemble of files
    sum_gen_weights = []
    efficiencies = []
    for file_name in file_names:
        with uproot.open(file_name) as file_:
            efficiencies.append(file_[tree_name][efficiency_branch_name].array()[0])
            if gen_weight_branch_name != "":
                sum_gen_weights.append(ak.sum(file_[gen_weight_branch_name].array(), axis=0))
            else:
                sum_gen_weights.append(1)

    sum_gen_weights_no_cut = [sum_gen_weight/efficiency for sum_gen_weight, efficiency in zip(sum_gen_weights, efficiencies)]
    tree_branches[efficiency_branch_name] = ak.Array([sum(sum_gen_weights) / sum(sum_gen_weights_no_cut)])

    return tree_branches


def efficiencies_merging(file_names, tree_name, weight_branch_name=""):
    """Return tree branches for efficiencies merging algorithm.

    Efficiencies from the different files are weighted by the sum of the
    weights read in the weight branch of the file.

    Args:
        file_names (list[str])
        tree_name (str)
        weight_branch_name (str, optional, default=""): "" for no weighting

    Returns:
        dict[str, awkward.highlevel.Array]:
           Keys are branch names, values are branch arrays
    """

    branches_to_read = get_all_branch_names(file_names[0], tree_name)
    branches = { branch_name: [] for branch_name in branches_to_read }
    weight_branches = []

    for file_name in file_names:
        file_ = uproot.open(file_name)
        tree = file_[tree_name]
        
        if weight_branch_name != "":
            weight_branches.append(ak.sum(file_[weight_branch_name].array(), axis=0))
        else:
            weight_branches.append(1)

        for branch_name in branches_to_read:
            branches[branch_name].append(tree[branch_name].array())

    tree_branches = { k: sum([w*b for w, b in zip(weight_branches, v)])/sum(weight_branches) for k, v in branches.items() }

    return tree_branches


def first_found_merging(file_names, tree_name):
    """Return tree branches for branches merging algorithm.

    Args:
        file_names (list[str])
        tree_name (str)

    Returns:
        dict[str, awkward.highlevel.Array]:
           Keys are branch names, values are branch arrays
    """

    tree_branches = {}
    for ifile, file_name in enumerate(file_names):
        file_ = uproot.open(file_name)
        if tree_name+";1" in file_.keys():
            tree = file_[tree_name].arrays()
            tree_branches = { k: tree[k] for k in tree.fields }
            break

    return tree_branches


def write_root_file(trees, merging_algos, output_file_name, debug):
    """Write histograms, stored in a coffea accumulator, to a ROOT file.

    Args:
        trees (dict[key, dict[key, awkward.highlevel.Array]])
            1st key is tree name
            2nd key is branch name
            Value is the branch (ak array)
        merging_algos (list[str])
        output_file_name (str)
        debug (bool)

    Returns:
        None
    """

    # Need to use compression=None because there is a bug with EFPs when the file is compressed:
    # More info at https://github.com/scikit-hep/uproot3/issues/506
    with uproot3.recreate(output_file_name, compression=None) as root_file:
        print("\nWriting down output ROOT file %s" %output_file_name)
        for (tree_name, tree), merging_algo in zip(trees.items(), merging_algos):

            # Making branches to write to tree
            branches, branches_init = uproot3utl.make_branches(tree, debug)

            if debug:
                print("\n%s branches keys:" %tree_name)
                print(branches.keys())
                print("\n%s branches_init keys:" %tree_name)
                print(branches_init.keys())

            # Save branches to ROOT file
            uproot3utl.write_tree_to_root_file(root_file, tree_name, branches, branches_init)
            print("TTree %s saved to output file" %tree_name)

    return


def main(file_names, tree_names, merging_algos, output_file, debug):

    print("Input files:")
    print(file_names)

    trees = {}
    for tree_name, merging_algo in zip(tree_names, merging_algos):

        if merging_algo == "branches":
            trees[tree_name] = branches_merging(file_names, tree_name)

        elif merging_algo == "events":
            trees[tree_name] = events_merging(file_names, tree_name)

        elif merging_algo == "sum":
            trees[tree_name] = sum_merging(file_names, tree_name)

        elif merging_algo == "metadata":
            trees[tree_name] = metadata_merging(file_names, tree_name)

        elif merging_algo.startswith("efficiencies"):
            # merging_algo can be efficiencies or efficiencies:<weight_branch_name>
            weight_branch_name = merging_algo.replace("efficiencies:", "")
            weight_branch_name = weight_branch_name.replace("efficiencies", "")
            trees[tree_name] = efficiencies_merging(file_names, tree_name, weight_branch_name)
           
        elif merging_algo == "firstFound":
            trees[tree_name] = first_found_merging(file_names, tree_name)

        else:
            print("ERROR: Unknown merging algorithm: %s" %merging_algo)
            sys.exit()

    write_root_file(trees, merging_algos, output_file, debug)



if __name__ == "__main__":

    tstart = time.time()

    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputfiles",
        help="Input ROOT file used to make new branches",
        required=True,
        )
    parser.add_argument(
        "-t", "--trees",
        help="Trees to merge",
        required=True,
        )
    parser.add_argument(
        "-m", "--mergingAlgo",
        help="Merging algorithm. Choices=['branches', 'events', 'efficiencies', 'efficiencies:<treeName/genWeightBranchName>', 'firstFound']",
        required=True,
        )
    parser.add_argument(
        "-o", "--outputfile",
        help="Output ROOT file name",
        required=True,
        )
    parser.add_argument(
        "-debug", "--debug",
        help="Debug mode",
        action="store_true",
        )
    
    
    args = parser.parse_args()

    main(make_files_list(args.inputfiles), args.trees.split(","), args.mergingAlgo.split(","), args.outputfile, args.debug)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


