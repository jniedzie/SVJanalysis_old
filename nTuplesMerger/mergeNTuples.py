from coffea import processor
from coffea.nanoevents import BaseSchema
import awkward as ak
import uproot
import uproot3
import numpy as np
import time
import argparse
import sys

sys.path.append("../utilities/")
import uproot3Utilities as uproot3utl


def make_files_list(files_arg):
    """ ... """

    # If list of ROOT files in txt file
    if not files_arg.endswith(".root"):
        with open (files_arg, "r") as txt_file:
            root_files = [ x.replace("\n", "") for x in txt_file.readlines() ]
    # Else we assume coma separated list of ROOT files
    else:
        root_files = files_arg.split(",")

    return root_files


def make_branches_to_read(file_names, tree_name):
    """ ... """

    branches_to_read = {}
    branches_names_added = []
    for ifile, file_name in enumerate(file_names):
        file_ = uproot3.open(file_name)
        tree = file_[tree_name]
        branches_this_file = [ branch.decode("utf-8") for branch in tree.keys() ]
        branches_to_read[file_name] = [ branch_name for branch_name in branches_this_file if branch_name not in branches_names_added ]
        branches_names_added = branches_names_added + branches_to_read[file_name]

    return branches_to_read

def get_all_branches_name(file_name, tree_name):
    """ ... """

    tree = uproot3.open(file_name)[tree_name]
    branches_to_read = [ branch.decode("utf-8") for branch in tree.keys() ]

    return branches_to_read


def write_root_file(trees, merging_algos, output_file_name, debug):
    """Write histograms, stored in a coffea accumulator, to a ROOT file.

    Args:
        trees (dict[dict[awkward.highlevel.Array]])
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
            print("\nTTree %s saved to output file %s" %(tree_name, output_file_name))

    return


def main(input_files, tree_names, merging_algos, output_file, debug):

    print("Input files:")
    print(input_files)

    trees = {}
    for tree_name, merging_algo in zip(tree_names, merging_algos):

        if merging_algo == "branches":
            branches_to_read = make_branches_to_read(input_files, tree_name)
            branches = {}
            for ifile, file_name in enumerate(input_files):
                tree = uproot.open(file_name + ":" + tree_name)
                for branch_name in branches_to_read[file_name]:
                    branches[branch_name] = tree[branch_name].array()

            trees[tree_name] = { k: v for k, v in branches.items() }

        elif merging_algo == "events":
            branches_to_read = get_all_branches_name(input_files[0], tree_name)
            #branches_to_read = get_all_branches_name(input_files[0], tree_name)[-1100:-1000]
            #branches_to_read = ["nJet"]
            branches = {}
            for ifile, file_name in enumerate(input_files):
                tree = uproot.open(file_name + ":" + tree_name)
                for branch_name in branches_to_read:
                    if ifile == 0:
                        branches[branch_name] = [tree[branch_name].array()]
                    else:
                        branches[branch_name].append(tree[branch_name].array())
 
            trees[tree_name] = { k: ak.concatenate(v, axis=0) for k, v in branches.items() }

        elif merging_algo.startswith("efficiencies"):
            branches_to_read = get_all_branches_name(input_files[0], tree_name)
            weight_branch_name = merging_algo.replace("efficiencies:", "")

            branches = {}
            weight_branches = {}

            for ifile, file_name in enumerate(input_files):
                tree = uproot.open(file_name + ":" + tree_name)
                if ifile == 0:
                    weight_branches = [ak.sum(uproot.open(file_name + ":" + weight_branch_name).array(), axis=0)]
                else:
                    weight_branches.append(ak.sum(uproot.open(file_name + ":" + weight_branch_name).array(), axis=0))

                for branch_name in branches_to_read:
                    if ifile == 0:
                        branches[branch_name] = [tree[branch_name].array()]
                    else:
                        branches[branch_name].append(tree[branch_name].array())

            trees[tree_name] = { k: sum([w*b for w, b in zip(weight_branches, v)])/sum(weight_branches) for k, v in branches.items() }
            
        elif merging_algo == "firstFound":
            for ifile, file_name in enumerate(input_files):
                file_ = uproot.open(file_name)
                print(file_.keys())
                if tree_name+";1" in file_.keys():
                    print(file_[tree_name])
                    tree = file_[tree_name].arrays()
                    trees[tree_name] = { k: tree[k] for k in tree.fields }
                    break

        else:
            print("ERROR: Unknown merging algorithm: %s" %merging_algo)
            sys.exit()

    ## Making output ROOT file
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
        help="Merging algorithm. Choices=['branches', 'events', 'efficiencies:<treeName/genWeightBranchName>', 'firstFound']",
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


