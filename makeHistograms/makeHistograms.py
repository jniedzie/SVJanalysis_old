import coffea as cf
import awkward as ak
import awkward0 as ak0
import uproot
import uproot3
import numpy as np
import time
import argparse
import os
import sys
import json

sys.path.append("../utilities/")
import utilities as utl
import histogramProcessors as histogram_processors


def get_cross_section(file_name, branch_name="Metadata/GenCrossSection"):
    """Return sample cross section written in ROOT file.

    Args:
        file_name (str): ROOT file name
        branch_name (str): Branch name with cross section

    Returns:
        float
    """

    return uproot.open(file_name + ":" + branch_name).array()[0]


def get_n_generated_events(file_names, branch_name="Cutflow/NoCut"):
    """Return total number of generated events for a list of ROOT files.

    Args:
        file_names (list[str]): List of ROOT file names
        branch_name (str): Branch name with number of generated events

    Returns:
        float
    """

    n_generated_events = 0
    for file_name in file_names:
        n_generated_events += uproot.open(file_name + ":" + branch_name).array()[0]

    return n_generated_events


## TODO: This function should be extracted in utilities
def print_cutflow(cutflow):
    """Print cut efficiencies for cuts needed for defining some of the variables.

    For instance, to compute deltaR between leading 2 jets, events needs
    to have at least 2 jets.

    Args:
        cutflow (dict[str, coffea.processor.defaultdict_accumulator(float)]):
            Keys are cut names
            Values are numbers of processed events (sum of gen weights)
 
    Returns:
        None
    """

    len_column1 = max([ len(k) for k in cutflow.keys() ])

    print("\nCutflow:")
    print("\tCut" + (len_column1-3)*" " + "  Abs. eff. [%]")
    nAll = cutflow["noCut"]
    for cut in cutflow.keys():
        print("\t%s%s  %.2f" %(cut, (len_column1-len(cut))*" ", 100*cutflow[cut]/nAll))

    return


# TODO: This function should be added to uproot3Utilities
def write_root_file(accumulator, output_file_name, lumi, cross_section, n_generated_events):
    """Write histograms, stored in a coffea accumulator, to a ROOT file.

    Number of generated events (cutflow), pre-selection cut efficiency (efficiency), luminosity (lumi)
    and cross-section (xSection) are needed to normalize the histograms.

    Args:
        accumulator (dict[str, coffea.hist.Hist])
        output_file_name (str): Name of the ROOT file to create
        lumi (float): integrated luminosity
        cross_section (float): Cross section
        n_generated_events (float): Number of generated events (sum of generator weights)

    Returns:
        None
    """

    # Need to use uproot3 because it is not implemented yet in uproot4 (Feb. 2021)
    print("\nWill recreate ROOT file %s\n" %output_file_name)

    written_hists = []
    with uproot3.recreate(output_file_name) as f:
        for variable, hist in accumulator.items():
            if isinstance(hist, cf.hist.Hist):
                hist.scale(lumi * cross_section / n_generated_events )
                f[variable] = cf.hist.export1d(hist)
                written_hists.append(variable)

    print("Histograms written to output ROOT file:")
    for written_hist in written_hists:
        print("\t%s" %written_hist)

    return


def main(input_files, output_file, binning, processor, chunksize, maxchunks, nworkers, lumi):
    """Produce a ROOT file filled with histograms.
    
    Produce a ROOT file filled with histograms read from input ROOT files or computed on the fly.

    Args:
        input_files (list)
        output_file (str)
        binning (dict): Binning information
        processor (str)
        chunksize (int)
        maxchunks (int)
        nworkers (int)
        lumi (float)

    Returns:
        None
    """

    # Make histograms
    fileset = { "fileset": input_files }
    output = cf.processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = getattr(histogram_processors, processor)(binning),
        executor = cf.processor.iterative_executor,
        executor_args = {"schema": cf.nanoevents.BaseSchema, "workers": nworkers},
        chunksize = chunksize,
        maxchunks = maxchunks,
    )

    # Get sample cross-section
    cross_section = get_cross_section(input_files[0])

    # Get number of generated events
    n_generated_events = get_n_generated_events(input_files)

    # Print out cutflow
    cutflow = output.pop("cutflow")
    print_cutflow(cutflow)

    # Save histograms to ROOT file
    write_root_file(output, output_file, lumi, cross_section, n_generated_events)

    return


if __name__ == "__main__":

    tstart = time.time()

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputFiles",
        help="Input file",
        )
    parser.add_argument(
        "-o", "--outputFile",
        help="Output file",
        )
    parser.add_argument(
        "-b", "--binning",
        help="json file describing binning of the histograms",
        required=True
        )
    parser.add_argument(
        "-p", "--processor",
        help="Coffea processor to be used, as defined in processorHistogram.py",
        required=True
        )
    parser.add_argument(
        "-c", "--chunksize",
        help="Size of the data chunks (default=100000)",
        default=100000, type=int
        )
    parser.add_argument(
        "-m", "--maxchunks",
        help="Maximum number of chunks to process, no flag means no maximum",
        type=int
        )
    parser.add_argument(
        "-n", "--nworkers",
        help="Number of worker nodes (default=4)",
        default=4, type=int
        )
    parser.add_argument(
        "-l", "--lumi",
        help="Total luminosity for normalization of the histograms (default=59725 pb-1)",
        default=59725.0,   # 21071.0+38654.0
        type=float,
        )

    args = parser.parse_args()

    input_files = utl.make_file_list(args.inputFiles)
    
    # Define the binning of the different variables to histogram
    with open(args.binning, 'r') as f:
        binning = json.load(f)["binning"]

    # Make histograms
    main(input_files, args.outputFile, binning, args.processor, args.chunksize, args.maxchunks, args.nworkers, args.lumi)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


