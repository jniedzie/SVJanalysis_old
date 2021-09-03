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


def calculate_efficiency(fileset, efficiency_branch, gen_weight_branch="Events/genWeight"):
    """
    Calculate efficiency using the effienciency written in input ROOT files,
    weighted by the sum of generator weights for the events in this ROOT file.

    Args:
        fileset (list[str]): list of ROOT files
        efficiency_branch (str): name of the efficiency branch
        gen_weight_branch (str): name of generator weights branch

    Returns:
        float
    """

    efficiencies = []
    sum_gen_weights = []

    for file_ in fileset:
        data = uproot.open(file_)
        efficiencies.append( data[efficiency_branch].arrays()[efficiency_branch.split("/")[-1]][0] )
        sum_gen_weights.append( ak.sum(data[gen_weight_branch].arrays()[gen_weight_branch.split("/")[-1]]) )

    efficiency = sum([ eff*sgw for eff, sgw in zip(efficiencies, sum_gen_weights)]) / sum(sum_gen_weights)
    
    return efficiency


# TODO: This function should be added to uproot3Utilities
def write_root_file(accumulator, output_file_name, n_processed_events, lumi, xsection, efficiency):
    """Write histograms, stored in a coffea accumulator, to a ROOT file.

    Number of generated events (cutflow), pre-selection cut efficiency (efficiency), luminosity (lumi)
    and cross-section (xSection) are needed to normalize the histograms.

    Args:
        accumulator (dict[str, coffea.hist.Hist])
        output_file_name (str): Name of the ROOT file to create
        n_generated_events (float): Number of processed events (sum of gen weights)
        lumi (float): integrated luminosity
        xsection (float): Cross section
        efficiency (float): Pre-selection cut efficiency 

    Returns:
        None
    """

    # Need to use uproot3 because it is not implemented yet in uproot4 (Feb. 2021)
    print("\nWill recreate ROOT file %s\n" %output_file_name)

    written_hists = []
    with uproot3.recreate(output_file_name) as f:    # Create/update output file
        for variable, hist in accumulator.items():
            if isinstance(hist, cf.hist.Hist):
                hist.scale(lumi * xsection * efficiency / n_processed_events )
                f[variable] = cf.hist.export1d(hist)
                written_hists.append(variable)

    print("Histograms written to output ROOT file:")
    for written_hist in written_hists:
        print("\t%s" %written_hist)

    return


def main(input_files, output_file, binning, processor, chunksize, maxchunks, nworkers, lumi, use_efficiencies, efficiency):
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
        use_efficiencies (bool)
        efficiency (str or float): Name of the efficiency branch or efficiency

    Returns:
        None
    """

    ## Get sample name and cross-section
    xsection = 1.  # TODO: Must be read from Metadata tree


    ## Fileset
    fileset = { "fileset": input_files }


    ## Make histograms
    output = cf.processor.run_uproot_job(
        fileset,
        treename = "Events",
        processor_instance = getattr(histogram_processors, processor)(binning),
        executor = cf.processor.iterative_executor,
        executor_args = {"schema": cf.nanoevents.BaseSchema, "workers": nworkers},
        chunksize = chunksize,
        maxchunks = maxchunks,
    )

    ## Efficiencies
    if use_efficiencies:
        if isinstance(efficiency, str):
            efficiency = calculate_efficiency(sample["fileset"], efficiency)
        # else efficiency is a float that already corresponds to the efficiency
    else:
        efficiency = 1.


    ## Print out cutflow
    cutflow = output.pop("cutflow")
    print_cutflow(cutflow)


    ## Save histograms to ROOT file
    write_root_file(output, output_file, cutflow["noCut"], lumi, xsection, efficiency)

    return


if __name__ == "__main__":

    tstart = time.time()

    ## Parse arguments
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
        )
    parser.add_argument(
        "-e", "--efficiency",
        help="Use efficiencies from the input files",
        nargs='?', const="Efficiencies/totalEfficiency", default="False",
        )

    args = parser.parse_args()

    
    ## Define the binning of the different variables to histogram
    with open(args.binning, 'r') as f:
        binning = json.load(f)["binning"]

    ## TODO: Read info from Cutflow and Metadata tree
    ## Efficiencies
    if args.efficiency == "False":
        use_efficiencies = False
        efficiency = None
    else:
        use_efficiencies = True
        if args.efficiency.replace(".", "").isdigit():
            efficiency = float(args.efficiency)
        else:
            efficiency = args.efficiency

    input_files = utl.make_file_list(args.inputFiles)


    # Make histograms
    main(input_files, args.outputFile, binning, args.processor, args.chunksize, args.maxchunks, args.nworkers, float(args.lumi), use_efficiencies, efficiency)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)


