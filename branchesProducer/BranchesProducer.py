import sys

from coffea import processor

sys.path.append("../utilities/")
import nameUtilities as nameutl
import coffeaUtilities as cfutl
import awkwardArrayUtilities as akutl
import physicsUtilities as phutl
import variablesComputation.awkwardArray.jetVariables as jetvars
import variablesComputation.awkwardArray.eventVariables as evtvars
import objects as obj


# TODO: Maybe these calculate_* functions can be moved to a separate script

def calculate_generalized_angularities(branches, jet_collection_name, jet_radius, jets, njets, jet_filter, jet_pf_cands, sum_pfcands_pt):
    """Compute generalized angularities.

    Args:
        branches (dict[str, awkward.highlevel.Array])
        jet_collection_name (str)
        jet_radius (float)
        jets (awkward.highlevel.Array)
        njets (awkward.highlevel.Array)
        jet_filter (awkward.highlevel.Array)
        jet_pf_cands (awkward.highlevel.Array)
        sum_pfcands_pt (awkward.highlevel.Array)

    Returns:
        None
    """

    branches[jet_collection_name + "_ptD"] = jetvars.calculate_ptD(jet_pf_cands, jet_filter=jet_filter, sum_pfcands_pt=sum_pfcands_pt)
    branches[jet_collection_name + "_LHA"] = jetvars.calculate_lha(jet_pf_cands, jets, jet_radius, jet_filter=jet_filter, njets=njets, sum_pfcands_pt=sum_pfcands_pt)
    branches[jet_collection_name + "_girth"] = jetvars.calculate_girth(jet_pf_cands, jets, jet_radius, jet_filter=jet_filter, njets=njets, sum_pfcands_pt=sum_pfcands_pt)
    branches[jet_collection_name + "_thrust"] = jetvars.calculate_thrust(jet_pf_cands, jets, jet_radius, jet_filter=jet_filter, njets=njets, sum_pfcands_pt=sum_pfcands_pt)
    branches[jet_collection_name + "_multiplicity"] = jetvars.calculate_multiplicity(jet_pf_cands, jet_filter=jet_filter)


def calculate_axes(branches, jet_collection_name, jets, njets, jet_filter, jet_pf_cands):
    """Compute jet axis minor, major and average.

    Args:
        branches (dict[str, awkward.highlevel.Array])
        jet_collection_name (str)
        jets (awkward.highlevel.Array)
        njets (awkward.highlevel.Array)
        jet_filter (awkward.highlevel.Array)
        jet_pf_cands (awkward.highlevel.Array)

    Returns:
        None
    """

    axis_major, axis_minor, axis_avg = jetvars.calculate_axes(jet_pf_cands, jets, jet_filter=jet_filter, njets=njets)
    branches[jet_collection_name + "_axisMajor"] = axis_major
    branches[jet_collection_name + "_axisMinor"] = axis_minor
    branches[jet_collection_name + "_axisAvg"] = axis_avg

    del axis_major, axis_minor, axis_avg


def calculate_energy_fractions(branches, jet_collection_name, jets, njets, jet_filter, jet_pf_cands):
    """Compute jet energy fractions.

    Args:
        branches (dict[str, awkward.highlevel.Array])
        jet_collection_name (str)
        jets (awkward.highlevel.Array)
        njets (awkward.highlevel.Array)
        jet_filter (awkward.highlevel.Array)
        jet_pf_cands (awkward.highlevel.Array)

    Returns:
        None
    """

    branches[jet_collection_name + "_chHEF"] = jetvars.calculate_chHEF(jet_pf_cands, jets, jet_filter=jet_filter, njets=njets)
    branches[jet_collection_name + "_neHEF"] = jetvars.calculate_neHEF(jet_pf_cands, jets, jet_filter=jet_filter, njets=njets)
   

def calculate_ecfs(branches, jet_collection_name, beta, jet_filter, jet_pf_cands, sum_pfcands_pt):
    """Compute jet energy correlation functions.

    Args:
        branches (dict[str, awkward.highlevel.Array])
        jet_collection_name (str)
        beta (float)
        jet_filter (awkward.highlevel.Array)
        jet_pf_cands (awkward.highlevel.Array)
        sum_pfcands_pt (awkward.highlevel.Array)

    Returns:
        None
    """

    # Some variables computation are commented out because it takes too much time to compute and
    # the code is kept for reference

    e2b1 = jetvars.calculate_ecf_e2(jet_pf_cands, beta, jet_filter=jet_filter, sum_pfcands_pt=sum_pfcands_pt)
    #e3b1 = jetvars.calculate_ecfs_e3(jet_pf_cands, beta, sum_pfcands_pt=sum_pfcands_pt, calculate_ecfgs=False)
    v1e3b1, v2e3b1, e3b1 = jetvars.calculate_ecfs_e3(jet_pf_cands, beta, jet_filter=jet_filter, sum_pfcands_pt=sum_pfcands_pt, calculate_ecfgs=True)
    #e4b1 = jetvars.calculate_ecfs_e4(jet_pf_cands, beta, sum_pfcands_pt=sum_pfcands_pt, calculate_ecfgs=False)
    #v1e4b1, v2e4b1, v3e4b1, v4e4b1, v5e4b1, e4b1 = jetvars.calculate_ecfs_e4(jet_pf_cands, beta, sum_pfcands_pt=sum_pfcands_pt, calculate_ecfgs=True)
    branches[jet_collection_name + "_e2b1"] = e2b1
    branches[jet_collection_name + "_v1e3b1"] = v1e3b1
    branches[jet_collection_name + "_v2e3b1"] = v2e3b1
    branches[jet_collection_name + "_e3b1"] = e3b1
    #branches[jet_collection_name + "_v1e4b1"] = v1e4b1
    #branches[jet_collection_name + "_v2e4b1"] = v2e4b1
    #branches[jet_collection_name + "_v3e4b1"] = v3e4b1
    #branches[jet_collection_name + "_v4e4b1"] = v4e4b1
    #branches[jet_collection_name + "_v5e4b1"] = v5e4b1
    #branches[jet_collection_name + "_e4b1"] = e4b1
    branches[jet_collection_name + "_c2b1"] = jetvars.calculate_ecf_c(1., e2b1, e3b1)
    branches[jet_collection_name + "_d2b1"] = jetvars.calculate_ecf_d(1., e2b1, e3b1)
    branches[jet_collection_name + "_m2b1"] = jetvars.calculate_ecf_m(e2b1, v1e3b1)
    branches[jet_collection_name + "_n2b1"] = jetvars.calculate_ecf_n(e2b1, v2e3b1)

    del e2b1, v1e3b1, v2e3b1, e3b1


def calculate_efps(branches, jet_collection_name, efp_degree, jet_filter, jet_pf_cands):
    """Compute jet energy flow polynomials.

    Args:
        branches (dict[str, awkward.highlevel.Array])
        jet_collection_name (str)
        efp_degree (float)
        jet_filter (awkward.highlevel.Array)
        jet_pf_cands (awkward.highlevel.Array)

    Returns:
        None
    """

    efps = jetvars.calculate_efps(jet_pf_cands, efp_degree, jet_filter=jet_filter)
    # Not storing EFP0, which is always 1, on purpose
    for idx, efp in enumerate(efps[1:]):
        branches[jet_collection_name + "_efp%dd%d" %(idx+1, efp_degree)] = efp

    del efps


def calculate_azimuthal_angle_variables(branches, jet_collection_name, jets, njets, met):
    """Compute azimuthal angle between jets and met.

    Args:
        branches (dict[str, awkward.highlevel.Array])
        jet_collection_name (str)
        jets (awkward.highlevel.Array)
        njets (awkward.highlevel.Array)
        met (awkward.highlevel.Array)

    Returns:
        None
    """

    jets_delta_phi = jetvars.calculate_delta_phi(jets, met)
    branches[jet_collection_name + "_deltaPhi"] = jets_delta_phi
    branches["DeltaPhiMin" + jet_collection_name] = evtvars.calculate_delta_phi_min(jets_delta_phi, njets=njets)

    del jets_delta_phi


def calculate_razor_variables(branches, jet_collection_name, jets, njets, met):
    """Compute Razor variables.

    Args:
        branches (dict[str, awkward.highlevel.Array])
        jet_collection_name (str)
        jets (awkward.highlevel.Array)
        njets (awkward.highlevel.Array)
        met (awkward.highlevel.Array)

    Returns:
        None
    """

    mr = evtvars.calculate_razor_MR(jets, njets=njets)
    mrt = evtvars.calculate_razor_MRT(jets, met, njets=njets)
    branches["RazorMR" + jet_collection_name] = mr
    branches["RazorMRT" + jet_collection_name] = mrt
    branches["RazorR" + jet_collection_name] = evtvars.calculate_razor_R(mr, mrt)

    del mr, mrt


def calculate_transverse_masses(branches, jet_collection_name, jets, njets, met):
    """Compute transverse mass for different pairs of jets.

    Args:
        branches (dict[str, awkward.highlevel.Array])
        jet_collection_name (str)
        jets (awkward.highlevel.Array)
        njets (awkward.highlevel.Array)
        met (awkward.highlevel.Array)

    Returns:
        None
    """

    branches["MT01" + jet_collection_name] = evtvars.calculate_MT(jets, met, jet_indices=[0,1], njets=njets)
    branches["MT02" + jet_collection_name] = evtvars.calculate_MT(jets, met, jet_indices=[0,2], njets=njets)
    branches["MT12" + jet_collection_name] = evtvars.calculate_MT(jets, met, jet_indices=[1,2], njets=njets)


def calculate_branches(jet_collection_name, jet_radius, jets, njets, jet_filter, jet_pf_cands, sum_pfcands_pt, met):
    """Compute all required new branches.

    Args:
        jet_collection_name (str)
        jet_radius (float)
        jets (awkward.highlevel.Array)
        njets (awkward.highlevel.Array)
        jet_filter (awkward.highlevel.Array)
        jet_pf_cands (awkward.highlevel.Array)
        sum_pfcands_pt (awkward.highlevel.Array)

    Returns:
        None
    """

    beta = 1        # Angular exponent for ECFs
    efp_degree = 3  # EFP graphs degree

    # Book histogram to store branches
    branches = {}

    # Always need to have the branch storing the size of jagged arrays
    branches["n" + jet_collection_name] = njets

    # Jet variables
    calculate_generalized_angularities(branches, jet_collection_name, jet_radius, jets, njets, jet_filter, jet_pf_cands, sum_pfcands_pt)
    calculate_axes(branches, jet_collection_name, jets, njets, jet_filter, jet_pf_cands)
    calculate_energy_fractions(branches, jet_collection_name, jets, njets, jet_filter, jet_pf_cands)
    calculate_ecfs(branches, jet_collection_name, beta, jet_filter, jet_pf_cands, sum_pfcands_pt)
    calculate_efps(branches, jet_collection_name, efp_degree, jet_filter, jet_pf_cands)
    calculate_azimuthal_angle_variables(branches, jet_collection_name, jets, njets, met)

    # Event variables depending on the jet flavor
    calculate_razor_variables(branches, jet_collection_name, jets, njets, met)
    calculate_transverse_masses(branches, jet_collection_name, jets, njets, met)

    return branches


def accumulate(branches, output):
    """Accumulate branches jagged arrays into coffea accumulator.
    
    Args:
        branches (dict[str, awkward.highlevel.Array])
        output (coffea.processor.dict_accumulator)

    Returns:
        None
    """

    for branch_name, branch in branches.items():
        output[branch_name] = cfutl.accumulate(branch)




class BranchesProducer(processor.ProcessorABC):
    """Processor for computing new branches to be added to a new ROOT file.

    The class must implement the following methods:
        * __init__: Define the branches to be computed
        * process: Calculation of the new quantities and accumulation into
                   the branches.
    """


    def __init__(self):
        """Define the branches to be computed.

        New branches to be computed are defined in the class attribute
        _accumulator (coffea.processor.dict_accumulator).

        Args:
            None
        
        Returns:
            None
        """

        self.jet_types = ["AK4", "AK8"]  # capitals used because capitals are used in PFNanoAOD

        ## Define accumulator
        #  By using dict_accumulator, the branches type do not have to be defined
        self._accumulator = processor.dict_accumulator({})


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        """Compute new quantities to be added to the ROOT NTuple.

        Args:
            events (awkward.Array): the Events TTree opened with uproot.

        Returns:
            coffea.processor.dict_accumulator
        """

        ## Define accumulator
        output = self.accumulator.identity()

        output["genWeight"] = cfutl.accumulate(events.genWeight)

        met = obj.get_met(events, "MET")

        for jet_type in self.jet_types:

            jet_collection_name = nameutl.jet_algo_name_to_jet_collection_name(jet_type)
            jet_radius = phutl.jet_name_to_jet_radius(jet_type)
            jets = obj.get_jets(events, jet_type)
            jet_pf_cands = obj.get_pf_cands(events, jet_type)
            njets = obj.get_collection_size(events, jet_collection_name)
            sum_pfcands_pt = jetvars.calculate_sum_pfcands_pt(jet_pf_cands, jagged=False)

            # Defining a jet_filter to make irregular ak array
            # e.g. for njets = [[1, 2], [1, 2, 3], [1]]
            # jet_filter = [[True, True, False], [True, True, True], [True, False, False]]
            jet_filter = jetvars.make_jet_filter(njets)

            # Calculate all required new branches...
            branches = calculate_branches(jet_collection_name, jet_radius, jets, njets, jet_filter, jet_pf_cands, sum_pfcands_pt, met)
            # ... and write them in coffea accumulator
            accumulate(branches, output)
 
            # Freeing the memory to avoid memory leaks
            del jets
            del jet_pf_cands
            del njets
            del sum_pfcands_pt
            del jet_filter
            del branches

        del met

        return output


    def postprocess(self, accumulator):
        return accumulator

