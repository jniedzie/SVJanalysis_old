from coffea import hist
from coffea import processor
import awkward as ak
import sys

sys.path.append("../pythonUtils/")
import utilities as utl


class HistogramDefault(processor.ProcessorABC):
    """Parent class for coffea processors accumulating histograms.

    Children classes must implement the following methods:
        * __init__
        * initialize_ak_arrays
        * initialize_gen_weights

    __init__ must define the following attributes:
        * self.variables
        * self.gen_weights_info
        * self.cuts
        * self._accumulator
        * self.file_type (will soon be obsolete when using 106X PFNanoAOD format only)

    initialize_ak_arrays must return the following objects:
        * var_arrays (dict): arrays for all variables
        * masks (dict): masks used for the necessary cuts to define some variables

    initialize_gen_weights must return the following objects:
        * gen_weights(dict): generator weights for all variables
    """


    def __init__(self):
        pass


    def check_dict_keys(self, dict1, dict2, dict1_name, dict2_name, error_keys_dict1_not_in_keys_dict2=True, error_keys_dict2_not_in_keys_dict1=True):
        """Check if dict1 and dict2 have the same keys. Interrupt program if not the case."""

        if isinstance(dict1, list):
            dict1 = { k: "" for k in dict1 }
        if isinstance(dict2, list):
            dict1 = { k: "" for k in dict2 }

        if error_keys_dict1_not_in_keys_dict2:
            message_level1 = "ERROR"
        else:
            message_level1 = "WARNING"
        if error_keys_dict2_not_in_keys_dict1:
            message_level2 = "ERROR"
        else:
            message_level2 = "WARNING"

        error = False
        dict1_not_none_keys = sorted([ k for k in dict1.keys() if dict1[k] is not None ])
        dict2_not_none_keys = sorted([ k for k in dict2.keys() if dict2[k] is not None ])
        if dict1_not_none_keys != dict2_not_none_keys:
            dict1_keys_not_in_dict2 = [ k for k in dict1_not_none_keys if k not in dict2_not_none_keys ]
            dict2_keys_not_in_dict1 = [ k for k in dict2_not_none_keys if k not in dict1_not_none_keys ]
            if len(dict1_keys_not_in_dict2) > 0:
                if error_keys_dict2_not_in_keys_dict1: error = True
                print("\n%s: Some %s have been defined but do not have a defined %s" %(message_level2, dict1_name, dict2_name))
                print("Missing %s:" %dict2_name)
                print(dict1_keys_not_in_dict2)
            if len(dict2_keys_not_in_dict1) > 0:
                if error_keys_dict1_not_in_keys_dict2: error = True
                print("\n%s: Some %s have a defined %s but have not been defined" %(message_level1, dict1_name, dict2_name))
                print("Missing %s:" %dict1_name)
                print(dict2_keys_not_in_dict1)
            if error:
                sys.exit()
            else:
                return False

        return True
     

    def make_variables_names(self, variables_descriptions):
        """
        Make variables based on a list of templates associated with a list of
        parameters to be replaced by a given value.

        Args:
            variables_descriptions (list[tuple[str, list[dict]]])

        Returns:
            list[str]

        Examples:
            >>> variables_descriptions =
                [ ( "{jet}Jet_n",  [ {"jet": "ak4"}, {"jet": "ak8"} ] ),
                    "{jet}Jet_pt", [ {"jet": "ak4"} ] ),
                ]
            >>> make_variables_names(variables_descriptions)
            ["ak4Jet_n", "ak8Jet_n", "ak4Jet_pt"]
        """

        variables = []

        for entry in variables_descriptions:
            variable_template = entry[0]
            list_of_params_dict = entry[1]
            for params_dict in list_of_params_dict:
                var = variable_template
                for k, v in params_dict.items():
                    var = var.replace("{"+k+"}", v)
                variables.append(var)

        return variables



    def get_binning(self):
        """Return binning for all variables.

        Form of self.binning:
            {
		"noregex": {
		    "ak8Jet1_ak8Jet2_mass"  :  [500  ,  0   , 5000 ],
		},
		"regex": {
		    "ak[48]Jet_n"           :  [20   ,  0   , 20   ],
		}
	    }
        """

        self.binning = {}

        # Loop over all variables to histogram
        for variable in self.variables:

            # List of regexes for binning definition
            regexes = list(self.binning_info["regex"].keys())

            if variable not in self.binning_info["noregex"]:
                indices = utl.inregex(variable, regexes)
                if len(indices) == 0:
                    print("ERROR: Binning of %s is not defined." %variable)
                    sys.exit()
                elif len(indices) > 1:
                    print("ERROR: %s matches several regexes. Binning cannot be defined." %variable)
                    print("regexes matched:")
                    for idx in indices: print("\t%s" %(regexes[idx]))
                    sys.exit()
                else: 
                    self.binning[variable] = self.binning_info["regex"][regexes[indices[0]]]
            else:
                self.binning[variable] = self.binning_info["noregex"][variable]

        return


    def define_histograms(self):
        """Define a dict with coffea histograms."""

        self.histograms = {}

        # Loop over all variables to histogram
        for variable in self.variables:
            # x axis label
            label = variable
            # Define coffea histogram for the variable
            binning = self.binning[variable]
            self.histograms[variable] = hist.Hist(variable, hist.Bin("axis", label, binning[0], binning[1], binning[2]))

        return


    def fill_histogram(self, output, variable, axis, weight):
        """Fill histogram if the variable has been defined."""

        type1 = str(ak.type(axis))
        type2 = str(ak.type(weight))
        type1 = utl.list2str(type1.split(" * ")[:-1], " * ")
        type2 = utl.list2str(type2.split(" * ")[:-1], " * ")
        if variable in output.keys():
            if type1 == type2:
                output[variable].fill(axis=axis, weight=weight)
            else:
                print("ERROR: Array and weight for variable %s have different sizes:" %variable)
                print("array  size: %s" %(type1))
                print("weight size: %s" %(type2))
                sys.exit()
        else:
            print("ERROR: Variable %s is not defined in the coffea accumulator!" %variable)
            sys.exit()
     
        return


    def fill_histograms(self, var_arrays, gen_weights, output):
        """Fill all histograms.
  
        Args:
            var_arrays (dict[str, ak.Array])
            gen_weights (dict[str, ak.Array])
            output (coffea.processor.dict_accumulator)

        Returns:
            None
        """

        size2gen = {}
        for key in gen_weights.keys():
            size = str(ak.size(gen_weights[key]))
            if size not in size2gen.keys():
                size2gen[size] = [key]
            else:
                size2gen[size].append(key)

        if len(size2gen.keys()) != len(gen_weights.keys()):
            print("\nWARNING: Some gen weights have the same length, automatic matching between variable and gen_weight arrays may be incorrect.")
            print("         Check the cutflow to see whether these variables have 100% efficiency.")
            print("         If not, the gen weights for the relevant variables can be defined in gen_weights_info.")
            print("         The following gen_weights have the same lengths:")
            for key, value in size2gen.items():
                if len(value) > 1:
                    print("%s with length %s" %(value, key))

        for variable in self.variables:
            # Find corresponding gen weights
            size = str(ak.size(var_arrays[variable]))
            if (size not in size2gen.keys()) and (variable not in self.gen_weights_info.keys()):
                print("\nWARNING: No gen weight for variable %s. Skipping." %variable)
            else:
                if variable in self.gen_weights_info.keys():
                    gen_weight = gen_weights[self.gen_weights_info[variable]]
                else:
                    gen_weight = gen_weights[size2gen[size][0]]
                
                ## For ECF, -1 means that the ECF could not be computed
                #  Do not add the -1 to the histograms, because this will
                # bias efficiency and ROC AUC calculations 
                if variable.endswith("_n2b1") or  variable.endswith("_n3b1"):
                    filter_ = var_arrays[variable] > -0.1
                    var_arrays[variable] = var_arrays[variable][filter_]
                    gen_weight = gen_weight[filter_]

                self.fill_histogram(output, variable, var_arrays[variable], gen_weight)


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        """Make arrays and gen_weights for all interesting variable, and fill histograms."""

        # Define accumulator
        output = self.accumulator.identity()

        # Make var_arrays that will be used to fill histograms
        # Also return masks used, will be useful for making gen weights
        var_arrays, masks = self.initialize_ak_arrays(events)

        # Running some sanity checks
        check_passed = self.check_dict_keys(self.variables, var_arrays, "variable", "ak array", error_keys_dict1_not_in_keys_dict2=False)
        if not check_passed:
            print("\nWARNING: Only variable defined in constructor will be histogrammed.")

        # Make gen weights
        gen_weights = self.initialize_gen_weights(events, masks)

        # Calculate the sum of gen weights for the cuts defined in contructors
        sum_gen_weights = {}
        for cut in self.cuts:
            sum_gen_weights[cut] = ak.sum(gen_weights[cut])

        # Save cutflow information
        for cut in self.cuts:
            output["cutflow"][cut] += sum_gen_weights[cut]

        # Fill histograms
        self.fill_histograms(var_arrays, gen_weights, output)

        return output


    def postprocess(self, accumulator):
        return accumulator


