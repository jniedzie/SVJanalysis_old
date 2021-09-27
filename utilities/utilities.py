from collections import Counter
import re
import json
import subprocess


def run_bash_command(bash_command):
    """Returns the output of a bash command.

    Args:
        bash_command (str)

    Returns:
        str
    """

    return subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")[:-1]


def make_file_list(files_arg):
    """Make list of ROOT files to merge.

    Args:
        files_arg (str): files argument from the argument parser.
            Comma separated ROOT file names e.g. file1.root,file2.root or
            text file name with a ROOT file name on each line.

    Returns:
        list[str]
    """

    # If the file argument has file name expansion (not really a regex, ls does not support regex)
    if "*" in files_arg or "!" in files_arg or "?" in files_arg or "^" in files_arg \
    or "[" in files_arg or "{" in files_arg:
        root_files = run_bash_command("ls %s 2>/dev/null" %files_arg).split("\n")
    # If list of ROOT files in txt file
    elif files_arg.endswith(".txt"):
        with open (files_arg, "r") as txt_file:
            root_files = [ x.replace("\n", "") for x in txt_file.readlines() ]
    # Else we assume coma separated list of ROOT files
    else:
        root_files = files_arg.split(",")

    # In case there were empty lines in txt file or extra comas, remove empty strings
    root_files = [ x for x in root_files if x != "" ]

    return root_files


def list2str(list_, str_for_concatenation=""):
    """Concatenate elements of a list into an str.

    Args:
        list_ (list): List of elements convertible to str
        str_for_concatenation (str, optional, default=""): Text to be placed
            in between concatenated elements

    Returns: 
        str

    """

    if len(list_)==0:
        text = ""
    else:
        text = list_[0]
        for el in list_[1:]: text = text + str_for_concatenation + str(el)

    return text


def most_common(list_):
    """Returns most common element of a list.

    Args:
        list_ (list[T])

    Returns:
        T
    """

    data = Counter(list_)
    return data.most_common(1)[0][0]


def inregex(name, regex_list):
    """
    Returns the list of the indices of the regexes which the text matches.
    An empty list is returned if the text matches none of the regexes.

    Args:
        name (str)
        regex_list (list): list or any type that can be casted to a list
    
    Returns:
        list[int]
    """

    if not isinstance(regex_list, list):
        regex_list = list(regex_list)
    indices = []
    for idx, regex in enumerate(regex_list):
        research = re.search(regex, name)
        if (research != None):
           indices.append(idx)

    return indices


def make_json_data(json_file_name):
    """
    Take as input the name of a json file and return json data as a dictionary
    where constants from the "CONSTANTS" key are replaced by their corresponding
    value. See doc string of makeDict.

    """

    with open(json_file_name, 'r') as f:
        json_data = json.load(f)
    return make_dict(json_data)


def make_dict(json_data_in, inplace=True):
    """
    Args:
        json_data_in (dict): Dictionary representing json data, e.g. from
            json.load(file_name), with a key "CONSTANTS" defining constants
            to be replaced in the json file.
            The key "CONSTANTS" will be removed from the dictionary.
            Constants must be written as {constant_name} in the rest of the json file.
        inplace (bool, optional, default=True): Modify json data inplace (True)
             or modify only a copy of the input dictionary (False)

    Returns:
        dict: modified copy of the input dictionary

    Examples:
        >>> json_data_in = {
                CONSTANTS": {"PATH": "/eos/cms/store/group/phys_exotica/svjets/t_channel"},
	        "files": ["{PATH}/1.root", "{PATH}/2.root"]
            }
        >>> make_dict(json_data_in)
        {'files': ['/eos/cms/store/group/phys_exotica/svjets/t_channel/1.root',
                   '/eos/cms/store/group/phys_exotica/svjets/t_channel/2.root']}
    """

    if "CONSTANTS" not in json_data_in.keys():
        print("WARNING: Cannot replace constants in json data because json data has no key \"CONSTANTS\".")
        print("         Returning json data as is.")
        return(json_data_in)

    if not inplace:
        json_data = json_data_in.copy()
    else:
        json_data = json_data_in

    constants = json_data["CONSTANTS"]
    json_data.pop("CONSTANTS", None)

    regexes = [ (re.compile(r"\{"+const+"\}"), constants[const]) for const in constants.keys() ]

    def replace(json_data):
        if isinstance(json_data, dict):
            for k, v in json_data.items():
                json_data[k] = replace(v)
            return json_data
        elif isinstance(json_data, list):
            for idx, v in enumerate(json_data):
                json_data[idx] = replace(v)
            return json_data
        elif isinstance(json_data, str):
            for regex, repl in regexes:
                json_data = re.sub(regex, repl, json_data)
            return json_data
        else:
            return json_data

    replace(json_data)
    return json_data
     

