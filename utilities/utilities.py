import re
import json


def list2str(L, strForConcatenation=""):
    """
    Takes a list of elements convertible to str and optionally a string.
    Returns the concatenation of elements in the list separated by the optional string.

    """

    if len(L)==0:
        return("")
    else:
        s = L[0]
        for el in L[1:]: s = s + strForConcatenation + str(el)
        return(s)


def inregex(name, regexList):
    """
    Takes a string text and either a regex or a list regexes.
    Returns the list of the indices of the regexes which the text matches.
    An empty list is returned if the text matches none of the regexes.
    
    """

    if type(regexList) != list:
        regexList = [regexList]
    indices = []
    for idx, regex in enumerate(regexList):
        research = re.search(regex, name)
        if (research != None):
           indices.append(idx)
    return(indices)


def makeJsonData(jsonFileName):
    """
    Take as input the name of a json file and return json data as a dictionary
    where constants from the "CONSTANTS" key are replaced by their corresponding
    value. See doc string of makeDict.

    """

    with open(jsonFileName, 'r') as f:
        jsonData = json.load(f)
    return(makeDict(jsonData))


def makeDict(jsonDataIn, inplace=True):
    """
    Take as input a dictionary representing json data, e.g. from json.load(filename)
    with a key "CONSTANTS" defining constants to be replaced in the json file.
    The key "CONSTANTS" will be removed from the dictionary.
    Constants must be written as {constantName} in the rest of the json file.
    The original dictionary is changed inplace by default. The original dictionary
    given as inputA modified copy of the is untouched if setting inplace=False.
    Return the modified dictionary in both cases (both inplace=True or False).

    Example:
	{"CONSTANTS": {
	    "PATH": "/eos/cms/store/group/phys_exotica/svjets/t_channel"
	    },
	 "files": ["{PATH}/1.root",
                   "{PATH}/2.root"]
        }
    will be changed into:
	{"files": ["/eos/cms/store/group/phys_exotica/svjets/t_channel/1.root",
                   "/eos/cms/store/group/phys_exotica/svjets/t_channel/2.root"]
        }

    """

    if "CONSTANTS" not in jsonDataIn.keys():
        print("WARNING: Cannot replace constants in json data because json data has no key \"CONSTANTS\".")
        print("         Returning json data as is.")
        return(jsonDataIn)

    if not inplace:
        jsonData = jsonDataIn.copy()
    else:
        jsonData = jsonDataIn

    constants = jsonData["CONSTANTS"]
    jsonData.pop("CONSTANTS", None)

    regexes = [ (re.compile(r"\{"+const+"\}"), constants[const]) for const in constants.keys() ]

    def replace(jsonData):
        if isinstance(jsonData, dict):
            for k, v in jsonData.items():
                jsonData[k] = replace(v)
            return(jsonData)
        elif isinstance(jsonData, list):
            for idx, v in enumerate(jsonData):
                jsonData[idx] = replace(v)
            return(jsonData)
        elif isinstance(jsonData, str):
            for regex, repl in regexes:
                jsonData = re.sub(regex, repl, jsonData)
            return(jsonData)
        else:
            return(jsonData)

    replace(jsonData)
    return(jsonData)
     

