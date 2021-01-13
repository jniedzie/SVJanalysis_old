import re


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

