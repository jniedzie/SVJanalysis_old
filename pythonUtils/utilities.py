import ROOT

def list2str(L, strForConcatenation=""):
    '''
    Takes a list of elements convertible to str and optionally a string.
    returns the concatenation of elements in the list separated by the optional string

    '''

    if len(L)==0:
        return("")
    else:
        s = L[0]
        for el in L[1:]: s = s + strForConcatenation + str(el)
        return(s)


