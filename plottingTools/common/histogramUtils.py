import math


def get_histogram(tfile, histName):
    """
    Return histogram from ROOT file

    Parameters
    ----------
    tfile   : ROOT.TFile
    histName: str
    """

    h = tfile.Get(histName)
    if not h:
        raise Exception("Failed to load histogram {}.".format(histName))
    return h


def get_hist_sum(hist, range_=None):
    """
    Return sum of all or some bins of a TH1 histogram.

    Parameters
    ----------
    hist: ROOT.TH1D
    range_: range
        If None, then sum of all bins, including underflow and overflow bins,
        is returned

    Returns
    -------
    histSum: float
    """

    histSum = 0

    if range_ is None:
        xaxis = histogram.GetXaxis()
        nbins = xaxis.GetNbins()
        range_ = range(0, nbins+2)

    for ibin in range_:
        histSum += hist.GetBinContent(ibin)

    return histSum


def calc_cumulative_sums(histogram):
    """
    Calculate cumulative sum of a TH1 histogram.
   
    Parameters
    ----------
    histogram: ROOT.TH1D
    
    Returns
    -------
    cumulativeSumLower: list < float >
    cumulativeSumUpper: list < float >
    histSum: float
    """

    # Calculate number of bins (GetNbins returns number of bins except underflow and overflow bins)
    xaxis = histogram.GetXaxis()
    nbins = xaxis.GetNbins()
    # Calculate < cut cumulative sum
    cumulativeSumLower = [ histogram.GetBinContent(0) ]
    for ibin in range(1, nbins+2):
        cumulativeSumLower.append(cumulativeSumLower[-1] + histogram.GetBinContent(ibin))
    # Calculate > cut cumulative sum
    histSum = get_hist_sum(histogram, range(0, nbins+2))
    cumulativeSumUpper = [ histSum - cumulativeSumLower[ibin] for ibin in range(nbins+2) ]

    return cumulativeSumLower, cumulativeSumUpper, histSum


def calc_efficiencies(cumulativeSum):
    """
    Calculate efficiencies from cumulative sums of the histogram.
 
    Parameters
    ----------
    cumulativeSum: list < float >

    Returns
    -------
    efficiencies: list < float >
    """

    efficiencies = {}
    sumHistogram = max(cumulativeSum)
    nbins = len(cumulativeSum) - 2 
 
    efficiencies = [ 100 * cumulativeSum[i]/sumHistogram for i in range(nbins+2) ]

    return efficiencies


def calc_significance(cumulativeSumSig, cumulativeSumBkg):
    """
    Calculate significance given cumulative sum of signal and background histograms.

    Parameters
    ----------
    cumulativeSumSig: list < float >
    cumulativeSumBkg: list < float >

    Returns
    -------
    significance: list < float >
    """

    nbins = len(cumulativeSumSig) - 2
    
    # If cumulative sums start with zeros
    if cumulativeSumSig[0] + cumulativeSumBkg[0] == 0:
        ibin0 = 1
        while cumulativeSumSig[ibin0] + cumulativeSumBkg[ibin0] == 0:
            ibin0 += 1
        part = ibin0 * [0.]
        range_ = range(ibin0, nbins+2)
        significance = part + [ cumulativeSumSig[ibin] / math.sqrt(cumulativeSumSig[ibin] + cumulativeSumBkg[ibin]) for ibin in range_ ]

    # If cumulative sums ends with zeros
    elif cumulativeSumSig[nbins+1] + cumulativeSumBkg[nbins+1] == 0:
        ibin0 = nbins
        while cumulativeSumSig[ibin0] + cumulativeSumBkg[ibin0] == 0:
            ibin0 -= 1
        part = (nbins-ibin0) * [0.]
        range_ = range(0, ibin0+1)
        significance = [ cumulativeSumSig[ibin] / math.sqrt(cumulativeSumSig[ibin] + cumulativeSumBkg[ibin]) for ibin in range_ ] + part
    
    else:
        significance = [ cumulativeSumSig[ibin] / math.sqrt(cumulativeSumSig[ibin] + cumulativeSumBkg[ibin]) for ibin in range(nbins+2) ]

    return significance
    

def calc_auc(x, y, percentage=False):
    """
    Return AUC given the x and y coordinates of the ROC curve.
    If percentage=True, then coordinates are cut efficiencies in percentage.
    x and y must be coordinates for increasing cut efficiencies, e.g. from (0, 0) to (1, 1)
    (from (0, 0) to (100, 100) if percentage=True).
    The presence of coordinates (0, 0) and (1, 1) in x and y is facultative.

    Parameters
    ----------
    x: list < float >
    y: list < float >
    percentage: boolean, optional

    Returns
    -------
    auc: float
    """

    if len(x) != len(y):
        print("ERROR: x and y must have the same size.")
        sys.exit()

    n = len(x)
    auc = x[0]*y[0]/2
    for i in range(n-1):
        auc += (x[i+1]-x[i])*(y[i+1]+y[i])/2
    if percentage:
        auc += (100-x[n-1])*(100+y[n-1])/2
        auc /= 100**2
    else:
        auc += (1-x[n-1])*(1+y[n-1])/2

    return auc


