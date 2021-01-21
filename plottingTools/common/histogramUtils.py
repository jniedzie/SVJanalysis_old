
def getHistogram(tfile, histName):
    h = tfile.Get(histName)
    if not h:
        raise Exception("Failed to load histogram {}.".format(histName))
    return h


def getHistSum(hist, range_):
    """Return sum of all bins of a TH1 histogram."""
    histSum = 0
    for ibin in range_:
        histSum += hist.GetBinContent(ibin)
    return histSum


def calcCumulativeSums(histogram):
    # Calculate number of bins (GetNbins returns number of bins except underflow and overflow bins)
    xaxis = histogram.GetXaxis()
    nbins = xaxis.GetNbins()
    # Calculate < cut cumulative sum
    cumulativeSumLower = [ histogram.GetBinContent(0) ]
    for ibin in range(1, nbins+1):
        cumulativeSumLower.append(cumulativeSumLower[-1] + histogram.GetBinContent(ibin))
    # Calculate > cut cumulative sum
    histSum = getHistSum(histogram, range(0, nbins+1))
    cumulativeSumUpper = [ histSum - cumulativeSumLower[ibin] for ibin in range(nbins+1) ]

    return cumulativeSumLower, cumulativeSumUpper, histSum
