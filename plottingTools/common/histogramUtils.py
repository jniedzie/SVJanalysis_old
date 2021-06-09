import math


def get_histogram(tfile, hist_name):
    """Return histogram from ROOT file.

    Args:
        tfile (ROOT.TFile)
        hist_name (str)

    Returns:
        ROOT.TH1
    """

    h = tfile.Get(hist_name)
    if not h:
        raise Exception("Failed to load histogram {}.".format(hist_name))
    return h


def normalize(histogram, norm=1.):
    """Normalize a TH1 histogram to 1."""

    integral = histogram.Integral()
    histogram.Scale(norm / integral)

    return
    

def rebin_histogram(histogram, nbins_max, verbose=False):
    """Rebin a TH1 histogram such that it has a less than nBinsMax bins.

    Args:
        histogram (ROOT.TH1)
        nbins_max (int): Maximum number of bins
        verbose (bool)

    Returns:
        int: New number of bins

    """

    xaxis = histogram.GetXaxis()
    nbins = xaxis.GetNbins()
    
    rebin_found = False
    ngroup = 2
    while not rebin_found:
        if nbins%ngroup == 0 and nbins/ngroup < nbins_max:
            rebin_found = True
        else:
            ngroup += 1

    nbins = int(nbins/ngroup)
    if verbose:
        print("Rebin histogram with %d groups. New numbers of bins = %d" %(ngroup, nbins))
    histogram.Rebin(ngroup)

    return nbins


def add_overflow_bin_to_last_bin(histogram):
    """Add content of the overflow bin to the last bin.
    
    Args:
        histogram (ROOT.TH1)

    Returns:
        None
    """

    xaxis = histogram.GetXaxis()
    nbins = xaxis.GetNbins()

    overflow = histogram.GetBinContent(nbins+1)
    last_bin_content = histogram.GetBinContent(nbins)
    histogram.SetBinContent(nbins+1, 0.)
    histogram.SetBinContent(nbins, last_bin_content+overflow)

    return


def get_hist_sum(hist, range_=None):
    """Return sum of all or some bins of a TH1 histogram.

    Args:
        hist (ROOT.TH1)
        range_ (range):
            If None, then sum of all bins, including underflow and overflow bins,
            is returned

    Returns:
        float
    """

    hist_sum = 0

    if range_ is None:
        xaxis = histogram.GetXaxis()
        nbins = xaxis.GetNbins()
        range_ = range(0, nbins+2)

    for ibin in range_:
        hist_sum += hist.GetBinContent(ibin)

    return hist_sum


def calc_cumulative_sums(histogram):
    """Calculate cumulative sum of a TH1 histogram.
    
    Args:
        histogram (ROOT.TH1D)
    
    Returns:
        tuple:
            cumulativeSumLower (list[float])
            cumulativeSumUpper (list[float])
            hist_sum (float)
    """

    # Calculate number of bins (GetNbins returns number of bins except underflow and overflow bins)
    xaxis = histogram.GetXaxis()
    nbins = xaxis.GetNbins()
    # Calculate < cut cumulative sum
    cumulative_sum_lower = [ histogram.GetBinContent(0) ]
    for ibin in range(1, nbins+2):
        cumulative_sum_lower.append(cumulative_sum_lower[-1] + histogram.GetBinContent(ibin))
    # Calculate > cut cumulative sum
    hist_sum = get_hist_sum(histogram, range(0, nbins+2))
    cumulative_sum_upper = [ hist_sum - cumulative_sum_lower[ibin] for ibin in range(nbins+2) ]

    return cumulative_sum_lower, cumulative_sum_upper, hist_sum


def calc_efficiencies(cumulative_sum):
    """Calculate efficiencies from cumulative sums of the histogram.
 
    Args:
        cumulative_sum (list[float])

    Returns:
        list[float]
    """

    efficiencies = {}
    sum_histogram = max(cumulative_sum)
    nbins = len(cumulative_sum) - 2 
 
    efficiencies = [ 100 * cumulative_sum[i]/sum_histogram for i in range(nbins+2) ]

    return efficiencies


def calc_significance(cumulative_sum_sig, cumulative_sum_bkg):
    """Calculate significance given cumulative sum of signal and background histograms.

    Args:
        cumulative_sum_sig (list[float])
        cumulative_sum_bkg (list[float])

    Returns:
        list[float]
    """

    nbins = len(cumulative_sum_sig) - 2
    
    # If cumulative sums start with zeros
    if cumulative_sum_sig[0] + cumulative_sum_bkg[0] == 0:
        ibin0 = 1
        while cumulative_sum_sig[ibin0] + cumulative_sum_bkg[ibin0] == 0:
            ibin0 += 1
        part = ibin0 * [0.]
        range_ = range(ibin0, nbins+2)
        significance = part + [ cumulative_sum_sig[ibin] / math.sqrt(cumulative_sum_sig[ibin] + cumulative_sum_bkg[ibin]) for ibin in range_ ]

    # If cumulative sums ends with zeros
    elif cumulative_sum_sig[nbins+1] + cumulative_sum_bkg[nbins+1] == 0:
        ibin0 = nbins
        while cumulative_sum_sig[ibin0] + cumulative_sum_bkg[ibin0] == 0:
            ibin0 -= 1
        part = (nbins-ibin0) * [0.]
        range_ = range(0, ibin0+1)
        significance = [ cumulative_sum_sig[ibin] / math.sqrt(cumulative_sum_sig[ibin] + cumulative_sum_bkg[ibin]) for ibin in range_ ] + part
    
    else:
        significance = [ cumulative_sum_sig[ibin] / math.sqrt(cumulative_sum_sig[ibin] + cumulative_sum_bkg[ibin]) for ibin in range(nbins+2) ]

    return significance
    

def calc_auc(x, y, percentage=False):
    """Return AUC given the x and y coordinates of the ROC curve.

    If percentage=True, then coordinates are cut efficiencies in percentage.
    x and y must be coordinates for increasing cut efficiencies, e.g. from (0, 0) to (1, 1)
    (from (0, 0) to (100, 100) if percentage=True).
    The presence of coordinates (0, 0) and (1, 1) in x and y is facultative.

    Args:
        x (list[float])
        y (list[float])
        percentage (bool, optional)

    Returns:
        float
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


