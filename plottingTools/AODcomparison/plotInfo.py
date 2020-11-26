from collections import OrderedDict


def getPlotInfo():
    '''
    Define all plot information:
       * comparisons
       * binning info
       * variables info

    1. Variable information
    1.1. Key synthax of the dictionary plotInfo["varInfo"]:
    Possible key syntax:
       * <var>
           e.g. MET, Jet_mass_0 etc...
       * <var>_common
           defines common information for all variables strating with <var>
       * !<var>
           defines a composite variable calculated from previously
           defined variables
         
    1.2. Description of the variable information
    handleName: edm::Handle to be used in Event::getByLabel
                corresponds to Type in edmDumpEventContent <AODfile>
    labelName : string to be used in Event::getByLabel
                corresponds to Module in edmDumpEventContent <AODfile>
    accessor  : the accessor to get the variable from the EDM object
                 e.g. ".pt()"
    leaf      : leaf name in the nanoAOD, e.g. "MET_pt"
    leadingIdx: 0 for leading, 1 for subleading etc...
    name      : the name you want to see in the plot title and file name
    min       : minimum of the histogram if plotInfo["plotInfo"][minMaxBinning"] == "minMaxKeys"
    max       : maximum ........................................................................
    logscale  : True/False  (y axis)
    plot      : True/False  make plot of the variable
    expr      : If the variable is a composite variable, defines its expression using
                #['var'] to refer to the variable var
                 e.g. for deltaPhi: (#['Jet_phi_0']-#['MET_phi'])%(2*np.pi)

    '''


    ## Dictionary storing all information about plots
    plotInfo = OrderedDict()


    ## Set here which comparison you want
    #  Use names as defined in AOD_FORMATS in compare.py
    plotInfo["comparisons"] = {
        "MiniAODvsPFNanoAOD": 1,
        "MiniAODvsNanoAOD"  : 0,
        "NanoAODvsPFNanoAOD": 0
        }


    ## Binning info
    plotInfo["plotInfo"] = {
        "minMaxBinning": "minMaxKeys"
        # minMaxKeys         : from min to max keys of the plotted variable
        # zeroMaxPerHistogram: from zero to max of the histogram 
        # minMaxPerHistogram : from min to max of the histogram 
    }


    plotInfo["varInfo"] = OrderedDict()

    # Create a pointer to the subdictionary plotInfo["varInfo"]
    varInfo = plotInfo["varInfo"]


    ## Define some defaults
    varInfo["defaults"] = {
        "plot"      : True,
        "logscale"  : True,
        "min"       : 0,
    }

    ## MET
    varInfo["MET_common"] = {
        "handleName": "vector<pat::MET>",
        "labelName" : "slimmedMETs",
        "leadingIdx": 0
        }
    varInfo["MET"] = {
        "name"      : "MET",
        #"accessor"  : ".energy()",
        "accessor"  : ".pt()",
        "leaf"      : "MET_pt",
        "max"       : 1500,
    #    "plot"      : False
    }
    varInfo["MET_phi"] = {
        "name"      : "MET_phi",
        "accessor"  : ".phi()",
        "leaf"      : "MET_phi",
        "min"       : -3.15,
        "max"       : +3.15,
        "plot"      : False
    }

    # Jet
    varInfo["Jet_common"] = {
        "handleName": "vector<pat::Jet>",
        "labelName" : "slimmedJets"
    }

    # Jet mass
    varInfo["Jet_mass_common"] = {
        "accessor"  : ".mass()",
        "leaf"      : "Jet_mass"
    }
    varInfo["Jet_mass_0"] = {
        "name"      : "Leading jet mass",
        "leadingIdx": 0,
        "max"       : 300,
    #    "plot"      : False
    }
    varInfo["Jet_mass_1"] = {
        "name"      : "Subleading jet mass",
        "leadingIdx": 1,
        "max"       : 300,
    #    "plot"      : False
    }
    varInfo["Jet_mass_2"] = {
        "name"      : "Subsubleading jet mass",
        "leadingIdx": 2,
        "max"       : 150,
    #    "plot"      : False
    }

    # Jet phi
    varInfo["Jet_phi_common"] = {
        "accessor"  : ".phi()",
        "leaf"      : "Jet_phi",
        "min"       : -3.15,
        "max"       : +3.15,
        "plot"      : False
    }
    varInfo["Jet_phi_0"] = {
        "name"      : "Leading jet phi",
        "leadingIdx": 0
    }
    varInfo["Jet_phi_1"] = {
        "name"      : "Subleading jet phi",
        "leadingIdx": 1
    }

    # Jet pt
    varInfo["Jet_pt_common"] = {
        "accessor"  : ".pt()",
        "leaf"      : "Jet_pt"
    }
    varInfo["Jet_pt_0"] = {
        "name"      : "Leading jet pt",
        "leadingIdx": 0,
        "max"       : 2000,
    #    "plot"      : False
    }
    varInfo["Jet_pt_1"] = {
        "name"      : "Subleading jet pt",
        "leadingIdx": 1,
        "max"       : 2000,
    #    "plot"      : False
    }
    varInfo["Jet_pt_2"] = {
        "name"      : "subsubleading jet pt",
        "leadingIdx": 2,
        "max"       : 1000,
    #    "plot"      : False
    }

    # Combined variables
    varInfo["!deltaPhi_common"] = {
        "min"       : 0,
        "max"       : 6.3,
    }
    varInfo["!deltaPhi_METj0"] = {
        "name"      : "delta phi MET jet0",
        "expr"      : "(#['Jet_phi_0']-#['MET_phi'])%(2*np.pi)",
    }

    varInfo["!deltaPhi_METj1"] = {
        "name"      : "delta phi MET jet1",
        "expr"      : "(#['Jet_phi_1']-#['MET_phi'])%(2*np.pi)",
    }

    varInfo["!deltaPhi_minMETj"] = {
        "name"      : "min delta phi MET jets",
        "expr"      : "min( (#['Jet_phi_0']-#['MET_phi'])%(2*np.pi), (#['Jet_phi_1']-#['MET_phi'])%(2*np.pi) )",
    }

    return(plotInfo)
