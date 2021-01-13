import sys
import json
from collections import OrderedDict
import time
import numpy as np

sys.path.append("../pythonUtils/")
import utilities as utl

import ROOT
ROOT.gROOT.SetBatch(True)


## Custom C++ function to be used in Define
ROOT.gInterpreter.Declare("""
Float_t DeltaPhiMinN(int NjetsMax, ROOT::VecOps::RVec<Float_t>& phi, Float_t& phi2) {
    ROOT::VecOps::RVec<Float_t> dPhi;
    int size = (int)phi.size();
    auto idxMax = std::min(size, NjetsMax);

    for (auto idx=0; idx<idxMax; idx++) {
        dPhi.push_back(ROOT::VecOps::DeltaPhi(phi[idx], phi2));
    }

    return(ROOT::VecOps::Min(dPhi));
}
""")


## Constants
HIST_OUTPUT_PATH = "/eos/user/f/fleble/SVJ/data/histograms/"
LUMI = 21071.0+38654.0

# Restrict histogram making to a couple of samples
# [] for all samples in sample.json file
# ["xxx", "yyy"] for xxx and yyy samples
LIST_OF_SAMPLES = ["tchannel_mMed-3000_mDark-20_rinv-0p3"]


## Useful functions

# Book a histogram for a specific variable
def bookHistogram(df, variable, range_, weight):
    return df.Histo1D(ROOT.ROOT.RDF.TH1DModel(variable, variable, range_[0], range_[1], range_[2]), variable, weight)

# Write a histogram with a given name to the output ROOT file
def writeHistogram(h, name):
    h.SetName(name)
    h.Write()


## Main function making histograms

def main(samples, variables, binning):

    # Set up multi-threading capability of ROOT
    ROOT.ROOT.EnableImplicitMT()

    # Loop through datasets and produce histograms of variables
    for sample in samples.keys():

        process = samples[sample]
        XSection = process["XSection"]

        # Create output file
        ROOTfileName = HIST_OUTPUT_PATH + sample + ".root"
        print("\nCreating ROOT file %s" %ROOTfileName)
        tfile = ROOT.TFile(ROOTfileName, "RECREATE")

        # Initialise event counter
        nGenEvts = 0

        # Make batches of files with a total of less than X million events
        # Need to do that because RDataFrame efficient features seems to break down
        # for too many events at once
        batches = [[]]
        nEvts = 0
        for file_ in process["fileset"]:

            # If file is on eos, add global redirector
            if file_.startswith("/store/mc/") or file_.startswith("/store/user/"):
                file_ = "root://cms-xrd-global.cern.ch/" + file_

            if nEvts < 2e6:
                batches[-1].append(file_)
                f = ROOT.TFile.Open(file_, "READ")
                t = f.Get("Events")
                nEvts += t.GetEntries()
            else:
                nEvts = 0
                batches.append([file_])

        print("%d batches of files made" %(len(batches)))

        # Object to store histograms from the different batches of file sets
        hists = []


        # Loop over all batches of files
        for fileset in batches:
            print("")

            # Object to store partial histograms from the different files in this batch
            histsBatch = []

            # Loop over all files
            nFiles = len(fileset)
            for iFile, fileName in enumerate(fileset):

                tstart = time.time()
                print("Reading file %d out of %d" %(iFile+1, nFiles))

                histsBatch.append({})

                # If file is on eos, add global redirector
                if fileName.startswith("/store/mc/") or fileName.startswith("/store/user/"):
                    fileName = "root://cms-xrd-global.cern.ch/" + fileName

                # Lists to store pointers to different RDataFrames (different filters)
                # and to Ordered dictionary storing the definition of new quantities
                dfList = []
                objDefinitionsList = []

                # Increment the number of event processed
                dfRun = ROOT.ROOT.RDataFrame("Runs", fileName)
                if "genEventSumw" in dfRun.GetColumnNames():
                    nGenEvts += dfRun.Sum("genEventSumw").GetValue()
                else:
                    nGenEvts += dfRun.Sum("genEventSumw_").GetValue()

                # In case one wants to recompute the sum of gen weights, do the following
                # instead of the previous block of instructions
                #nGenEvts += df.Sum("genWeight").GetValue()  


                ## Now we define the variables needed to fill in the histograms

                ## RDF without any cut
                dfIdx = 0

                # Load dataset
                dfList.append(ROOT.ROOT.RDataFrame("Events", fileName))

                # Define objects and variables
                objDefinitionsList.append(OrderedDict())
                o = objDefinitionsList[dfIdx]
                o["GoodFatJet"           ] = "abs(FatJet_eta) < 2.4 && FatJet_pt > 200"
                o["nGoodFatJet"          ] = "Sum(GoodFatJet)"
                o["GoodFatJet_pt"        ] = "FatJet_pt[GoodFatJet]"
                o["GoodFatJet_phi"       ] = "FatJet_phi[GoodFatJet]"
                o["GoodFatJet_eta"       ] = "FatJet_eta[GoodFatJet]"
                o["GoodFatJet_mass"      ] = "FatJet_mass[GoodFatJet]"
                o["GoodFatJet_msoftdrop" ] = "FatJet_msoftdrop[GoodFatJet]"
                o["GoodFatJet_tau1"      ] = "FatJet_tau1[GoodFatJet]"
                o["GoodFatJet_tau2"      ] = "FatJet_tau2[GoodFatJet]"
                o["GoodFatJet_tau3"      ] = "FatJet_tau3[GoodFatJet]"
                o["GoodFatJet_tau21"     ] = "GoodFatJet_tau2/GoodFatJet_tau1"
                o["GoodFatJet_tau32"     ] = "GoodFatJet_tau3/GoodFatJet_tau2"
                o["HT_AK8"               ] = "Sum(GoodFatJet_pt)"
                o["ST_AK8"               ] = "HT_AK8 + MET_pt"
                o["METrHT_AK8"           ] = "MET_pt / HT_AK8"
                o["METrST_AK8"           ] = "MET_pt / ST_AK8"

                o["GoodJet"              ] = "abs(Jet_eta) < 2.4 && Jet_pt > 30"
                o["nGoodJet"             ] = "Sum(GoodJet)"
                o["GoodJet_pt"           ] = "Jet_pt[GoodJet]"
                o["GoodJet_phi"          ] = "Jet_phi[GoodJet]"
                o["GoodJet_eta"          ] = "Jet_eta[GoodJet]"
                o["GoodJet_mass"         ] = "Jet_mass[GoodJet]"
                o["HT_AK4"               ] = "Sum(GoodJet_pt)"
                o["ST_AK4"               ] = "HT_AK4 + MET_pt"
                o["METrHT_AK4"           ] = "MET_pt / HT_AK4"
                o["METrST_AK4"           ] = "MET_pt / ST_AK4"


                for obj, definition in objDefinitionsList[dfIdx].items():
                    dfList[dfIdx] = dfList[dfIdx].Define(obj, definition)


                ## RDF with at least 1 good FatJet
                dfIdx = 1
                dfList.append(dfList[0].Filter("nGoodFatJet >= 1", "Greater than 1 good FatJets"))
                 
                # Define objects and variables
                objDefinitionsList.append(OrderedDict())
                o = objDefinitionsList[dfIdx]
                o["J1_pt"                ] = "GoodFatJet_pt[0]"
                o["J1_eta"               ] = "GoodFatJet_eta[0]"
                o["J1_phi"               ] = "GoodFatJet_phi[0]"
                o["J1_mass"              ] = "GoodFatJet_mass[0]"
                o["J1_msoftdrop"         ] = "GoodFatJet_msoftdrop[0]"
                o["J1_tau1"              ] = "GoodFatJet_tau1[0]"
                o["J1_tau2"              ] = "GoodFatJet_tau2[0]"
                o["J1_tau3"              ] = "GoodFatJet_tau3[0]"
                o["J1_tau21"             ] = "GoodFatJet_tau21[0]"
                o["J1_tau32"             ] = "GoodFatJet_tau32[0]"
                o["dPhi_J1MET"           ] = "abs(ROOT::VecOps::DeltaPhi(GoodFatJet_phi[0], MET_phi))"
                o["dPhiMin_JMET"         ] = "ROOT::VecOps::Min(abs(ROOT::VecOps::DeltaPhi(GoodFatJet_phi, MET_phi)))"
                o["dPhiMinUpTo2_JMET"    ] = "DeltaPhiMinN(2, GoodFatJet_phi, MET_phi)"
                o["dPhiMinUpTo4_JMET"    ] = "DeltaPhiMinN(4, GoodFatJet_phi, MET_phi)"

                for obj, definition in objDefinitionsList[dfIdx].items():
                    dfList[dfIdx] = dfList[dfIdx].Define(obj, definition)


                ## RDF with at least 1 good AK4 Jet
                dfIdx = 2
                dfList.append(dfList[0].Filter("nGoodJet >= 1", "Greater than 1 good Jets"))

                # Define objects and variables
                objDefinitionsList.append(OrderedDict())
                o = objDefinitionsList[dfIdx]
                o["j1_pt"                ] = "GoodJet_pt[0]"
                o["j1_eta"               ] = "GoodJet_eta[0]"
                o["j1_phi"               ] = "GoodJet_phi[0]"
                o["j1_mass"              ] = "GoodJet_mass[0]"
                o["dPhi_j1MET"           ] = "abs(ROOT::VecOps::DeltaPhi(GoodJet_phi[0], MET_phi))"
                o["dPhiMin_jMET"         ] = "ROOT::VecOps::Min(abs(ROOT::VecOps::DeltaPhi(GoodJet_phi, MET_phi)))"
                o["dPhiMinUpTo2_jMET"    ] = "DeltaPhiMinN(2, GoodJet_phi, MET_phi)"
                o["dPhiMinUpTo4_jMET"    ] = "DeltaPhiMinN(4, GoodJet_phi, MET_phi)"

                for obj, definition in objDefinitionsList[dfIdx].items():
                    dfList[dfIdx] = dfList[dfIdx].Define(obj, definition)


                ## RDF with at least 2 good FatJets
                dfIdx = 3

                # Define dataframe
                dfList.append(dfList[1].Filter("nGoodFatJet >= 2", "Greater than 2 good FatJets"))
                 
                # Define objects and variables
                objDefinitionsList.append(OrderedDict())
                o = objDefinitionsList[dfIdx]
                o["J2_pt"             ] = "GoodFatJet_pt[1]"
                o["J2_eta"            ] = "GoodFatJet_eta[1]"
                o["J2_phi"            ] = "GoodFatJet_phi[1]"
                o["J2_mass"           ] = "GoodFatJet_mass[1]"
                o["J2_msoftdrop"      ] = "GoodFatJet_msoftdrop[1]"
                o["J2_tau1"           ] = "GoodFatJet_tau1[1]"
                o["J2_tau2"           ] = "GoodFatJet_tau2[1]"
                o["J2_tau3"           ] = "GoodFatJet_tau3[1]"
                o["J2_tau21"          ] = "GoodFatJet_tau21[1]"
                o["J2_tau32"          ] = "GoodFatJet_tau32[1]"

                o["dR_J1J2"           ] = "ROOT::VecOps::DeltaR(J1_eta, J2_eta, J1_phi, J2_phi)"
                o["dPhiMin2_JMET"     ] = "DeltaPhiMinN(2, GoodFatJet_phi, MET_phi)"
                o["dEta_J1J2"         ] = "std::abs(GoodFatJet_eta[0] - GoodFatJet_eta[1])"
                o["dPhi_J2MET"        ] = "std::abs(ROOT::VecOps::DeltaPhi(GoodFatJet_phi[1], MET_phi))"
                o["PtEtaPhiM_J1"      ] = "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > (GoodFatJet_pt[0], GoodFatJet_eta[0], GoodFatJet_phi[0], GoodFatJet_mass[0])"
                o["PtEtaPhiM_J2"      ] = "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > (GoodFatJet_pt[1], GoodFatJet_eta[1], GoodFatJet_phi[1], GoodFatJet_mass[1])"
                o["PtEtaPhiM_J1J2"    ] = "PtEtaPhiM_J1 + PtEtaPhiM_J2"

                o["J1J2_mass"         ] = "PtEtaPhiM_J1J2.M()"
                o["J1J2_mass2"        ] = "std::pow(J1J2_mass, 2)"
                o["J1J2_pt"           ] = "PtEtaPhiM_J1J2.Pt()"
                o["J1J2_pt2"          ] = "std::pow(J1J2_pt, 2)"
                o["J1J2_phi"          ] = "PtEtaPhiM_J1J2.Phi()"
                o["dPhi_J1J2MET"      ] = "std::abs(ROOT::VecOps::DeltaPhi(J1J2_phi, MET_phi))"

                o["MT_AK8"            ] = "std::sqrt( J1J2_mass2  +  2 * ( std::sqrt(J1J2_mass2 + J1J2_pt2) * MET_pt - MET_pt * J1J2_pt * std::cos(dPhi_J1J2MET) ) )"
                o["RT_AK8"            ] = "MET_pt / MT_AK8"
             
                for obj, definition in objDefinitionsList[dfIdx].items():
                    dfList[dfIdx] = dfList[dfIdx].Define(obj, definition)


                ## RDF with at least 2 good Jets
                dfIdx = 4

                # Define dataframe
                dfList.append(dfList[2].Filter("nGoodJet >= 2", "Greater than 2 good Jets"))

                # Define objects and variables
                objDefinitionsList.append(OrderedDict())
                o = objDefinitionsList[dfIdx]
                o["j2_pt"             ] = "GoodJet_pt[1]"
                o["j2_eta"            ] = "GoodJet_eta[1]"
                o["j2_phi"            ] = "GoodJet_phi[1]"
                o["j2_mass"           ] = "GoodJet_mass[1]"

                o["dR_j1j2"           ] = "ROOT::VecOps::DeltaR(j1_eta, j2_eta, j1_phi, j2_phi)"
                o["dPhiMin2_jMET"     ] = "DeltaPhiMinN(2, GoodJet_phi, MET_phi)"
                o["dEta_j1j2"         ] = "std::abs(GoodJet_eta[0] - GoodJet_eta[1])"
                o["dPhi_j2MET"        ] = "std::abs(ROOT::VecOps::DeltaPhi(GoodJet_phi[1], MET_phi))"
                o["PtEtaPhiM_j1"      ] = "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > (GoodJet_pt[0], GoodJet_eta[0], GoodJet_phi[0], GoodJet_mass[0])"
                o["PtEtaPhiM_j2"      ] = "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > (GoodJet_pt[1], GoodJet_eta[1], GoodJet_phi[1], GoodJet_mass[1])"
                o["PtEtaPhiM_j1j2"    ] = "PtEtaPhiM_j1 + PtEtaPhiM_j2"

                o["j1j2_mass"         ] = "PtEtaPhiM_j1j2.M()"
                o["j1j2_mass2"        ] = "std::pow(j1j2_mass, 2)"
                o["j1j2_pt"           ] = "PtEtaPhiM_j1j2.Pt()"
                o["j1j2_pt2"          ] = "std::pow(j1j2_pt, 2)"
                o["j1j2_phi"          ] = "PtEtaPhiM_j1j2.Phi()"
                o["dPhi_j1j2MET"      ] = "std::abs(ROOT::VecOps::DeltaPhi(j1j2_phi, MET_phi))"

                o["MT_AK4"            ] = "std::sqrt( j1j2_mass2  +  2 * ( std::sqrt(j1j2_mass2 + j1j2_pt2) * MET_pt - MET_pt * j1j2_pt * std::cos(dPhi_j1j2MET) ) )"
                o["RT_AK4"            ] = "MET_pt / MT_AK4"
                
                for obj, definition in objDefinitionsList[dfIdx].items():
                    dfList[dfIdx] = dfList[dfIdx].Define(obj, definition)


                ## Book histograms
                for variable in variables:

                    # Check if binning is defined for this variable
                    if variable not in binning["noregex"]:
                        regexes = list(binning["regex"].keys())
                        indices = utl.inregex(variable, regexes)
                        if len(indices) == 0:
                            print("Binning of %s is not defined. Skipping." %variable)
                            continue
                        elif len(indices) > 1:
                            print("%s matches several regexes. Binning cannot be defined. Skipping." %variable)
                            continue
                        else: 
                            binning_ = binning["regex"][regexes[indices[0]]]
                    else:
                        binning_ = binning["noregex"][variable]

                    # Find from which RDataFrame the variable has been defined
                    for idx in range(len(objDefinitionsList)):
                        if variable in objDefinitionsList[idx].keys():
                            break
                        else:
                            # If variable is not defined, then it's taken from the uncut dataframe
                            if idx == len(objDefinitionsList)-1:
                                idx = 0
                    
                    # Book histogram
                    # Histograms should NOT be added together yet as RDataFrame does not proceed in one loop
                    # Should be done at the very end
                    weight = "genWeight"
                    histsBatch[iFile][variable] = bookHistogram(dfList[idx], variable, binning_, weight)

                elapsed = time.time() - tstart
                print("Elapsed time: %d s" %elapsed)


            ## Adding histograms of the batch together
            print("Adding together histograms of batch %d..." %(len(hists)+1))
            tstrat = time.time()
            hists.append({})
            for variable in histsBatch[0].keys():
                hists[-1][variable] = histsBatch[0][variable]
                for iFile in range(1, len(histsBatch)):
                    hists[-1][variable].Add(histsBatch[iFile][variable].GetValue())
            elapsed = time.time() - tstart
            print("Elapsed time: %d s" %elapsed)

 
        ## Normalize histograms and write to output ROOT file
        print("\nNormalizing histograms and writing to output ROOT file...")
        tstart = time.time()
        tfile.cd()
        for variable in hists[0].keys():
            hist = hists[0][variable]
            for iFile in range(1, len(hists)):
                hist.Add(hists[iFile][variable].GetValue())
            hist.Scale( XSection*LUMI/nGenEvts )
            writeHistogram(hist, "{}".format(variable))
            print("%s histogram saved" %variable)
        elapsed = time.time() - tstart
        print("Elapsed time: %d s" %elapsed)

        tfile.Close()



if __name__ == "__main__":

    tstart = time.time()

    ## Define samples for which to make histograms
    #  key is the sample name
    #  fileset denotes the path to the files composing this sample dataset
    #  XSection is the cross-section of the sample process
    with open("samples.json", 'r') as f:
        samples0 = json.load(f)

    # Restrict to list of sample
    if LIST_OF_SAMPLES == []:
        samples = samples0
    else:
        samples = { sample: samples0[sample] for sample in samples0.keys() if sample in LIST_OF_SAMPLES }


    ## Variables to histogram
    with open("variables.json", 'r') as f:
        variables = json.load(f)["variables"]


    ## Define the binning of the different variables to histogram
    with open("binning.json", 'r') as f:
        binning = json.load(f)


    ## Make histograms
    main(samples, variables, binning)


    elapsed = time.time() - tstart
    print("\nTotal elapsed time: %d s" %elapsed)
