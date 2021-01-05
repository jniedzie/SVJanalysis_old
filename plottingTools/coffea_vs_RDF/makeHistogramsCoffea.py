import coffea as cf
import uproot
import awkward1 as ak

import time


LUMI = 21071.0+38654.0
XSection = 0.01891


# t-channel files
samplesPath = "/eos/home-f/fleble/SVJ/data/production/102X/tchannel/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR/"
fileset = {
    't-channel': [ samplesPath+str(i)+".root" for i in range(1, 101) ]
}


class histogramProcessor(cf.processor.ProcessorABC):
    def __init__(self):
        pt_axis = cf.hist.Bin("pt", r"$p_{T}$ [GeV]",   200,    0.0,    2000.0)

        self._accumulator = cf.processor.dict_accumulator({
            'h_FatJet_pt': cf.hist.Hist("h_FatJet_pt", pt_axis),
            'cutflow': cf.processor.defaultdict_accumulator(float)
            })


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        output = self.accumulator.identity()

        # Define variables to histogram
        FatJet_pt = ak.flatten(events["FatJet_pt"])
        FatJet_eta = ak.flatten(events["FatJet_eta"])
        goodFatJet = ( (FatJet_pt>200) & (abs(FatJet_eta)<2.4) )
        goodFatJet_pt = FatJet_pt[goodFatJet]

        # Recompute sum of genWeights (not to have to open the Runs tree)
        sumGenWeight = sum(events["genWeight"])

        # Reshape genWeight array to match FatJet awkward-array
        genWeight = ak.flatten(ak.broadcast_arrays(events["genWeight"], events["FatJet_pt"])[0])[goodFatJet]

        # Fill histogram
        output['h_FatJet_pt'].fill(pt=goodFatJet_pt, weight=genWeight)

        # Increment number of generated event
        output["cutflow"]["all"] += sumGenWeight

        return output


    def postprocess(self, accumulator):
        return accumulator


def main():

    ## Run the processor
    output = cf.processor.run_uproot_job(
        fileset,
        treename='Events',
        processor_instance=histogramProcessor(),
        executor=cf.processor.futures_executor,
        executor_args={"workers": 8},
        chunksize=100000,
        maxchunks=None,
    )


    # Export the histograms to root files
    fout = uproot.recreate("hist.root")
    nGenEvts = output["cutflow"]["all"]

    for key, hist in output.items():
        if isinstance(hist, cf.hist.Hist):
            hist.scale(LUMI * XSection / nGenEvts)
            fout[key] = cf.hist.export1d(hist)
    fout.close()


if __name__ == "__main__":
    tstart = time.time()
    main()
    elapsed = time.time() - tstart
    print("Elapsed time: %d s" %elapsed)
