from coffea import processor
from coffea.nanoevents import BaseSchema
import awkward as ak
import numpy as np
import pandas as pd


## Coffea processor
class Processor(processor.ProcessorABC):
    """Computes triggers efficiency"""

    def __init__(self, triggers):
        self.triggers = triggers
        self._accumulator = processor.dict_accumulator({})
        for trigger in self.triggers:
            self._accumulator[trigger] = processor.value_accumulator(int, initial=0)
        self._accumulator["nEvt"] = processor.value_accumulator(int, initial=0)


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        output = self.accumulator.identity()

        for trigger in self.triggers:
            output[trigger] += ak.sum(getattr(events, trigger))
        output["nEvt"] += len(events.nJet)

        return output

    def postprocess(self, accumulator):
        return accumulator


## Opening file with all good trigger paths
dfTrigger = pd.read_csv("goodTriggerPath_2018.csv", delimiter=",")
triggers = dfTrigger["triggerPath"].to_list()
datasets = dfTrigger["dataset"].to_list()
triggers = [ x[:-2] if x.endswith("_v") else x for x in triggers ]


## Running coffea processor
#fname = "/afs/cern.ch/user/f/fleble/svjets/t_channel/samples/102X/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR/1.root"
fname = "/afs/cern.ch/user/f/fleble/svjets/t_channel/samples/102X/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR/merged.root"
fileset = {
    "baseline": [fname]
}

output = processor.run_uproot_job(
    fileset,
    treename = "Events",
    processor_instance = Processor(triggers),
    executor = processor.iterative_executor,
    executor_args = {"schema": BaseSchema, "workers": 4},
)


## Print out all trigger efficiencies
maxLen = max([len(trigger) for trigger in triggers])
col1Name = "Trigger path"
nSpaces = maxLen - len(col1Name)

print("\nProcessed %d events" %(output["nEvt"].value))
print(col1Name + nSpaces*" " + " | Efficiency [%]")

maxEff = 0.
maxEffTrigger = "--"
triggerListEfficiency = []
triggerList = []
datasetList = []

for trigger, dataset in zip(triggers, datasets):
    if trigger not in ("HLT_HcalNZS", "HLT_HcalCalibration", "HLT_HcalPhiSym"):
        eff = 100*output[trigger].value/output["nEvt"].value
        triggerList.append(trigger)
        triggerListEfficiency.append(eff)
        datasetList.append(dataset)

        nSpaces = maxLen - len(trigger)
        print(trigger + nSpaces*" " + " | %.1f" %eff)
        if eff > maxEff:
            maxEff = eff
            maxEffTrigger = trigger


## Print out max efficiency
print("")
print("Maximum efficiency: %.3f" %maxEff)
print("for trigger path %s" %maxEffTrigger)


## Save trigger effienciency to csv
df = pd.DataFrame({"trigger": triggerList, "dataset": datasetList, "efficiency [%]": triggerListEfficiency})
df.sort_values("efficiency [%]", inplace=True, ascending=False)
df.to_csv("triggerEfficiencies.csv", index=False)
