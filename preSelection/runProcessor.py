from coffea import processor
from coffea.nanoevents import BaseSchema
import awkward as ak
import awkward0 as ak0
import uproot3
import numpy as np

from preSelectionProcessor import Preselection


## Fileset
fname = "/afs/cern.ch/user/f/fleble/svjets/t_channel/samples/102X/mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/NANOAODJMAR/1.root"
fileset = {
    "baseline": [fname]
}


## Make pre-selections
output = processor.run_uproot_job(
    fileset,
    treename = "Events",
    processor_instance = Preselection(),
    executor = processor.iterative_executor,
    executor_args = {"schema": BaseSchema, "workers": 4},
)


## Print out cutflow
cutflow = output.pop("cutflow")
lenCol1 = max([ len(k) for k in cutflow.keys() ])

print("\nCutflow:")
print("\tCut" + (lenCol1-3)*" " + "  Abs. eff. [%]   Rel. eff. [%]")
nAll = cutflow["all"].value
for cut, n in cutflow.items():
    if cut != "all":
        print("\t%s%s  %.2f           %.2f" %(cut, (lenCol1-len(cut))*" ", 100*n.value/nAll, 100*n.value/nPreviousCut))
    nPreviousCut = n.value


## Making output ROOT file
branches = {}
branches_init = {}
lenKeys = []

# Finding keys giving the length of jagged arrays
for k, v in output.items():
    nKey = "n"+str(k.split("_")[0])
    if (not k.startswith("n")) and (nKey in output.keys()):
        lenKeys.append(nKey)

# Making branches
# Need to use ak0 because it is not yet implemented in uproot4 i.e. ak1 (Feb. 2021)
for k, v in output.items():
    if not k.startswith("n"):
        nKey = "n"+str(k.split("_")[0])
        if nKey in lenKeys:
            branches[k] = ak.to_awkward0(ak.Array(v.value))
            branches_init[k] = uproot3.newbranch(np.dtype("f8"), size=nKey)
        else:
            branches[k] = ak.to_awkward0(ak.Array(v.value))
            branches_init[k] = uproot3.newbranch(v.value.dtype)
    else:
        if k in lenKeys:
            branches[k] = ak.to_awkward0(ak.Array(v.value))
        else:
            branches[k] = ak.to_awkward0(ak.Array(v.value))
            branches_init[k] = uproot3.newbranch(v.value.dtype)

# Save branches to ROOT file
# Need to use uproot3 because it is not implemented yet in uproot4 (Feb. 2021)
with uproot3.recreate("test.root") as f:
    f["Events"] = uproot3.newtree(branches_init)
    f["Events"].extend(branches)

