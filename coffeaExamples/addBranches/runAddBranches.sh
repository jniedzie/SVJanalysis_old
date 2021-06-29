processor=BranchesMaker

python addBranches.py --input test_QCD.root --processor ${processor} --schema Base --pfnanoaodversion 106X --chunksize 5000 --debug
python addBranches.py --input test_QCD.root --processor ${processor} --schema PFNanoAOD --pfnanoaodversion 106X --chunksize 5000 --debug
python addBranches.py --input test_tchannel.root --processor ${processor} --schema Base --pfnanoaodversion 102X --chunksize 5000 --debug
python addBranches.py --input test_tchannel.root --processor ${processor} --schema PFNanoAOD --pfnanoaodversion 102X --chunksize 5000 --debug
