#!/usr/bin/env python3
"""This is an example of configuring and using the module to derive the jetId variable."""

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.jetId import jetId

testVersion = 15

if testVersion == 12 :
    jetId = jetId(None, nanoVersion=12)
    defaultFile = "root://cms-xrd-global.cern.ch//store/mc/Run3Summer22NanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/30000/542e888a-24b0-4e51-a494-d2690f894c0d.root"
else:
    json = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-07-17/jetid.json.gz"
    jetId = jetId(json, nanoVersion=15)
    defaultFile = "root://cms-xrd-global.cern.ch//store/data/Run2024C/MuonEG/NANOAOD/MINIv6NANOv15-v1/2530000/694ee48f-c22e-4ae3-83f8-19f98edb481b.root"

# Settings for post-processor
from argparse import ArgumentParser
parser = ArgumentParser(description="Process nanoAOD files and add branches",epilog="Good luck!")
parser.add_argument('-i', '--infiles', nargs='+')
parser.add_argument('-m', '--maxevts', type=int, default=10000) # limit number of events (per file) for testing
args = parser.parse_args()
branchsel = None
fnames = args.infiles or [defaultFile]

# Process nanoAOD file
p = PostProcessor(".", fnames, "nJet>0 && Jet_pt>30",
                  branchsel,
                  [jetId],
                  maxEntries=args.maxevts,
                  postfix="_jetId",
                  outputbranchsel=branchsel,
                  provenance=True, prefetch=True, longTermCache=True)
p.run()
