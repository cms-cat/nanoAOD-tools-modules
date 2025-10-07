#!/usr/bin/env python3
"""This is an example of configuring and using the module for jet veto maps."""

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.jetVetoMap import jetVMAP

json_JVMAP = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/2025-09-23/jetvetomaps.json.gz"
corrName = "Summer22_23Sep2023_RunCD_V1"
veto_map_name = "jetvetomap"

jetVetoMaps = jetVMAP(json_JVMAP, corrName, veto_map_name)

# Settings for post-processor
from argparse import ArgumentParser
parser = ArgumentParser(description="Process nanoAOD files and add branches",epilog="Good luck!")
parser.add_argument('-i', '--infiles', nargs='+')
parser.add_argument('-m', '--maxevts', type=int, default=10000) # limit number of events (per file) for testing
args = parser.parse_args()
branchsel = None
fnames = args.infiles or [
  "root://cms-xrd-global.cern.ch//store/mc/Run3Summer22NanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/30000/542e888a-24b0-4e51-a494-d2690f894c0d.root",
]

# Process nanoAOD file
p = PostProcessor(".", fnames, None,
                  branchsel,
                  [jetVetoMaps],
                  maxEntries=args.maxevts,
                  postfix="_jetVeto",
                  outputbranchsel=branchsel,
                  provenance=True, prefetch=True, longTermCache=True)
p.run()
