#!/usr/bin/env python3
"""This is an example of configuring and using the module for deriving the jetId variable."""

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.jetId import jetId

json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2024_Winter24/jetid.json.gz"

jetId = jetId(json, jetType="AK4PUPPI")

# Settings for post-processor
from argparse import ArgumentParser
parser = ArgumentParser(description="Process nanoAOD files and add branches",epilog="Good luck!")
parser.add_argument('-i', '--infiles', nargs='+')
parser.add_argument('-m', '--maxevts', type=int, default=10000) # limit number of events (per file) for testing
args = parser.parse_args()
branchsel = None
fnames = args.infiles or [
      "root://cms-xrd-global.cern.ch//store/data/Run2024C/MuonEG/NANOAOD/MINIv6NANOv15-v1/2530000/694ee48f-c22e-4ae3-83f8-19f98edb481b.root",
]

# Process nanoAOD file
p = PostProcessor(".", fnames, "nJet>0 && Jet_pt>30",
                  branchsel,
                  [jetId],
                  maxEntries=args.maxevts,
                  postfix="_jetId",
                  outputbranchsel=branchsel,
                  provenance=True, prefetch=True, longTermCache=True)
p.run()
