#!/usr/bin/env python3
"""This is an example of configuring and using the muon scale/smearing module."""

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.muonScaleRes import *

# Set up the correction module.
# json files are not included in the central cvmfs area at the time of release.
# see https://gitlab.cern.ch/cms-muonPOG/muonscarekit/-/tree/master/corrections for details.

json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22EE/muon_scalesmearing.json.gz"
muSS = muonScaleRes(json, overwritePt=True, is_mc=True, minPt=3.)


# Settings for post-processor
from argparse import ArgumentParser
parser = ArgumentParser(description="Process nanoAOD files and add branches",epilog="Good luck!")
parser.add_argument('-i', '--infiles', nargs='+')
parser.add_argument('-m', '--maxevts', type=int, default=10000) # limit number of events (per file) for testing
args = parser.parse_args()
branchsel = None
fnames = args.infiles or [
  "root://cms-xrd-global.cern.ch//store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root"
]

# Process nanoAOD file
p = PostProcessor(".", fnames, "nMuon>=2 && Muon_pt>5",branchsel,
                  [muSS], maxEntries=args.maxevts, postfix="_MuonSS",
                  outputbranchsel=branchsel, provenance=True,
                  prefetch=True, longTermCache=True)
p.run()
