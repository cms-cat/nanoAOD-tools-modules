#!/usr/bin/env python3
"""This is an example of configuring and using the b-tagging SF module."""

import os
from argparse import ArgumentParser

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.jetBtag import jetBtag


tagger = "particleNet"
tagger_name = "btagPNetB"
WP = "M"
is_mc = True
json_SF = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22EE/btagging.json.gz"
json_eff = os.path.join(os.path.dirname(__file__), "btag_2022EE.json.gz")


parser = ArgumentParser(description="Process nanoAOD files and add b-tagging branches", epilog="Good luck!")
parser.add_argument("-i", "--infiles", nargs="+")
parser.add_argument("-m", "--maxevts", type=int, default=1000)
args = parser.parse_args()

jet_btag = jetBtag(
    is_mc=is_mc,
    tagger=tagger,
    tagger_name=tagger_name,
    WP=WP,
    json_SF=json_SF,
    json_eff=json_eff,
)

branchsel = None
fnames = args.infiles or [
    "root://cms-xrd-global.cern.ch//store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root",
]

p = PostProcessor(
    ".",
    fnames,
    "nJet>0 && Jet_pt>30 && abs(Jet_eta)<2.5",
    branchsel,
    [jet_btag],
    maxEntries=args.maxevts,
    postfix="_jetBtag",
    outputbranchsel=branchsel,
    provenance=True,
    prefetch=True,
    longTermCache=True,
)
p.run()
