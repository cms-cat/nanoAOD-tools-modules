#!/usr/bin/env python3
"""This is an example of configuring and using the b-tagging SF module."""

import os
from argparse import ArgumentParser

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.jetBtag import jetBtag


def get_btag_config(era, tag):
    if era != 2022:
        raise ValueError("example_jetBtag.py currently has local efficiency JSONs only for 2022/2022EE")

    tagger = "particleNet"
    tagger_name = "btagPNetB"
    cmssw_base = os.environ.get("CMSSW_BASE")
    if cmssw_base is None:
        cmssw_base = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
    data_dir = os.path.join(cmssw_base, "src", "ZZAnalysis", "NanoAnalysis", "data", "btagEff")

    if "pre_EE" in tag:
        json_SF = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22/btagging.json.gz"
        json_eff = os.path.join(data_dir, "btag_2022.json.gz")
    else:
        json_SF = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22EE/btagging.json.gz"
        json_eff = os.path.join(data_dir, "btag_2022EE.json.gz")

    return tagger, tagger_name, json_SF, json_eff


parser = ArgumentParser(description="Process nanoAOD files and add b-tagging branches", epilog="Good luck!")
parser.add_argument("-i", "--infiles", nargs="+")
parser.add_argument("-m", "--maxevts", type=int, default=1000)
parser.add_argument("--era", type=int, default=2022)
parser.add_argument("--tag", default="post_EE", choices=["pre_EE", "post_EE"])
parser.add_argument("--wp", default="M", choices=["L", "M", "T", "XT", "XXT"])
parser.add_argument("--data", default=False, action="store_true", help="Run without applying MC SF reshuffling branches")
args = parser.parse_args()

tagger, tagger_name, json_SF, json_eff = get_btag_config(args.era, args.tag)
jet_btag = jetBtag(
    is_mc=not args.data,
    tagger=tagger,
    tagger_name=tagger_name,
    WP=args.wp,
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
