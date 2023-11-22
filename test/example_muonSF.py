#!/usr/bin/env python3
###
# This is an example of configuring and calling a generic correctionlib module.
# Notes:
# - For the sake of speed, it is important to apply a preselection whenever possible.
# - Also, it is advisable to insert correction modules like this after modules
#   that apply selections, so that corrections are not computed for events that
#   are then discarded.
# - Find all available corrections and parameters by running in the command
#     correction summary POG/MUO/2016postVFP_UL/muon_Z.json.gz
###
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.muonSF import *

# Set up the muon correction module
muSF = MuonSF("data/POG/MUO/2016postVFP_UL/muon_Z.json.gz")
era = '2016postVFP_UL'

# Add TrackerMuon reconstruction scale factor
muSF.addCorrection('NUM_TrackerMuons_DEN_genTracks', era, 'sf')
# Add Medium ID scale factor
muSF.addCorrection('NUM_MediumID_DEN_TrackerMuons', era, 'sf')
# Add Medium ID scale factor, down variation
muSF.addCorrection('NUM_MediumID_DEN_TrackerMuons', era, 'systdown', 'sfsysdn')
# Add Medium ID scale factor, up variation
muSF.addCorrection('NUM_MediumID_DEN_TrackerMuons', era, 'systup', 'sfsysup')

# Settings for post-processor
from argparse import ArgumentParser
parser = ArgumentParser(description="Process nanoAOD files and add branches",epilog="Good luck!")
parser.add_argument('-i', '--infiles', nargs='+')
parser.add_argument('-m', '--maxevts', type=int, default=10000) # limit number of events (per file) for testing
args = parser.parse_args()
branchsel = None #"keepElectron.txt" # keep only Electron branches for speed
fnames = args.infiles or [
  "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16NanoAODv9/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/20UL16JMENano_106X_mcRun2_asymptotic_v17-v1/2820000/11061525-9BB6-F441-9C12-4489135219B7.root"
]

# Process nanoAOD file
p = PostProcessor(".", fnames, "nMuon>=2 && Muon_pt>18", branchsel, [muSF], maxEntries=args.maxevts, postfix="_MuonSFs",
                  outputbranchsel=branchsel, provenance=True, prefetch=True, longTermCache=True)
p.run()