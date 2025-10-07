#!/usr/bin/env python3
"""This is an example of configuring and calling a generic correctionlib module.

Notes:
- For the sake of speed, it is important to apply a preselection whenever possible.
- Also, it is advisable to insert correction modules like this after modules
  that apply selections, so that corrections are not computed for events that
  are then discarded.
- Find all available corrections and parameters by running in the command
    correction summary POG/MUO/2016postVFP_UL/muon_Z.json.gz
"""

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.muonSF import *

# Set up the muon correction module
muSF = MuonSF("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-22EFGSep23-Summer22EE-NanoAODv12/2025-08-14/muon_Z.json.gz")

# Add Medium ID scale factor
muSF.addCorrection('NUM_MediumID_DEN_TrackerMuons', 'nominal', "SF")
# Add Medium ID scale factor, syst
muSF.addCorrection('NUM_MediumID_DEN_TrackerMuons', 'syst', 'SFsyst')

# Settings for post-processor
from argparse import ArgumentParser
parser = ArgumentParser(description="Process nanoAOD files and add branches",epilog="Good luck!")
parser.add_argument('-i', '--infiles', nargs='+')
parser.add_argument('-m', '--maxevts', type=int, default=10000) # limit number of events (per file) for testing
args = parser.parse_args()
branchsel = None #"keepElectron.txt" # keep only Electron branches for speed
fnames = args.infiles or [
    "root://cms-xrd-global.cern.ch//store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root", # 2022EE

]

# Process nanoAOD file
p = PostProcessor(".", fnames, "nMuon>=2 && Muon_pt>18", branchsel, [muSF], maxEntries=args.maxevts, postfix="_MuonSFs",
                  outputbranchsel=branchsel, provenance=True, prefetch=True, longTermCache=True)
p.run()
