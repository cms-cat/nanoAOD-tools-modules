#!/usr/bin/env python3
"""This is an example of configuring and using the module for correcting jets."""

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.jetCorr import jetJERC

json_JERC = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/2025-09-23/jet_jerc.json.gz"
json_JERsmear = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz"

jes_systematics_11split = [
        "Regrouped_Absolute",
        "Regrouped_Absolute_{year}",
        "Regrouped_BBEC1",
        "Regrouped_BBEC1_{year}",
        "Regrouped_EC2",
        "Regrouped_EC2_{year}",
        "Regrouped_FlavorQCD",
        "Regrouped_HF",
        "Regrouped_HF_{year}",
        "Regrouped_RelativeBal",
        "Regrouped_RelativeSample_{year}",
    ]

L1Key = "Summer22_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi"
L2Key = "Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi"
L3Key = "Summer22_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi"
L2L3Key = "Summer22_22Sep2023_V3_MC_L2L3Residual_AK4PFPuppi"
scaleTotalKey = "Summer22_22Sep2023_V3_MC_Total_AK4PFPuppi"
scaleKeyRegrouped11 = [f"Summer22_22Sep2023_V3_MC_{label.format(year='2022')}_AK4PFPuppi" for label in jes_systematics_11split] # 11-source JES uncertainties
smearKey = "JERSmear" # 1 Total source JES uncertainties
JERKey = "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
JERsfKey = "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
overwritePt = True
usePhiDependentJEC = False
useRunDependentJEC = False
useJesSplittingScheme11 = False
scaleKey = scaleKeyRegrouped11 if useJesSplittingScheme11 else scaleTotalKey

jetCorrected = jetJERC(json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key, scaleKey, smearKey, JERKey, JERsfKey, overwritePt, usePhiDependentJEC, useRunDependentJEC)

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
p = PostProcessor(".", fnames, "nJet>0 && Jet_pt>30 && abs(Jet_eta)<4.7",
                  branchsel,
                  [jetCorrected],
                  maxEntries=args.maxevts,
                  postfix="_jetCorr",
                  outputbranchsel=branchsel,
                  provenance=True, prefetch=True, longTermCache=True)
p.run()
