#!/usr/bin/env python3
#
# Example of using the PU reweight module.
# Creates a skimmed tree and a histogram file with original and reweighed PU profiles.
#
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Set up the PU reweighting module, 2022EE
from PhysicsTools.NATModules.modules.puWeightProducer import puWeightProducer 
puWeight = puWeightProducer(json="/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2022_Summer22EE/puWeights.json.gz", key="Collisions2022_359022_362760_eraEFG_GoldenJson")


# Settings for post-processor
from argparse import ArgumentParser
parser = ArgumentParser(description="Process nanoAOD files and add branches",epilog="Good luck!")
parser.add_argument('-i', '--infiles', nargs='+')
parser.add_argument('-m', '--maxevts', type=int, default=10000) # limit number of events (per file) for testing
args = parser.parse_args()
branchsel = None
fnames = args.infiles or [
  "root://cms-xrd-global.cern.ch//store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root",
]

class plotPU(Module):
    def __init__(self):
        self.writeHistFile = True

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.h_nTrueInt = ROOT.TH1F('N_TrueInt_original', 'N_TrueInt_original', 99, 0, 99)
        self.h_nTrueIntW = ROOT.TH1F('N_TrueInt_weighed', 'N_TrueInt_weighed', 99, 0, 99)
        self.addObject(self.h_nTrueInt)
        self.addObject(self.h_nTrueIntW)
        

    def analyze(self, event):

        # fill histograms
        self.h_nTrueInt.Fill(event.Pileup_nTrueInt)
        self.h_nTrueIntW.Fill(event.Pileup_nTrueInt, event.puWeight)

        return True


# Process nanoAOD file
p = PostProcessor(".", fnames, None,
                  branchsel,
                  modules=[puWeight, plotPU()],
                  maxEntries=args.maxevts,
                  postfix="_puWeight",
                  outputbranchsel=['drop *', 'keep Pileup*' , 'keep puWeight*'], 
                  histFileName="histPU.root",
                  histDirName="plots",
                  provenance=True, prefetch=True, longTermCache=True)
p.run()

