"""
Add jetId variable, following the example in https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/jetidExample.py .
See example in test/example_jetId.py for usage.
"""
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib

class jetId(Module):
    def __init__(self, json,  jetType="AK4PUPPI"):
        """Module to determine jetID variables (passTight, passTightLepVeto), 
        packed in Jet_jetId as in nanoAODv12
        Parameters:
        - json: jetID json file
        - jetType: "AK4PUPPI" or "AK4CHS"
        """
        self.evaluator = correctionlib.CorrectionSet.from_file(json)
        self.key_tight = f"{jetType}_Tight"
        self.key_tightLeptonVeto = f"{jetType}_TightLeptonVeto"
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_jetId", "B", title="Jet ID flag: bit2 is tight, bit3 is tightLepVeto (recomputed using JSON)")


    def analyze(self, event):
        jets = Collection(event, "Jet")

        for jet in jets:
            multiplicity = jet.chMultiplicity + jet.neMultiplicity

            passTight = self.evaluator[self.key_tight].evaluate(
                jet.eta,
                jet.chHEF,
                jet.neHEF,
                jet.chEmEF,
                jet.neEmEF,
                jet.muEF,
                jet.chMultiplicity,
                jet.neMultiplicity,
                multiplicity
            )

            passTightLepVeto = self.evaluator[self.key_tightLeptonVeto].evaluate(
                jet.eta,
                jet.chHEF,
                jet.neHEF,
                jet.chEmEF,
                jet.neEmEF,
                jet.muEF,
                jet.chMultiplicity,
                jet.neMultiplicity,
                multiplicity
            )

        # Fill the branch with the veto result
        self.out.fillBranch("Jet_jetId", int(passTight)*2 + int(passTightLepVeto)*4)
        return True
