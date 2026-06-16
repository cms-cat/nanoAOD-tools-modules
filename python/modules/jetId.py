"""
Add or update Jet_jetId variable, following the JME recipes.
See example in test/example_jetId.py for usage.
"""
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib
from array import array

class jetId(Module):
    def __init__(self, json,  nanoVersion=15):
        """Module to fill/update jetID variables (passTight, passTightLepVeto), packed in Jet_jetId.
        Parameters:
        - json: jetID json file (relevant only for nanoVersion>9)
        - nanoVersion: determines the type of update: 
          9   = no correction (https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL#NanoAODv9)
          12  = update existing jetId according to recipe at https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV#Jet_ID_implementation_with_NanoA
          >12 = compute using JME correctionlib json as per https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/jetidExample.py
        """
        self.nanoVersion=nanoVersion
        if self.nanoVersion > 12:
            jetType="AK4PUPPI" # Note: Jet type is "AK4CHS" only for nano v9, which requires no correction
            self.evaluator = correctionlib.CorrectionSet.from_file(json)
            self.key_tight = f"{jetType}_Tight"
            self.key_tightLeptonVeto = f"{jetType}_TightLeptonVeto"
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.nanoVersion > 9 :
            docstr = "(recomputed using JSON)" if self.nanoVersion > 12 else "(updated according to nanoAODv12 recipe)"
            self.out.branch("Jet_jetId", "b", lenVar="nJet", title=f"Jet ID flag: bit2 is tight, bit3 is tightLepVeto {docstr}") # Note: b=UChar_t


    def analyze(self, event):
        if self.nanoVersion < 12 : return True # no correction needed
        
        jets = Collection(event, "Jet")
        jet_Ids = array('B', event.nJet*[0]) # Note: UChar_t is uppercase 'B' in python array

        if self.nanoVersion == 12 : # Update existing jetId variable 
            passTight = False
            passTightLepVeto = False
            for ijet, jet in enumerate(jets):
                passTightOrig = bool(jet.jetId & (1 << 1))
                # Jet-passJetIdTight based on eta conditions
                if abs(jet.eta) <= 2.7:
                    passTight = passTightOrig
                elif 2.7 < abs(jet.eta) <= 3.0:
                    passTight = passTightOrig and (jet.neHEF < 0.99)
                elif abs(jet.eta) > 3.0:
                    passTight = passTightOrig and (jet.neEmEF < 0.4)

                # Jet-passJetIdTightLepVeto based on additional lepton veto conditions
                if abs(jet.eta) <= 2.7:
                    passTightLepVeto = passTight and (jet.muEF < 0.8) and (jet.chEmEF < 0.8)
                else:
                    passTightLepVeto = passTight

                jet_Ids[ijet] = int(passTight)*2 + int(passTightLepVeto)*4

        elif self.nanoVersion >= 13 : # jetId recomputed using correctionlib files
            for ijet, jet in enumerate(jets):
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

                jet_Ids[ijet] = int(passTight)*2 + int(passTightLepVeto)*4

        self.out.fillBranch("Jet_jetId", jet_Ids)
        return True
