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
    def __init__(self, json, nanoVersion=15, era=None):
        """Module to fill/update jetID variables (passTight, passTightLepVeto), packed in Jet_jetId, following recipes from https://cms-jme-jmar.docs.cern.ch/recommendations/jet-identification/#version-specific-instructions
        Parameters:
        - json: jetID json file (None for cases where ad-hoc code must be used) 
        - nanoVersion, era: determine the type of update: 
          9 (old Run2 samples) = do nothing
          12 (2022 and 2023 samples) = update existing jetId with NanoAODv12 recipe
          >12: compute using correctionlib json if json is not None; otherwise, if era corresponds to Run2, use the Run2 UltraLegacy recipe, https://cms-jme-jmar.docs.cern.ch/recommendations/jet-identification/#criteria)
        """
        self.nanoVersion=nanoVersion
        self.era=era
        if self.nanoVersion >= 15 and json is None and era in [2016, 2017, 2018] :
            # Run2 v15, use JME code
            self.evaluator = None
        elif self.nanoVersion == 12 and json is None:
            # Run3 v12, fix existing jetId with JME code
            self.evaluator = None            
        elif self.nanoVersion > 12 and json is not None :
            # Run3 v15, compute from JSON
            jetType="AK4PUPPI" # Note: Jet type is "AK4CHS" only for nano v9, which requires no correction
            self.evaluator = correctionlib.CorrectionSet.from_file(json)
            self.key_tight = f"{jetType}_Tight"
            self.key_tightLeptonVeto = f"{jetType}_TightLeptonVeto"
        elif self.nanoversion == 9:
            print("WARNING: jetId: nothing to be done for nanoAODv9; please checl the section NanoAODv9 at https://cms-jme-jmar.docs.cern.ch/recommendations/jet-identification/#version-specific-instructions")
        else :
            # all other cases are invalid
            print("ERROR: jetId: unsupported ombination of json, nanoVersion, era:", json, nanoVersion, era)
            exit(1)

    def jetId_Run2_v15(self, jet) :
        passJetIdTight = False
        passJetIdTightLepVeto = False
        abseta = abs(jet.eta)
        if self.era == 2016 :
            if abseta <= 2.4 :
                passJetIdTight = jet.neHEF < 0.9 and jet.neEmEF < 0.9 and (jet.chMultiplicity+jet.neMultiplicity) > 1 and jet.chHEF > 0.0 and jet.chMultiplicity > 0
            elif abseta > 2.4 and abseta <= 2.7 :
                passJetIdTight = jet.neHEF < 0.98 and jet.neEmEF < 0.99
            elif abseta > 2.7 and abseta <= 3.0 :
                passJetIdTight = jet.neMultiplicity >= 1
            elif abseta > 3.0 and abseta <= 5.0 :
                passJetIdTight = jet.neMultiplicity > 2 and jet.neEmEF < 0.9

            if abseta <= 2.4 :
                passJetIdTightLepVeto = passJetIdTight and jet.muEF < 0.8 and jet.chEmEF < 0.8
            else :
                passJetIdTightLepVeto = passJetIdTight              

        else: #2017, 2018
            if abseta <= 2.6 :
                passJetIdTight = jet.neHEF < 0.9 and jet.neEmEF < 0.9 and (jet.chMultiplicity+jet.neMultiplicity) > 1 and jet.chHEF > 0.0 and jet.chMultiplicity > 0
            elif abseta > 2.6 and abseta <= 2.7 :
                passJetIdTight = jet.neHEF < 0.90 and jet.neEmEF < 0.99
            elif abseta > 2.7 and abseta <= 3.0 :
                passJetIdTight = jet.neHEF < 0.9999
            elif abseta > 3.0 and abseta <= 5.0 :
                passJetIdTight = jet.neMultiplicity > 2 and jet.neEmEF < 0.9

            passJetIdTightLepVeto = False
            if abseta <= 2.7 :
                passJetIdTightLepVeto = passJetIdTight and jet.muEF < 0.8 and jet.chEmEF < 0.8
            else :
                passJetIdTightLepVeto = passJetIdTight

        return (int(passJetIdTight)*2 + int(passJetIdTightLepVeto)*4)

    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.nanoVersion >= 12 :
            docstr = "(computed using JSON)" if self.evaluator is not None else "(updated according to JME recipe)"
            self.out.branch("Jet_jetId", "b", lenVar="nJet", title=f"Jet ID flag: bit2 is tight, bit3 is tightLepVeto {docstr}") # Note: b=UChar_t


    def analyze(self, event):
        if self.nanoVersion < 12 : return True # no correction applied
        
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

        elif self.nanoVersion >= 13 and self.evaluator is not None : # jetId recomputed using correctionlib files
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

        else : # v15, evaluator is None -> Run2 recipe
            for ijet, jet in enumerate(jets):
                jet_Ids[ijet] = self.jetId_Run2_v15(jet)
                
        self.out.fillBranch("Jet_jetId", jet_Ids)
        return True
