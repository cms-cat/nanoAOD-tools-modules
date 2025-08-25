"""Add jet veto map flag to jets.
See example in test/example_jetVeto.py for usage.
"""
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib

class jetVMAP(Module):
    def __init__(self, json_JVMAP, corrName=None, veto_map_name = "jetvetomap"):
        """Module to apply veto maps that veto out events with important jets in the "hot" or 
        "cold" zones. These maps should be applied similarly both on Data and MC, to keep the
        phase-spaces equal. 
        Parameters:
        """
        self.veto_map_name = veto_map_name
        self.evaluator_JVMAP= correctionlib.CorrectionSet.from_file(json_JVMAP)
        self.evaluator_VETO = self.evaluator_JVMAP[corrName]
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Flag_JetVetoed", "O", title="Event veto flag from Jet Veto Map") #FIXME
    
    def fixPhi(self, phi):
        epsilon = 1e-6  # Small offset to avoid boundary issues
        if phi > np.pi:
            #print(f"phi {phi} is greater than pi. Setting phi to pi - epsilon.")
            phi = np.pi - epsilon
        elif phi < -np.pi:
            #print(f"phi {phi} is less than -pi. Setting phi to -pi + epsilon.")
            phi = -np.pi + epsilon
        return phi

    def analyze(self, event):
        '''nominal “loose selection”
        - jet pT > 15 GeV
        - tight jet ID
        - jet EM fraction (charged + neutral) < 0.9
        - jets that don't overlap with PF muon (dR < 0.2)
        '''
        jets = Collection(event, "Jet")
        veto_flag = False


        for i, jet in enumerate(jets):
            if (jet.pt> 15 and (jet.jetId ==2 or jet.jetId ==6) and (jet.chEmEF + jet.neEmEF)<0.9 and jet.muonIdx1 == -1 and jet.muonIdx2 == -1):

                # Correct phi and evaluate veto map
                phi = self.fixPhi(jet.phi)
                veto_map_value = self.evaluator_VETO.evaluate(self.veto_map_name, jet.eta, phi)

                # Check if the jet is vetoed
                if veto_map_value > 0:
                    veto_flag = True  # Set flag if a vetoed jet is found
                    break  # Break out of the loop since we only need one veto to trigger

        # Fill the branch with the veto result
        self.out.fillBranch("Flag_JetVetoed", veto_flag)
        return True
