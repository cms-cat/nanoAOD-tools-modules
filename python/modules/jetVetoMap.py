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
        #print(f"Initialized jetVMAP module with veto map '{self.veto_map_name}' and correction '{corrName}'")
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_vetoedEvent", "I", title="Flag for VetoMaps")
        #print("Branch 'Jet_vetoedEvent' initialized to store event veto flags.")
    
    def fixPhi(self, phi):
        if phi > np.pi:
            phi -= 2*np.pi
        elif phi < -np.pi:
            phi += 2*np.pi
        return phi

    def analyze(self, event):
        '''nominal “loose selection”
        - jet pT > 15 GeV
        - tight jet ID
        - jet EM fraction (charged + neutral) < 0.9
        - jets that don't overlap with PF muon (dR < 0.2)
        '''
        jets = Collection(event, "Jet")
        veto_flag = 0


        for i, jet in enumerate(jets):
            if (jet.pt> 15 and (jet.jetId ==2 or jet.jetId ==6) and jet.chEmEF < 0.9 and jet.neEmEF <0.9 and jet.muonIdx1 == -1):
                #print(f"Jet {i} passed selection criteria.")
                # Correct phi and evaluate veto map
                phi = self.fixPhi(jet.phi)
                veto_map_value = self.evaluator_VETO.evaluate(self.veto_map_name, jet.eta, phi)

                # Check if the jet is vetoed
                if veto_map_value > 0:
                    veto_flag = 1  # Set flag if a vetoed jet is found
                    break  # Break out of the loop since we only need one veto to trigger

        # Fill the branch with the veto result
        self.out.fillBranch("Jet_vetoedEvent", veto_flag)
        #print(f"Event veto status filled with value: {veto_values}")
        return True