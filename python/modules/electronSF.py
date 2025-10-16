"""
Compute electron SFs using correctionlib, and store in new branches.
Load as:
 eleSF = electronSF("POG/EGM/2016postVFP_UL/electron.json.gz")
 eleSF.addCorrection("NUM_TrackerMuons_DEN_genTracks", "2016postVFP", "sf")
 eleSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "2016postVFP", "sfdown", "sfsysdn")
 eleSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "2016postVFP", "sfup", "sfsysup")

See example in test/example_electronSF.py for details.
"""

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from correctionlib import CorrectionSet

class ElectronSF(Module):
    def __init__(self, json, collection="Electron"):
        """Electron SF correction module.
        Parameters:
            json: the correction file
            collection: name of the collection to be corrected
        Use addCorrection() to set which factors should be added.
        """
        self.collection = collection
        self.names = [ ]
        self.scenarios = [ ]
        self.wps = [ ]
        self.valtypes = [ ]
        self.varnames = [ ]
        self.evaluators = [ ]
        self.evaluator = CorrectionSet.from_file(json)
    
    def addCorrection(self, name, scenario, wp, valtype, varname=None):
        """
        Call this method to add a correction factor.
        Parameters:
            name: name of the corrections, e.g. 'UL-Electron-ID-SF'
            scenario: year/scenario, e.g. '2016postVFP'
            wp: working point
                - string, e.g. 'Loose', 'Medium', 'wp80iso', ...
                - function of pT, e.g. lambda pt: 'RecoAbove20' if pt>=20 else 'RecoBelow20'
            valtype: type of factor, e.g. 'sf', 'sfup', 'sfdown', ...
            varname: branch name suffix (defaults to {wp}_{valtype})
        """
        if varname==None: # default suffix
          assert isinstance(wp,str), "Please use the varname option to name the branch"
          varname = wp+'_'+valtype
        self.names.append(name)
        self.scenarios.append(scenario)
        self.valtypes.append(valtype)
        self.wps.append(wp)
        self.varnames.append(f"{self.collection}_{varname}") # branch name
        self.evaluators.append(self.evaluator[name])
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        """Add branch for every correction to output file."""
        self.out = wrappedOutputTree        
        for varname in set(self.varnames): # avoid duplicates
            self.out.branch(varname, 'F', lenVar='nElectron')

    def getSF(self, evaluator, valtype, scenario, wp, eta, pt, phi=None):
        varlist = [i.name for i in evaluator.inputs]

        if "phi" in varlist and phi is not None:
            return evaluator.evaluate(scenario, valtype, wp, eta, pt, phi)
        else:
            return evaluator.evaluate(scenario, valtype, wp, eta, pt)

    def analyze(self, event):
        pts = [max(10.001,event.Electron_pt[i]) for i in range(event.nElectron)]
        etas = [event.Electron_eta[i] for i in range(event.nElectron)]
        phis = [event.Electron_phi[i] for i in range(event.nElectron)]
        for ic in range(len(self.evaluators)):
            # We cannot make a single call to evaluate passing eta, pt as arrays
            # since POG JSONS are currently provided with flow="error", so we
            # have to loop to protect for values out of binning range.
            sfs = [1.]*event.nElectron
            for iEle in range(event.nElectron):
                try:
                  if isinstance(self.wps[ic],str): # WP is a simple string
                      wp = self.wps[ic]
                  else: # assume WP is a function of pT
                      wp = self.wps[ic](pts[iEle]) # evaluate WP in pT
                      if wp==None: # evaluate correction only if WP is defined
                         continue
                     
                      sfs[iEle] = self.getSF(self.evaluators[ic], self.valtypes[ic], self.scenarios[ic], wp, etas[iEle], pts[iEle], phis[iEle])
                     
                except: 
                    print(f"ElectronSF.analyze: Exception for {self.scenarios[ic]}, {self.valtypes[ic]}, wp={self.wps[ic]}, eta={etas[iEle]:6.4f}, pt={pts[iEle]:6.4f}, phi={phis[iEle]:6.4f}")
                    pass # default sf = 1
            self.out.fillBranch(self.varnames[ic], sfs)
        return True
