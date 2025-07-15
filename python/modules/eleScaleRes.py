"""Add branches for electron scale and resolution corrections.
See example in test/example_electronSS.py for usage.
"""

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib
from math import pi

class eleScaleRes(Module):
    def __init__(self, json, scaleKey=None, smearKey=None, overwritePt=False, EtDependent=False):
        """Add branches for electron scale and resolution corrections.
        Parameters:
            json: full path of json file
            scaleKey: key for the scale correction for data and scale uncertainty in MC
            smearKey: key for the MC smearing correction (None for Data)
            overwritePt: replace value in the pt branch, and store the old one as "uncorrected_pt"
        """
        if (scaleKey == None and smearKey == None) :
            raise ValueError("eleScaleResProducer: either scaleKey or smearKey should be set")
        if (smearKey != None and scaleKey == None ) :
            raise ValueError("eleScaleResProducer: scaleKey is required when smearKey is set (ie in MC)")
        self.overwritePt = overwritePt
        self.EtDependent = EtDependent
        
        evaluator = correctionlib.CorrectionSet.from_file(json)
        self.evaluator_scale = None
        self.evaluator_smear = None
        self.is_mc = False

        if self.EtDependent :
            if scaleKey != None :
                self.evaluator_scale = evaluator.compound[scaleKey]
        else:
            if scaleKey != None :
                self.evaluator_scale = evaluator[scaleKey]

        if smearKey != None :
            self.evaluator_smear = evaluator[smearKey]
            self.is_mc = True

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Electron_pt", "F", lenVar="nElectron", title="pT (with scale/smearing corrections)")
            self.out.branch("Electron_uncorrected_pt", "F", lenVar="nElectron", title="original (uncorrected) pT")
        else:
            self.out.branch("Electron_corrected_pt", "F", lenVar="nElectron", title="pT (with scale/smearing corrections)")
        if self.is_mc :
            self.out.branch("Electron_scaleUp_pt", "F", lenVar="nElectron", title="scale uncertainty")
            self.out.branch("Electron_scaleDn_pt", "F", lenVar="nElectron", title="scale uncertainty")
            self.out.branch("Electron_smearUp_pt", "F", lenVar="nElectron", title="smearing uncertainty")
            self.out.branch("Electron_smearDn_pt", "F", lenVar="nElectron", title="smearing uncertainty")

    def analyze(self, event):
        electrons = Collection(event, "Electron")

        pt_corr = []
        pt_smear_up = []
        pt_smear_dn = []
        pt_scale_up = []
        pt_scale_dn = []
        
        for ele in electrons:
            if self.is_mc :
                # Set up a deterministic random seed.
                # The seed is unique by event and electron.
                # A fixed entropy value is also included to decorrelate different modules doing similar things.
                rng = np.random.default_rng(seed=np.random.SeedSequence([event.luminosityBlock, event.event, int(abs((ele.phi/pi)%1)*1e12), 5402201385]))

                if self.EtDependent:
                    rho = self.evaluator_smear.evaluate("smear", ele.pt, ele.r9, abs(ele.eta))
                else:
                    rho = self.evaluator_smear.evaluate("rho", ele.eta, ele.r9)
                smearing = rng.normal(loc=1., scale=rho)
                pt_corr.append(smearing * ele.pt)

                if self.EtDependent:
                    unc_rho = self.evaluator_smear.evaluate("esmear", ele.pt, ele.r9, abs(ele.eta))
                else:
                    unc_rho = self.evaluator_smear.evaluate("err_rho", ele.eta, ele.r9)
                smearing_up = rng.normal(loc=1., scale=rho + unc_rho)
                smearing_dn = rng.normal(loc=1., scale=rho - unc_rho)
                pt_smear_up.append(smearing_up * ele.pt)
                pt_smear_dn.append(smearing_dn * ele.pt)

                if self.EtDependent :
                    scale_MC_unc = self.evaluator_scale.evaluate("escale", float(event.run), ele.eta, ele.r9, abs(ele.eta), ele.pt, float(ele.seedGain))
                else:
                    scale_MC_unc = self.evaluator_scale.evaluate("total_uncertainty", ele.seedGain, float(event.run), ele.eta, ele.r9, ele.pt)
                pt_scale_up.append((1+scale_MC_unc) * ele.pt)
                pt_scale_dn.append((1-scale_MC_unc) * ele.pt)
            else :
                if self.EtDependent :
                    scale = self.evaluator_scale.evaluate("scale", float(event.run), ele.eta, ele.r9, abs(ele.eta), ele.pt, float(ele.seedGain))
                else:
                    scale = self.evaluator_scale.evaluate("total_correction", ele.seedGain, float(event.run), ele.eta, ele.r9, ele.pt)
                pt_corr.append(scale * ele.pt)

        if self.overwritePt :
            pt_uncorr = list(ele.pt for ele in electrons)
            self.out.fillBranch("Electron_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Electron_pt", pt_corr)
        else :
            self.out.fillBranch("Electron_corrected_pt", pt_corr)

        if self.is_mc :
            self.out.fillBranch("Electron_smearUp_pt", pt_smear_up)
            self.out.fillBranch("Electron_smearDn_pt", pt_smear_dn)
            self.out.fillBranch("Electron_scaleUp_pt", pt_scale_up)
            self.out.fillBranch("Electron_scaleDn_pt", pt_scale_dn)

        return True
