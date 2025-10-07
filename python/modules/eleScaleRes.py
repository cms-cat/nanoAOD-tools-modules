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
    def __init__(self, json, scaleKey=None, smearKey=None, overwritePt=False):
        """Add branches for electron scale and resolution corrections.
        Parameters:
            json: full path of json file
            scaleKey: key for the scale correction for data and scale uncertainty in MC
            smearKey: key for the MC smearing correction (None for Data; any other value implies MC)
            overwritePt: replace value in the pt branch, and store the old one as "uncorrected_pt"
        """
        if (scaleKey == None and smearKey == None) :
            raise ValueError("eleScaleResProducer: either scaleKey or smearKey should be set")
        if (smearKey != None and scaleKey == None ) :
            raise ValueError("eleScaleResProducer: scaleKey is required when smearKey is set (ie in MC)")
        self.overwritePt = overwritePt
        
        evaluator = correctionlib.CorrectionSet.from_file(json)
        self.evaluator_scale = None
        self.evaluator_smear = None
        self.is_mc = False

        if scaleKey != None :
            if scaleKey in evaluator.compound : # Et-dependent corrections. We keep track of this only as a cross-check.
                self.EtDependent = True
                self.evaluator_scale = evaluator.compound[scaleKey]
            else : # "Standard" corrections
                self.EtDependent = False 
                self.evaluator_scale = evaluator[scaleKey]

        if smearKey != None :
            self.evaluator_smear = evaluator[smearKey]
            self.is_mc = True

        # Detect which set of variables should be passed for scale (data) and smearing+uncertainties (MC)
        self.getSmear = None;
        self.getSmearUnc = None;
        self.getScale = None;
        self.getScaleUnc = None;
        varlist = [i.name for i in self.evaluator_smear.inputs]
        if varlist == ['valtype', 'eta', 'r9'] and not self.EtDependent : # Standard
            # Note: SC eta is available as ele.superclusterEta starting from nanoAODv15. It needs to be recomputed to support older versions.
            self.getSmear    = lambda ele: self.evaluator_smear.evaluate("rho",     ele.eta+ele.deltaEtaSC, ele.r9)
            self.getSmearUnc = lambda ele: self.evaluator_smear.evaluate("err_rho", ele.eta+ele.deltaEtaSC, ele.r9)
            # Note: scale corrections have to be taken from the scale module for "standard" corrections (?); see below.
        elif varlist == ['syst', 'pt', 'r9', 'AbsScEta'] and self.EtDependent :# 2022 ET-dependent
            self.getSmear    = lambda ele: self.evaluator_smear.evaluate("smear",  ele.pt, ele.r9, abs(ele.eta+ele.deltaEtaSC))
            self.getSmearUnc = lambda ele: self.evaluator_smear.evaluate("esmear", ele.pt, ele.r9, abs(ele.eta+ele.deltaEtaSC))
            self.getScaleUnc = lambda ele: self.evaluator_smear.evaluate("escale", ele.pt, ele.r9, abs(ele.eta+ele.deltaEtaSC))
        elif varlist == ['syst', 'pt', 'r9', 'ScEta'] and self.EtDependent :# 2024 ET-dependent
            self.getSmear    = lambda ele: self.evaluator_smear.evaluate("smear",  ele.pt, ele.r9, ele.eta+ele.deltaEtaSC)
            self.getSmearUnc = lambda ele: self.evaluator_smear.evaluate("esmear", ele.pt, ele.r9, ele.eta+ele.deltaEtaSC)
            self.getScaleUnc = lambda ele: self.evaluator_smear.evaluate("escale", ele.pt, ele.r9, ele.eta+ele.deltaEtaSC)
        else :
            raise ValueError(f"ERROR: eleScaleRes: unexpected inputs for key={smearKey} and EtDependent={self.EtDependent}: {varlist}")

        varlist = [i.name for i in self.evaluator_scale.inputs]        
        if varlist ==  ['valtype', 'gain', 'run', 'eta', 'r9', 'et'] and not self.EtDependent : # 2022-2023 standard
            self.getScale    = lambda event,ele: self.evaluator_scale.evaluate("total_correction",  ele.seedGain, float(event.run), ele.eta+ele.deltaEtaSC, ele.r9, ele.pt)
            self.getScaleUnc = lambda ele: self.evaluator_scale.evaluate("total_uncertainty", ele.seedGain, 1., ele.eta+ele.deltaEtaSC, ele.r9, ele.pt) # Note: Run is always 1 in MC
        elif varlist ==  ['syst', 'run', 'ScEta', 'r9', 'AbsScEta', 'pt', 'seedGain'] and self.EtDependent : # 2022-2023 ET-dependent
            self.getScale    = lambda event,ele : self.evaluator_scale.evaluate("scale",  float(event.run), ele.eta+ele.deltaEtaSC, ele.r9, abs(ele.eta+ele.deltaEtaSC), ele.pt, float(ele.seedGain))
        elif varlist ==  ['syst', 'run', 'ScEta', 'r9', 'pt', 'seedGain'] and self.EtDependent : # 2024 ET-dependent
            self.getScale    = lambda event,ele : self.evaluator_scale.evaluate("scale",  float(event.run), ele.eta+ele.deltaEtaSC, ele.r9, ele.pt, float(ele.seedGain))
            
        else :
            raise ValueError(f"ERROR: eleScaleRes: unexpected inputs for key={scaleKey} and EtDependent={self.EtDependent}: : {varlist}")

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
            if self.is_mc : # smearing and uncertainties for MC
                # Set up a deterministic random seed.
                # The seed is unique by event and electron.
                # A fixed entropy value is also included to decorrelate different modules doing similar things.
                rng = np.random.default_rng(seed=np.random.SeedSequence([event.luminosityBlock, event.event, int(abs((ele.phi/pi)%1)*1e12), 5402201385]))
                rnd = rng.normal(loc=0., scale=1) # The same random number should be used for central value and up/down variation

                smear = self.getSmear(ele)
                smeared_pt = ele.pt*(1.+smear*rnd)
                pt_corr.append(smeared_pt)

                unc_smear = self.getSmearUnc(ele)
                pt_smear_up.append(ele.pt*(1.+(smear+unc_smear)*rnd))
                pt_smear_dn.append(ele.pt*(1.+max(0.,smear-unc_smear)*rnd))

                # Scale uncertainty is evaluated from original variables, but applied to the smeared pT.
                scale_MC_unc = self.getScaleUnc(ele)
                pt_scale_up.append((1+scale_MC_unc) * smeared_pt)
                pt_scale_dn.append((1-scale_MC_unc) * smeared_pt)
                
            else : # Scale correction for data
                scale = self.getScale(event, ele)
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
