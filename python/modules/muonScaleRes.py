"""Add branches for muon scale and resolution corrections.
See example in test/example_muonSS.py for usage.
"""

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
from math import pi
from ROOT import MuonScaRe

class muonScaleRes(Module):
    def __init__(self, json, is_mc, overwritePt=False, maxPt=200., minPt=26.):
        """Add branches for muon scale and resolution corrections.
        Parameters:
            json: full path of json file
            is_mc: True for MC (smear pt+add uncertainties), False for data (scale pt)
            overwritePt: replace value in the pt branch, and store the old one as "uncorrected_pt"
            maxPt: # Do not correct for muons above this pT (200 GeV according to current recipe, 8/24)
        """

        self.is_mc = is_mc
        self.overwritePt = overwritePt
        self.maxPt = maxPt
        self.minPt = minPt
        
        self.corrModule = MuonScaRe(json)


    def getPtCorr(self, muon):
        if muon.pt < self.minPt or muon.pt > self.maxPt:
            return muon.pt, muon.pt  # no correction above maxPt

        isData = int(not self.is_mc)
        scale_corr = self.corrModule.pt_scale(isData, muon.pt, muon.eta, muon.phi, muon.charge)
        #print(f"[getPtCorr] Muon pt {muon.pt:.2f} scale corrected to {scale_corr:.2f} (is_mc={self.is_mc})")

        if self.is_mc:
            smear_corr = self.corrModule.pt_resol(scale_corr, muon.eta, muon.nTrackerLayers)
            #print(f"muon.eta, muon.nTrackerLayers", muon.eta, muon.nTrackerLayers)
            #if abs(smear_corr) > 10*muon.pt:
                #print(f"muon.eta, muon.nTrackerLayers", muon.eta, muon.nTrackerLayers)
                #print(f"[getPtCorr] Muon pt {scale_corr:.2f} smeared to {smear_corr:.2f} (MC smear applied)")
            return scale_corr, smear_corr  # MC: return both
        else:
            #print(f"[getPtCorr] Data: no smearing, scale_corr used twice")
            return scale_corr, scale_corr  # Data: no smearing, return scale_corr twice

    def getPtVarRes(self, muon, pt_corr_scale, pt_corr_scaleres, updn):
        """Handle Resolution variations"""
        if muon.pt < self.minPt or muon.pt > self.maxPt:
            return muon.pt
        
        return self.corrModule.pt_resol_var(pt_corr_scale, pt_corr_scaleres, muon.eta, updn)

    def getPtVarScale(self, muon, pt_corr_scaleres, updn):
        """Handle Scale variations"""
        if muon.pt < self.minPt or muon.pt > self.maxPt:
            return muon.pt

        return self.corrModule.pt_scale_var(pt_corr_scaleres, muon.eta, muon.phi, muon.charge, updn)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Muon_pt", "F", lenVar="nMuon", title="pT (with scale/smearing corrections)")
            self.out.branch("Muon_uncorrected_pt", "F", lenVar="nMuon", title="original (uncorrected) pT")
        else:
            self.out.branch("Muon_corrected_pt", "F", lenVar="nMuon", title="pT (with scale/smearing corrections)")
        if self.is_mc:
            self.out.branch("Muon_scaleUp_pt", "F", lenVar="nMuon", title="scale uncertainty")
            self.out.branch("Muon_scaleDn_pt", "F", lenVar="nMuon", title="scale uncertainty")
            self.out.branch("Muon_smearUp_pt", "F", lenVar="nMuon", title="smearing uncertainty")
            self.out.branch("Muon_smearDn_pt", "F", lenVar="nMuon", title="smearing uncertainty")


    def analyze(self, event):
        if event.nMuon == 0 :
            return True

        muons = Collection(event, "Muon")

        pt_corr = [0.]*len(muons)
        pt_corr_scale = [0.]*len(muons)

        if self.is_mc:
            scaleUp_pt = [0.]*len(muons)
            scaleDn_pt = [0.]*len(muons)
            smearUp_pt = [0.]*len(muons)
            smearDn_pt = [0.]*len(muons)
                
        for imu, muon in enumerate(muons):
            # Set up a deterministic random seed.
            # The seed is unique by event and muon.
            # A fixed entropy value is also included to decorrelate different modules doing similar things.
            seedSeq = np.random.SeedSequence([event.luminosityBlock, event.event, int(abs((muon.phi/pi*100.)%1)*1e10), 351740215])
            self.corrModule.setSeed(int(seedSeq.generate_state(1,np.uint64)[0]))

            #print(f"\n[analyze] Muon {imu} seed set to {seed}")
            #print(f"[analyze] Muon {imu} original pt: {muon.pt:.2f}")

            #pt_corr[imu] = self.getPtCorr(muon)

            pt_corr_scale[imu], pt_corr[imu] = self.getPtCorr(muon)

            #print(f"[analyze] Muon {imu} scale corrected pt: {pt_corr_scale[imu]:.2f}")
            #print(f"[analyze] Muon {imu} final corrected pt: {pt_corr[imu]:.2f}")

            if self.is_mc:
                scaleUp_pt[imu] = self.getPtVarScale(muon, pt_corr[imu], "up")
                scaleDn_pt[imu] = self.getPtVarScale(muon, pt_corr[imu], "dn")

                smearUp_pt[imu] = self.getPtVarRes(muon, pt_corr_scale[imu], pt_corr[imu],"up")
                smearDn_pt[imu] = self.getPtVarRes(muon, pt_corr_scale[imu], pt_corr[imu],"dn")

        if self.overwritePt :
            pt_uncorr = list(mu.pt for mu in muons)
            self.out.fillBranch("Muon_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Muon_pt", pt_corr)
        else :
            self.out.fillBranch("Muon_corrected_pt", pt_corr)

        if self.is_mc:
            self.out.fillBranch("Muon_scaleUp_pt", scaleUp_pt)
            self.out.fillBranch("Muon_scaleDn_pt", scaleDn_pt)
            self.out.fillBranch("Muon_smearUp_pt", smearUp_pt)
            self.out.fillBranch("Muon_smearDn_pt", smearDn_pt)


        return True

