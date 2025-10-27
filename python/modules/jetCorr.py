"""Add branches for muon scale and resolution corrections.
See example in test/example_jetCorr.py for usage.
"""

###
# Useful CMStalk thread on how implementing jet corrections: https://cms-talk.web.cern.ch/t/jes-for-2022-re-reco-cde-and-prompt-fg/32873/3
# Minimal demo provided from JME POG: https://github.com/cms-jet/JECDatabase/blob/master/scripts/JERC2JSON/minimalDemo.py
# TODO:
# - Implementing JES uncertainties. TBD which scheme we want to implement
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib

class jetJERC(Module):
    def __init__(self, json_JERC, json_JERsmear, L1Key=None, L2Key=None, L3Key=None, L2L3Key=None, scaleTotalKey=None, scaleKeyRegrouped11=None, smearKey=None, JERKey=None, JERsfKey=None, overwritePt=False, usePhiDependentJEC=False, useRunDependentJEC=False, useJesSplittingScheme11=False):
        """Correct jets following recommendations of JME POG.
        Parameters:
            json_JERC: full path of json file with JERC corrections
            json_JERsmear: full path of json file with smearing terms for JER
            L1Key: key for L1 corrections
            L2Key: key for L2 corrections
            L3Key: key for L3 corrections
            L2L3Key: key for residual corrections
            scaleTotalKey : key for JES uncertainties (1 source, Total)
            scaleKeyRegrouped11: key for JES uncertainties splitted in 11 sources
            smearKey: key for smearing formula (None for Data)
            JERKey: key for JER (None for Data)
            JERsfKey: key for JER scale factor (None for Data)
            overwritePt: replace value in the pt branch, and store the old one as "uncorrected_pt"
        """
        self.overwritePt = overwritePt
        self.usePhiDependentJEC = usePhiDependentJEC
        self.useRunDependentJEC = useRunDependentJEC
        self.useJesSplittingScheme11 = useJesSplittingScheme11

        self.evaluator_JERC = correctionlib.CorrectionSet.from_file(json_JERC)
        self.evaluator_jer = correctionlib.CorrectionSet.from_file(json_JERsmear)

        ## Starting from Run 3, the use of PUPPI jets eliminates the need for L1 corrections. To maintain compatibility with Run 2 scripts, a dummy file is provided.
        self.evaluator_L1 = self.evaluator_JERC[L1Key]
        ## For the next step, there is a mismatch between the terminology in the twiki and the tags in correctionlib
        ## MC-truth = L2 + L3
        self.evaluator_L2 = self.evaluator_JERC[L2Key]
        self.evaluator_L3 = self.evaluator_JERC[L3Key]
        ## L2L3residuals should be applied only to data
        ## A dummy file with all entries equal to 1 is provided for MC to have a common script for both data and MC
        self.evaluator_L2L3 = self.evaluator_JERC[L2L3Key]

        self.is_mc = False
        self.evaluator_JERsmear = None
        self.evaluator_JER = None
        self.evaluator_JERsf = None
        self.evaluator_JES = None
        if smearKey != None: ## JER is applied only to MC
            self.evaluator_JERsmear = self.evaluator_jer[smearKey]
            self.evaluator_JER = self.evaluator_JERC[JERKey]
            self.evaluator_JERsf = self.evaluator_JERC[JERsfKey]
            if self.useJesSplittingScheme11:
                # For 11-source splitting, create a dictionary of evaluators
                self.evaluator_JES = {}
                for full_key in scaleKeyRegrouped11:
                    label = full_key.split("_Regrouped_")[-1]
                    if label.endswith("_AK4PFPuppi"):
                        label = label[:-len("_AK4PFPuppi")]
                    label = label.replace("-", "_").replace(".", "_")
                    self.evaluator_JES[label] = self.evaluator_JERC[full_key]
            else:
                self.evaluator_JES = self.evaluator_JERC[scaleTotalKey]
            self.is_mc = True

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Jet_pt", "F", lenVar="nJet", title="pT (with JES/JER corrections)")
            self.out.branch("Jet_mass", "F", lenVar="nJet", title="mass (with JES/JER corrections)", limitedPrecision=12)
            self.out.branch("Jet_uncorrected_pt", "F", lenVar="nJet", title="original (uncorrected) pT")
            self.out.branch("Jet_uncorrected_mass", "F", lenVar="nJet", title="original (uncorrected) mass", limitedPrecision=12)
        else:
            self.out.branch("Jet_corrected_pt", "F", lenVar="nJet", title="pT (with JES/JER corrections)")
            self.out.branch("Jet_corrected_mass", "F", lenVar="nJet", title="mass (with JES/JER corrections)", limitedPrecision=12)
        if self.is_mc:
            self.out.branch("Jet_smearUp_pt", "F", lenVar="nJet", title="smearing uncertainty")
            self.out.branch("Jet_smearDn_pt", "F", lenVar="nJet", title="smearing uncertainty")
            self.out.branch("Jet_smearUp_mass", "F", lenVar="nJet", title="smearing uncertainty", limitedPrecision=12)
            self.out.branch("Jet_smearDn_mass", "F", lenVar="nJet", title="smearing uncertainty", limitedPrecision=12)
            if self.useJesSplittingScheme11:  # 11-source splitting
                for label in self.evaluator_JES.keys():
                    self.out.branch(f"Jet_{label}_ScaleUp_pt", "F", lenVar="nJet")
                    self.out.branch(f"Jet_{label}_ScaleDn_pt", "F", lenVar="nJet")
                    self.out.branch(f"Jet_{label}_ScaleUp_mass", "F", lenVar="nJet", limitedPrecision=12)
                    self.out.branch(f"Jet_{label}_ScaleDn_mass", "F", lenVar="nJet", limitedPrecision=12)
            else:
                self.out.branch("Jet_scaleUp_pt", "F", lenVar="nJet", title="scale uncertainty")
                self.out.branch("Jet_scaleDn_pt", "F", lenVar="nJet", title="scale uncertainty")
                self.out.branch("Jet_scaleUp_mass", "F", lenVar="nJet", title="scale uncertainty", limitedPrecision=12)
                self.out.branch("Jet_scaleDn_mass", "F", lenVar="nJet", title="scale uncertainty", limitedPrecision=12)

    def fixPhi(self, phi):
        if phi > np.pi:
            phi -= 2*np.pi
        elif phi < -np.pi:
            phi += 2*np.pi
        return phi

    def analyze(self, event):
        jets = Collection(event, "Jet")
        if self.is_mc: ## genJet info is necessary for JER
            gen_jets = Collection(event, "GenJet")
            gen_jets_pt = np.array([gen_jet.pt for gen_jet in gen_jets])
            gen_jets_eta = np.array([gen_jet.eta for gen_jet in gen_jets])
            gen_jets_phi = np.array([gen_jet.phi for gen_jet in gen_jets])

        pt_corr = []
        pt_uncorr = []
        mass_corr = []
        mass_uncorr = []
        pt_smear_up = []
        pt_smear_dn = []
        mass_smear_up = []
        mass_smear_dn = []
        pt_scale_up = []
        pt_scale_dn = []
        mass_scale_up = []
        mass_scale_dn = []

        # For 11-source JES, create lists of dicts
        if self.useJesSplittingScheme11:
            pt_JES_up_list = []
            pt_JES_dn_list = []
            mass_JES_up_list = []
            mass_JES_dn_list = []

        for jet in jets:
            #### JEC ####
            ## To be applied to both data and MC
            ## Jet in NanoAOD are already corrected
            ## The correction should be removed and the latest one available should be applied
            pt_raw = jet.pt * (1 - jet.rawFactor)
            mass_raw = jet.mass * (1 - jet.rawFactor)
            ## The three steps of JEC corrections are provided separately
            pt_L1 = pt_raw * self.evaluator_L1.evaluate(jet.area, jet.eta, pt_raw, event.Rho_fixedGridRhoFastjetAll)

            if self.usePhiDependentJEC:
                pt_L2 = pt_L1 * self.evaluator_L2.evaluate(jet.eta, jet.phi, pt_L1)
            else:
                pt_L2 = pt_L1 * self.evaluator_L2.evaluate(jet.eta, pt_L1)

            pt_L3 = pt_L2 * self.evaluator_L3.evaluate(jet.eta, pt_L2)

            if self.useRunDependentJEC:
                pt_JEC = pt_L3 * self.evaluator_L2L3.evaluate(float(event.run), jet.eta, pt_L3)
            else:
                pt_JEC = pt_L3 * self.evaluator_L2L3.evaluate(jet.eta, pt_L3)

            JEC = pt_JEC / pt_raw
            mass_JEC = mass_raw * JEC

            if self.is_mc:
                #### JER ####
                ## Hybrid method is implemented [https://cms-jerc.web.cern.ch/JER/#smearing-procedures]
                ## except in the eta region 2.5 < |eta| < 3.0 where the Scaling method is implemented
                ## Scaling method was introduced following JME recommendations to handle the 'Horns issue'
                ## Ref: https://gitlab.cern.ch/cms-jetmet/coordination/coordination/-/issues/113
                ## TODO: This should be removed once a JSON-level fix is implemented.

                JER = self.evaluator_JER.evaluate(jet.eta, pt_JEC, event.Rho_fixedGridRhoFastjetAll)
                ## GenMatching with genJet
                delta_eta = jet.eta - gen_jets_eta
                fixPhi = np.vectorize(self.fixPhi, otypes=[float])
                delta_phi = fixPhi(jet.phi - gen_jets_phi)
                pt_gen = np.where((np.abs(pt_JEC - gen_jets_pt) < 3 * pt_JEC * JER) & (np.sqrt(delta_eta**2 + delta_phi**2)<0.2), gen_jets_pt, -1.0)
                pt_gen = pt_gen[pt_gen > 0][0] if np.any(pt_gen > 0) else -1. ## If no gen-matching, simply -1

                JERsf = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "nom")
                JERsf_up = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "up")
                JERsf_dn = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "down")
                # FIXME: prevent error where event is outside the int32 range. cf: github.com/cms-nanoAOD/correctionlib/issues/298
                JERsmear = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, int(event.event&0x7FFFFFFF), JER, JERsf)
                JERsmear_up = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, int(event.event&0x7FFFFFFF), JER, JERsf_up)
                JERsmear_dn = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, int(event.event&0x7FFFFFFF), JER, JERsf_dn)

                #JES
                if self.useJesSplittingScheme11:  # 11-source JES
                    JES_up = {}
                    JES_dn = {}
                    JES_up_mass = {}
                    JES_dn_mass = {}
                    for label, evaluator in self.evaluator_JES.items():
                        u = evaluator.evaluate(jet.eta, pt_JEC)
                        JES_up[label] = pt_JEC * (1 + u)
                        JES_dn[label] = pt_JEC * (1 - u)
                        JES_up_mass[label] = mass_JEC * (1 + u)
                        JES_dn_mass[label] = mass_JEC * (1 - u)
                    pt_JES_up_list.append(JES_up)
                    pt_JES_dn_list.append(JES_dn)
                    mass_JES_up_list.append(JES_up_mass)
                    mass_JES_dn_list.append(JES_dn_mass)
                else:  # single total JES
                    u = self.evaluator_JES.evaluate(jet.eta, pt_JEC)
                    pt_JES_up = pt_JEC * (1 + u)
                    pt_JES_dn = pt_JEC * (1 - u)
                    mass_JES_up = mass_JEC * (1 + u)
                    mass_JES_dn = mass_JEC * (1 - u)

                if pt_gen < 0 and (2.5 < abs(jet.eta) < 3):
                    JERsmear_nominal = 1.0
                    JERsmear_up = JERsmear_up / JERsmear
                    JERsmear_dn = JERsmear_dn / JERsmear
                else:
                    JERsmear_nominal = JERsmear

                pt_JEC_JER = pt_JEC * JERsmear_nominal
                pt_JEC_JER_up = pt_JEC * JERsmear_up
                pt_JEC_JER_dn = pt_JEC * JERsmear_dn
                mass_JEC_JER = mass_JEC * JERsmear_nominal
                mass_JEC_JER_up = mass_JEC * JERsmear_up
                mass_JEC_JER_dn = mass_JEC * JERsmear_dn

                pt_corr.append(pt_JEC_JER)
                pt_uncorr.append(pt_raw)
                mass_corr.append(mass_JEC_JER)
                mass_uncorr.append(mass_raw)
                pt_smear_up.append(pt_JEC_JER_up)
                pt_smear_dn.append(pt_JEC_JER_dn)
                mass_smear_up.append(mass_JEC_JER_up)
                mass_smear_dn.append(mass_JEC_JER_dn)

                if not self.useJesSplittingScheme11:
                    pt_scale_up.append(pt_JES_up)
                    pt_scale_dn.append(pt_JES_dn)
                    mass_scale_up.append(mass_JES_up)
                    mass_scale_dn.append(mass_JES_dn)
            else:
                ## Data
                ## No JER for Data
                pt_corr.append(pt_JEC)
                pt_uncorr.append(pt_raw)
                mass_corr.append(mass_JEC)
                mass_uncorr.append(mass_raw)

        if self.overwritePt :
            self.out.fillBranch("Jet_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Jet_pt", pt_corr)
            self.out.fillBranch("Jet_uncorrected_mass", mass_uncorr)
            self.out.fillBranch("Jet_mass", mass_corr)
        else :
            self.out.fillBranch("Jet_corrected_pt", pt_corr)
            self.out.fillBranch("Jet_corrected_mass", mass_corr)

        if self.is_mc:
            self.out.fillBranch("Jet_smearUp_pt", pt_smear_up)
            self.out.fillBranch("Jet_smearDn_pt", pt_smear_dn)
            self.out.fillBranch("Jet_smearUp_mass", mass_smear_up)
            self.out.fillBranch("Jet_smearDn_mass", mass_smear_dn)
            # Fill JES branches
            if self.useJesSplittingScheme11:  # 11-source JES
                for label in self.evaluator_JES.keys():
                    self.out.fillBranch(f"Jet_{label}_ScaleUp_pt", [JES_up[label] for JES_up in pt_JES_up_list])
                    self.out.fillBranch(f"Jet_{label}_ScaleDn_pt", [JES_dn[label] for JES_dn in pt_JES_dn_list])
                    self.out.fillBranch(f"Jet_{label}_ScaleUp_mass", [JES_up[label] for JES_up in mass_JES_up_list])
                    self.out.fillBranch(f"Jet_{label}_ScaleDn_mass", [JES_dn[label] for JES_dn in mass_JES_dn_list])
            else:  # single total JES
                self.out.fillBranch("Jet_scaleUp_pt", pt_scale_up)
                self.out.fillBranch("Jet_scaleDn_pt", pt_scale_dn)
                self.out.fillBranch("Jet_scaleUp_mass", mass_scale_up)
                self.out.fillBranch("Jet_scaleDn_mass", mass_scale_dn)

        return True
