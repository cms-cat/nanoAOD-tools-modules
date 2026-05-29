from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import correctionlib
import numpy as np
from array import array
from math import pi

class jetBtag(Module):

    def __init__(self, is_mc, tagger, tagger_name, WP, json_SF, json_eff, uncKey=None):
        """
        Module to determine if jets are b-tagged and compute b-tagging SFs. Following BTV recommendations.
        Explanation of the parameters:
       
        - is_mc is self-explanary: True or False.
        
        - tagger is the tagger name as part of the key in the b-tagging SF json.gz. Options are: "particleNet", "robustParticleTransformer" for NanoAODv12; "UParTAK4" for NanoAODv15. They are the key to obtain the working point, scale factor values, and systematic variations from the corresponding correctionlib json.gz. Check https://cms-analysis-corrections.docs.cern.ch/corrections/BTV/.
        
        - tagger_name describes the branches of the tagger score in NanoAOD. Basically, the branch will be Jet_{tagger_name}. Options are: "btagPNetB", "btagRobustParTAK4B" for NanoAODv12; "Jet_btagUParTAK4B" for NanoAODv15.
        
        - WP is the b-tagging working point. Options are "L", "M", "T", "XT", "XXT".
        
        - json_SF is the path to the json.gz of the b-tagging scale factors.
        
        - json_eff is the path to the json.gz of the b-tagging efficiency in simulation. It is required to apply b-tagging SF for untagged events. See https://btv-wiki.docs.cern.ch/PerformanceCalibration/fixedWPSFRecommendations/. This has to be computed individually based on the phase space of a certain analysis. Procedure to compute it:
            - Take all simulated processes in the signal region.
            - Apply all event-level weights properly.
            - Collect all jets in them.
            - Keep only jets passing jet ID and other jet selections applied.
            - Compute the b-tagging efficiency in bins of pt, eta, for different hadron flavors (justified by NanoAOD branch Jet_hadronFlavour).
            - Preferably compute at least separately for signal and background processes.

        - uncKey is the (list of) systematic input key(s) to be included.
            - If it is a list, a list of systematic variations namely [up_{key}, down_{key} for key in uncKey] will be computed.
            - If it is a string, then up_{uncKey} and down_{uncKey} will be computed.
            - If it is None, then total uncertainties are computed, i.e. up and down.
        """
        self.tagger = tagger
        self.tagger_name = tagger_name
        self.is_mc = is_mc
        self.uncKey = uncKey
        self.WP_name = WP

        evaluator = correctionlib.CorrectionSet.from_file(json_SF)
        self.WP = evaluator["%s_wp_values" % self.tagger].evaluate(WP)
        if self.is_mc:
            self.evaluator_SF = evaluator["%s_comb" % self.tagger]
            try:
                self.evaluator_SF_light = evaluator["%s_light" % self.tagger]
            except Exception:
                self.evaluator_SF_light = self.evaluator_SF
            self.evaluator_eff = correctionlib.CorrectionSet.from_file(json_eff)["btag_efficiency"]

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_isBtagged", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger" % self.tagger)
        
        if self.is_mc:
            self.out.branch("Jet_btagSF", "F", lenVar="nJet", title="B-tagging scale factor for the jet")
            self.out.branch("Jet_isBtaggedwithSF", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF" % self.tagger)
            if self.uncKey and isinstance(self.uncKey, list):
                for label in self.uncKey:
                    self.out.branch(f"Jet_{label}Up_btagSF", "F", lenVar="nJet", title="B-tagging scale factor for the jet, up variation, %s" % label)
                    self.out.branch(f"Jet_{label}Dn_btagSF", "F", lenVar="nJet", title="B-tagging scale factor for the jet, down variation, %s" % label)
                    self.out.branch(f"Jet_{label}Up_isBtaggedwithSF", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF up variation, %s" % (self.tagger, label))
                    self.out.branch(f"Jet_{label}Dn_isBtaggedwithSF", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF down variation, %s" % (self.tagger, label))
            elif self.uncKey:
                self.out.branch(f"Jet_{self.uncKey}Up_btagSF", "F", lenVar="nJet", title="B-tagging scale factor for the jet, up variation, %s" % self.uncKey)
                self.out.branch(f"Jet_{self.uncKey}Dn_btagSF", "F", lenVar="nJet", title="B-tagging scale factor for the jet, down variation, %s" % self.uncKey)
                self.out.branch(f"Jet_{self.uncKey}Up_isBtaggedwithSF", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF up variation, %s" % (self.tagger, self.uncKey))
                self.out.branch(f"Jet_{self.uncKey}Dn_isBtaggedwithSF", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF down variation, %s" % (self.tagger, self.uncKey))
            else:
                self.out.branch(f"Jet_Up_btagSF", "F", lenVar="nJet", title="B-tagging scale factor for the jet, up variation")
                self.out.branch(f"Jet_Dn_btagSF", "F", lenVar="nJet", title="B-tagging scale factor for the jet, down variation")
                self.out.branch(f"Jet_Up_isBtaggedwithSF", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF up variation" % (self.tagger))
                self.out.branch(f"Jet_Dn_isBtaggedwithSF", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF down variation" % (self.tagger))

    def analyze(self, event):
        jets = Collection(event, "Jet")

        isBtagged = array('B', event.nJet*[0])
        if self.is_mc:
            btagSF = array('f', event.nJet*[1.0])
            isBtaggedWithSF = array('B', event.nJet*[0])

            if self.uncKey and isinstance(self.uncKey, list):
                btagSFUp = {label: array('f', event.nJet*[1.0]) for label in self.uncKey}
                btagSFDn = {label: array('f', event.nJet*[1.0]) for label in self.uncKey}
                isBtaggedWithSFUp = {label: array('B', event.nJet*[0]) for label in self.uncKey}
                isBtaggedWithSFDn = {label: array('B', event.nJet*[0]) for label in self.uncKey}
            else:
                btagSFUp = array('f', event.nJet*[1.0])
                btagSFDn = array('f', event.nJet*[1.0])
                isBtaggedWithSFUp = array('B', event.nJet*[0])
                isBtaggedWithSFDn = array('B', event.nJet*[0])

        for ijet, jet in enumerate(jets):
            isBtagged[ijet] = getattr(jet, self.tagger_name) > self.WP

            if not self.is_mc:
                continue

            flavor = int(jet.hadronFlavour)
            abseta = float(abs(jet.eta))
            pt = float(jet.pt)

            sfEvaluator = self.evaluator_SF if abs(flavor) in [4, 5] else self.evaluator_SF_light
            isBtaggedWithSF[ijet] = isBtagged[ijet]
            if self.uncKey and isinstance(self.uncKey, list):
                for label in self.uncKey:
                    isBtaggedWithSFUp[label][ijet] = isBtagged[ijet]
                    isBtaggedWithSFDn[label][ijet] = isBtagged[ijet]
            else:
                isBtaggedWithSFUp[ijet] = isBtagged[ijet]
                isBtaggedWithSFDn[ijet] = isBtagged[ijet]

            try:
                btagSF[ijet] = sfEvaluator.evaluate("central", self.WP_name, flavor, abseta, pt)
                if self.uncKey and isinstance(self.uncKey, list):
                    for label in self.uncKey:
                        btagSFUp[label][ijet] = sfEvaluator.evaluate("up_%s" % label, self.WP_name, flavor, abseta, pt)
                        btagSFDn[label][ijet] = sfEvaluator.evaluate("down_%s" % label, self.WP_name, flavor, abseta, pt)
                elif self.uncKey:
                    btagSFUp[ijet] = sfEvaluator.evaluate("up_%s" % self.uncKey, self.WP_name, flavor, abseta, pt)
                    btagSFDn[ijet] = sfEvaluator.evaluate("down_%s" % self.uncKey, self.WP_name, flavor, abseta, pt)
                else:
                    btagSFUp[ijet] = sfEvaluator.evaluate("up", self.WP_name, flavor, abseta, pt)
                    btagSFDn[ijet] = sfEvaluator.evaluate("down", self.WP_name, flavor, abseta, pt)

                btagEff = self.evaluator_eff.evaluate(flavorToLabel(flavor), pt, abseta)

            except Exception:
                continue

            random = makeRandomValue(event, jet)

            isBtaggedWithSF[ijet] = applySF(isBtagged[ijet], btagSF[ijet], btagEff, random)
            if self.uncKey and isinstance(self.uncKey, list):
                for label in self.uncKey:
                    isBtaggedWithSFUp[label][ijet] = applySF(isBtagged[ijet], btagSFUp[label][ijet], btagEff, random)
                    isBtaggedWithSFDn[label][ijet] = applySF(isBtagged[ijet], btagSFDn[label][ijet], btagEff, random)
            else:
                isBtaggedWithSFUp[ijet] = applySF(isBtagged[ijet], btagSFUp[ijet], btagEff, random)
                isBtaggedWithSFDn[ijet] = applySF(isBtagged[ijet], btagSFDn[ijet], btagEff, random)

        self.out.fillBranch("Jet_isBtagged", isBtagged)
        if self.is_mc:
            self.out.fillBranch("Jet_btagSF", btagSF)
            self.out.fillBranch("Jet_isBtaggedwithSF", isBtaggedWithSF)
            if self.uncKey and isinstance(self.uncKey, list):
                for label in self.uncKey:
                    self.out.fillBranch(f"Jet_{label}Up_btagSF", btagSFUp[label])
                    self.out.fillBranch(f"Jet_{label}Dn_btagSF", btagSFDn[label])
                    self.out.fillBranch(f"Jet_{label}Up_isBtaggedwithSF", isBtaggedWithSFUp[label])
                    self.out.fillBranch(f"Jet_{label}Dn_isBtaggedwithSF", isBtaggedWithSFDn[label])

        return True

def flavorToLabel(hadronFlavor):
    hadronFlavor = abs(int(hadronFlavor))
    if hadronFlavor == 5:
        return "b"
    if hadronFlavor == 4:
        return "c"
    return "light"

def makeRandomValue(event, jet):
    # The seed is unique by event and jet, with a fixed entropy value to decorrelate this module.
    rng = np.random.default_rng(seed=np.random.SeedSequence([event.luminosityBlock, event.event, int(abs((jet.phi/pi)%1)*1e12), 3141592653]))
    return rng.uniform()

def applySF(isTagged, SF, eff, random):
    """Apply b-tagging SF to determine if the jet is tagged after SF application, using a random number."""
    if SF == 1:
        return isTagged
    elif SF < 1:
        if not isTagged:
            return False
        else:
            return random > 1 - SF
    else:
        if isTagged:
            return True
        else:
            if eff <= 0 or eff >= 1:
                return isTagged
            return random < (1 - SF) / (1 - 1/eff)
