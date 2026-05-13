from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import correctionlib
import numpy as np
from array import array
from math import pi

class jetBtag(Module):

    def __init__(self, is_mc, tagger, tagger_name, WP, json_SF, json_eff):
        """
        Module to determine if jets are b-tagged and compute b-tagging SFs."""
        self.tagger = tagger
        self.tagger_name = tagger_name
        self.is_mc = is_mc
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
        self.out.branch("Jet_isBtagged", "b", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger" % self.tagger)
        
        if self.is_mc:
            self.out.branch("Jet_btagSF", "F", lenVar="nJet", title="B-tagging scale factor for the jet")
            self.out.branch("Jet_btagSF_up_correlated", "F", lenVar="nJet", title="B-tagging scale factor for the jet, up variation, correlated across years")
            self.out.branch("Jet_btagSF_down_correlated", "F", lenVar="nJet", title="B-tagging scale factor for the jet, down variation, correlated across years")
            self.out.branch("Jet_btagSF_up_uncorrelated", "F", lenVar="nJet", title="B-tagging scale factor for the jet, up variation, uncorrelated across years")
            self.out.branch("Jet_btagSF_down_uncorrelated", "F", lenVar="nJet", title="B-tagging scale factor for the jet, down variation, uncorrelated across years")
            
            self.out.branch("Jet_isBtaggedwithSF", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF" % self.tagger)
            self.out.branch("Jet_isBtaggedwithSF_up_correlated", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF up variation, correlated across years" % self.tagger)
            self.out.branch("Jet_isBtaggedwithSF_down_correlated", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF down variation, correlated across years" % self.tagger)
            self.out.branch("Jet_isBtaggedwithSF_up_uncorrelated", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF up variation, uncorrelated across years" % self.tagger)
            self.out.branch("Jet_isBtaggedwithSF_down_uncorrelated", "O", lenVar="nJet", title="Whether the jet is b-tagged according to the %s tagger, after applying SF down variation, uncorrelated across years" % self.tagger)

    def analyze(self, event):
        jets = Collection(event, "Jet")

        isBtagged = array('B', event.nJet*[0])
        if self.is_mc:
            btagSF = array('f', event.nJet*[1.0])
            btagSF_up_correlated = array('f', event.nJet*[1.0])
            btagSF_down_correlated = array('f', event.nJet*[1.0])
            btagSF_up_uncorrelated = array('f', event.nJet*[1.0])
            btagSF_down_uncorrelated = array('f', event.nJet*[1.0])

            isBtaggedWithSF = array('B', event.nJet*[0])
            isBtaggedWithSF_up_correlated = array('B', event.nJet*[0])
            isBtaggedWithSF_down_correlated = array('B', event.nJet*[0])
            isBtaggedWithSF_up_uncorrelated = array('B', event.nJet*[0])
            isBtaggedWithSF_down_uncorrelated = array('B', event.nJet*[0])

        for ijet, jet in enumerate(jets):
            isBtagged[ijet] = getattr(jet, self.tagger_name) > self.WP

            if not self.is_mc:
                continue

            flavor = int(jet.hadronFlavour)
            abseta = float(abs(jet.eta))
            pt = float(jet.pt)

            sfEvaluator = self.evaluator_SF if abs(flavor) in [4, 5] else self.evaluator_SF_light
            isBtaggedWithSF[ijet] = isBtagged[ijet]
            isBtaggedWithSF_up_correlated[ijet] = isBtagged[ijet]
            isBtaggedWithSF_down_correlated[ijet] = isBtagged[ijet]
            isBtaggedWithSF_up_uncorrelated[ijet] = isBtagged[ijet]
            isBtaggedWithSF_down_uncorrelated[ijet] = isBtagged[ijet]

            try:
                btagSF[ijet] = sfEvaluator.evaluate("central", self.WP_name, flavor, abseta, pt)
                btagSF_up_correlated[ijet] = sfEvaluator.evaluate("up_correlated", self.WP_name, flavor, abseta, pt)
                btagSF_down_correlated[ijet] = sfEvaluator.evaluate("down_correlated", self.WP_name, flavor, abseta, pt)
                btagSF_up_uncorrelated[ijet] = sfEvaluator.evaluate("up_uncorrelated", self.WP_name, flavor, abseta, pt)
                btagSF_down_uncorrelated[ijet] = sfEvaluator.evaluate("down_uncorrelated", self.WP_name, flavor, abseta, pt)
                btagEff = self.evaluator_eff.evaluate(flavorToLabel(flavor), pt, abseta)
            except Exception:
                continue

            random = makeRandomValue(event, jet)

            isBtaggedWithSF[ijet] = applySF(isBtagged[ijet], btagSF[ijet], btagEff, random)
            isBtaggedWithSF_up_correlated[ijet] = applySF(isBtagged[ijet], btagSF_up_correlated[ijet], btagEff, random)
            isBtaggedWithSF_down_correlated[ijet] = applySF(isBtagged[ijet], btagSF_down_correlated[ijet], btagEff, random)
            isBtaggedWithSF_up_uncorrelated[ijet] = applySF(isBtagged[ijet], btagSF_up_uncorrelated[ijet], btagEff, random)
            isBtaggedWithSF_down_uncorrelated[ijet] = applySF(isBtagged[ijet], btagSF_down_uncorrelated[ijet], btagEff, random)

        self.out.fillBranch("Jet_isBtagged", isBtagged)
        if self.is_mc:
            self.out.fillBranch("Jet_btagSF", btagSF)
            self.out.fillBranch("Jet_btagSF_up_correlated", btagSF_up_correlated)
            self.out.fillBranch("Jet_btagSF_down_correlated", btagSF_down_correlated)
            self.out.fillBranch("Jet_btagSF_up_uncorrelated", btagSF_up_uncorrelated)
            self.out.fillBranch("Jet_btagSF_down_uncorrelated", btagSF_down_uncorrelated)

            self.out.fillBranch("Jet_isBtaggedwithSF", isBtaggedWithSF)
            self.out.fillBranch("Jet_isBtaggedwithSF_up_correlated", isBtaggedWithSF_up_correlated)
            self.out.fillBranch("Jet_isBtaggedwithSF_down_correlated", isBtaggedWithSF_down_correlated)
            self.out.fillBranch("Jet_isBtaggedwithSF_up_uncorrelated", isBtaggedWithSF_up_uncorrelated)
            self.out.fillBranch("Jet_isBtaggedwithSF_down_uncorrelated", isBtaggedWithSF_down_uncorrelated)

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
