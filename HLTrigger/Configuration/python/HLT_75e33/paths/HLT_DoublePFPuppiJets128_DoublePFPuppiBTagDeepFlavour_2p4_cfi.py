import FWCore.ParameterSet.Config as cms

from ..modules.hltBTagPFPuppiDeepFlavour0p935DoubleEta2p4_cfi import *
from ..modules.hltDoublePFPuppiJets128Eta2p4MaxDeta1p6_cfi import *
from ..modules.hltDoublePFPuppiJets128MaxEta2p4_cfi import *
from ..modules.l1tDoublePFPuppiJet112offMaxEta2p4_cfi import *
from ..modules.l1tDoublePFPuppiJets112offMaxDeta1p6_cfi import *
from ..sequences.hgcalLocalRecoTask_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTBtagDeepFlavourSequencePFPuppiModEta2p4_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTTrackingV61Task_cfi import *
from ..sequences.localrecoTask_cfi import *
from ..sequences.RawToDigiTask_cfi import *

HLT_DoublePFPuppiJets128_DoublePFPuppiBTagDeepFlavour_2p4 = cms.Path(HLTBeginSequence+l1tDoublePFPuppiJet112offMaxEta2p4+l1tDoublePFPuppiJets112offMaxDeta1p6+RawToDigiTask+hgcalLocalRecoTask+localrecoTask+HLTTrackingV61Task+HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+hltDoublePFPuppiJets128MaxEta2p4+hltDoublePFPuppiJets128Eta2p4MaxDeta1p6+HLTBtagDeepFlavourSequencePFPuppiModEta2p4+hltBTagPFPuppiDeepFlavour0p935DoubleEta2p4+HLTEndSequence)
