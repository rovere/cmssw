import FWCore.ParameterSet.Config as cms

from ..sequences.hgcalLocalRecoTask_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTBtagDeepCSVSequencePFPuppi_cfi import *
from ..sequences.HLTBtagDeepFlavourSequencePFPuppi_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTTrackingV61Task_cfi import *
from ..sequences.localrecoTask_cfi import *
from ..sequences.RawToDigiTask_cfi import *

MC_BTV = cms.Path(HLTBeginSequence+RawToDigiTask+hgcalLocalRecoTask+localrecoTask+HLTTrackingV61Task+HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+HLTBtagDeepCSVSequencePFPuppi+HLTBtagDeepFlavourSequencePFPuppi)
