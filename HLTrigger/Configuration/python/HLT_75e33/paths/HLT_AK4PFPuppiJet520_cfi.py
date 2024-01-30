import FWCore.ParameterSet.Config as cms

from ..modules.hltSingleAK4PFPuppiJet520_cfi import *
from ..modules.l1tSinglePFPuppiJet230off_cfi import *
from ..sequences.hgcalLocalRecoTask_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTTrackingV61Task_cfi import *
from ..sequences.localrecoTask_cfi import *
from ..sequences.RawToDigiTask_cfi import *

HLT_AK4PFPuppiJet520 = cms.Path(HLTBeginSequence+l1tSinglePFPuppiJet230off+RawToDigiTask+hgcalLocalRecoTask+localrecoTask+HLTTrackingV61Task+HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+hltSingleAK4PFPuppiJet520+HLTEndSequence)
