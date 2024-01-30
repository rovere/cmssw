import FWCore.ParameterSet.Config as cms

from ..modules.hltAK4PFJetsForTaus_cfi import *
from ..modules.hltHpsDoublePFTau40TrackPt1MediumChargedIsolation_cfi import *
from ..modules.hltHpsSelectedPFTausTrackPt1MediumChargedIsolation_cfi import *
from ..modules.hltPreDoublePFTauHPS_cfi import *
from ..sequences.hgcalLocalRecoTask_cfi import *
from ..sequences.HLTAK4PFJetsReconstruction_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTHPSMediumChargedIsoPFTauSequence_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTPFTauHPS_cfi import *
from ..sequences.HLTTrackingV61Task_cfi import *
from ..sequences.localrecoTask_cfi import *
from ..sequences.RawToDigiTask_cfi import *

HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1 = cms.Path(HLTBeginSequence+hltPreDoublePFTauHPS+RawToDigiTask+hgcalLocalRecoTask+localrecoTask+HLTTrackingV61Task+HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFJetsReconstruction+hltAK4PFJetsForTaus+HLTPFTauHPS+HLTHPSMediumChargedIsoPFTauSequence+hltHpsSelectedPFTausTrackPt1MediumChargedIsolation+hltHpsDoublePFTau40TrackPt1MediumChargedIsolation+HLTEndSequence)
