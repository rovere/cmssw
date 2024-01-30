import FWCore.ParameterSet.Config as cms

from ..modules.hltPFPuppiHT_cfi import *
from ..modules.hltPFPuppiMHT_cfi import *
from ..sequences.hgcalLocalRecoTask_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTHgcalTiclPFClusteringForEgamma_cfi import *
from ..sequences.HLTJMESequence_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTTrackingV61Task_cfi import *
from ..sequences.localrecoTask_cfi import *
from ..sequences.RawToDigiTask_cfi import *

MC_JME = cms.Path(HLTBeginSequence+RawToDigiTask+hgcalLocalRecoTask+localrecoTask+HLTTrackingV61Task+HLTMuonsSequence+HLTParticleFlowSequence+HLTHgcalTiclPFClusteringForEgamma+HLTJMESequence+hltPFPuppiHT+hltPFPuppiMHT)
