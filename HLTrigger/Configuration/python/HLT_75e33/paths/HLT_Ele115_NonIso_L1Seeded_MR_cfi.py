import FWCore.ParameterSet.Config as cms

from ..modules.hltPreEle115NonIsoL1Seeded_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEle115NonIsoL1SeededSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *

HLT_Ele115_NonIso_L1Seeded_MR = cms.Path(HLTBeginSequence+hltPreEle115NonIsoL1Seeded+HLTEle115NonIsoL1SeededSequence+HLTEndSequence)
