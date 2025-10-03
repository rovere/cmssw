import FWCore.ParameterSet.Config as cms

from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTHGCalLocalRecoOnlySequence_cfi import *

MC_HGCalLocalRecoOnly = cms.Path(
    HLTBeginSequence
    + HLTHGCalLocalRecoOnlySequence
)

