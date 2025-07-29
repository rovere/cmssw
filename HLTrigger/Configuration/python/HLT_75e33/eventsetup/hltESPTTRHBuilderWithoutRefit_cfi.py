import FWCore.ParameterSet.Config as cms

def _addProcessTTRHBuilderWithoutRefit(process):
    process.hltESPTTRHBuilderWithoutRefit = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
        ComponentName = cms.string('hltESPTTRHBuilderWithoutRefit'),
        ComputeCoarseLocalPositionFromDisk = cms.bool(False),
        Matcher = cms.string('Fake'),
        Phase2StripCPE = cms.string(''),
        PixelCPE = cms.string('Fake'),
        StripCPE = cms.string('Fake')
    )

from Configuration.ProcessModifiers.trackingLST_cff import trackingLST
from Configuration.ProcessModifiers.ngtScouting_cff import ngtScouting
modifyConfigurationForTrackingLSTTTRHBuilderWithoutRefit_ = (trackingLST | ngtScouting).makeProcessModifier(_addProcessTTRHBuilderWithoutRefit)
