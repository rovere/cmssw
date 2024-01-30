import FWCore.ParameterSet.Config as cms

SimTrackerRAW = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_allTrackMCMatch_*_*',
        'keep *_prunedTrackingParticles_*_*',
        'keep *_prunedDigiSimLinks_*_*'
    )
)