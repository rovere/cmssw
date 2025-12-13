import FWCore.ParameterSet.Config as cms

hltPhase2OtRecHitsSoA = cms.EDProducer('Phase2OTRecHitsSoAConverter@alpaka',
  pixelRecHitSoASource = cms.InputTag('hltPhase2SiPixelRecHitsSoA'),
  otRecHitSource = cms.InputTag('hltSiPhase2RecHits'),
  beamSpot = cms.InputTag('hltOnlineBeamSpot'),
  mightGet = cms.optional.untracked.vstring,
  alpaka = cms.untracked.PSet(
    backend = cms.untracked.string('')
  )
)

_hltPhase2OtRecHitsFullSoA = cms.EDProducer('Phase2OTRecHitsFullSoAConverter@alpaka',
  pixelRecHitSoASource = cms.InputTag('hltPhase2SiPixelRecHitsSoA'),
  otRecHitSource = cms.InputTag('hltSiPhase2RecHits'),
  beamSpot = cms.InputTag('hltOnlineBeamSpot'),
  mightGet = cms.optional.untracked.vstring,
  alpaka = cms.untracked.PSet(
    backend = cms.untracked.string('')
  )
)

from Configuration.ProcessModifiers.phase2CAExtensionFull_cff import phase2CAExtensionFull
phase2CAExtensionFull.toReplaceWith(hltPhase2OtRecHitsSoA, _hltPhase2OtRecHitsFullSoA)
