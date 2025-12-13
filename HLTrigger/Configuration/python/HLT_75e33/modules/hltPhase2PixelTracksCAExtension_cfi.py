import FWCore.ParameterSet.Config as cms

hltPhase2PixelTracksCAExtension = cms.EDProducer("PixelTrackProducerFromSoAAlpaka",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    minNumberOfHits = cms.int32(0),
    minQuality = cms.string('tight'),
    pixelRecHitLegacySrc = cms.InputTag("hltSiPixelRecHits"),
    trackSrc = cms.InputTag("hltPhase2PixelTracksSoA"),
    outerTrackerRecHitSrc = cms.InputTag("hltSiPhase2RecHits"),
    outerTrackerRecHitSoAConverterSrc = cms.InputTag("hltPhase2OtRecHitsSoA"),
    useOTExtension = cms.bool(True),
    requireQuadsFromConsecutiveLayers = cms.bool(True)
)

_hltPhase2PixelTracksCAFullExtension = cms.EDProducer("PixelTrackProducerFromSoAAlpaka",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    minNumberOfHits = cms.int32(0),
    minQuality = cms.string('tight'),
    pixelRecHitLegacySrc = cms.InputTag("hltSiPixelRecHits"),
    trackSrc = cms.InputTag("hltPhase2PixelTracksSoA"),
    outerTrackerRecHitSrc = cms.InputTag("hltSiPhase2RecHits"),
    outerTrackerRecHitSoAConverterSrc = cms.InputTag("hltPhase2OtRecHitsSoA"),
    useOTExtension = cms.bool(False),
    useOTFullExtension = cms.bool(True),
    requireQuadsFromConsecutiveLayers = cms.bool(True)
)

from Configuration.ProcessModifiers.phase2CAExtensionFull_cff import phase2CAExtensionFull
phase2CAExtensionFull.toReplaceWith(hltPhase2PixelTracksCAExtension, _hltPhase2PixelTracksCAFullExtension)
