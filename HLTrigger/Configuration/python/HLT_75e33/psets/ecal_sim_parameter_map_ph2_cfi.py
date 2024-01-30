import FWCore.ParameterSet.Config as cms

ecal_sim_parameter_map_ph2 = cms.PSet(
    binOfMaximum = cms.int32(6),
    doPhotostatistics = cms.bool(True),
    photoelectronsToAnalogBarrel = cms.double(0.000444444),
    samplingFactor = cms.double(1.0),
    simHitToPhotoelectronsBarrel = cms.double(2250.0),
    syncPhase = cms.bool(True),
    timePhase = cms.double(0.0)
)