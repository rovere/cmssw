import FWCore.ParameterSet.Config as cms

HFNose_noise_fC = cms.PSet(
    doseMap = cms.string(''),
    scaleByDose = cms.bool(False),
    scaleByDoseAlgo = cms.uint32(0),
    scaleByDoseFactor = cms.double(1),
    values = cms.vdouble(0.336430626, 0.336430626, 0.256328096)
)