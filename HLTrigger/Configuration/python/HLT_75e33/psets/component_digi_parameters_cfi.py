import FWCore.ParameterSet.Config as cms

component_digi_parameters = cms.PSet(
    componentAddToBarrel = cms.bool(False),
    componentDigiTag = cms.string('Component'),
    componentSeparateDigi = cms.bool(False),
    componentTimePhase = cms.double(0.0),
    componentTimeTag = cms.string('Component')
)