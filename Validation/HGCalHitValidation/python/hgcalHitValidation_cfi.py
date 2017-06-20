import FWCore.ParameterSet.Config as cms

hgcalHitValidation = cms.EDAnalyzer("HGCalHitCalibration",
                                    detector = cms.string("all"),
                                    rawRecHits = cms.bool(True))

hgcalHitValidationSequece = cms.Sequence(hgcalHitValidation)
