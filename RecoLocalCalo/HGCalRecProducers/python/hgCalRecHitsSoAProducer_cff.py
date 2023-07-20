import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.hgCalRecHitsSoAProducer_cfi import hgCalRecHitsSoAProducer as hgCalRecHitsSoAProducer_
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import HGCalRecHit

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import fC_per_ele, HGCAL_noises, hgceeDigitizer, hgchebackDigitizer, hfnoseDigitizer


hgCalRecHitsSoAProducerEE = hgCalRecHitsSoAProducer_.clone(
  detector = cms.string('EE'),
  recHits = cms.InputTag('HGCalRecHit', 'HGCEERecHits'),
  fcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value(),
  thicknessCorrection = HGCalRecHit.thicknessCorrection.value(),
  noises = HGCAL_noises.values.value() + HGCAL_noises.values.value(),
  dEdXweights = HGCalRecHit.layerWeights.value(),
  alpaka = cms.untracked.PSet(
    backend = cms.untracked.string('cuda_async')
  )
)

hgCalRecHitsSoAProducerHSi = hgCalRecHitsSoAProducer_.clone(
  detector = cms.string('FH'),
  recHits = cms.InputTag('HGCalRecHit', 'HGCHEFRecHits'),
  fcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value(),
  thicknessCorrection = HGCalRecHit.thicknessCorrection.value(),
  noises = HGCAL_noises.values.value() + HGCAL_noises.values.value(),
  dEdXweights = HGCalRecHit.layerWeights.value(),
  alpaka = cms.untracked.PSet(
    backend = cms.untracked.string('cuda_async')
  )
)