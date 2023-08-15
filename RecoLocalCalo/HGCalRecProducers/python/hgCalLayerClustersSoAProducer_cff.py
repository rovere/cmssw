import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.hgCalLayerClustersSoAProducer_cfi import hgCalLayerClustersSoAProducer as hgCalLayerClustersSoAProducer_
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import HGCalRecHit

hgCalLayerClustersSoAProducerEE = hgCalLayerClustersSoAProducer_.clone(
  hgcalRecHitsSoA = cms.InputTag('hgCalRecHitsSoAProducerEE'),
  alpaka = cms.untracked.PSet(
    backend = cms.untracked.string('cuda_async')
    # backend = cms.untracked.string('serial_sync')
  )
)

hgCalLayerClustersSoAProducerHSi = hgCalLayerClustersSoAProducer_.clone(
  hgcalRecHitsSoA = cms.InputTag('hgCalRecHitsSoAProducerHSi'),
  alpaka = cms.untracked.PSet(
    backend = cms.untracked.string('cuda_async')
    # backend = cms.untracked.string('serial_sync') 
  )
)