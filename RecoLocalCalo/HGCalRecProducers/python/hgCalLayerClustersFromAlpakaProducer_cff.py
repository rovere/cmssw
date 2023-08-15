import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.hgCalLayerClustersFromAlpakaProducer_cfi import hgCalLayerClustersFromAlpakaProducer as hgCalLayerClustersFromAlpakaProducer_
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import HGCalRecHit

hgCalLayerClustersFromAlpakaProducerEE = hgCalLayerClustersFromAlpakaProducer_.clone(
  hgcalOutSoA = cms.InputTag('hgCalLayerClustersSoAProducerEE'),
  hgcalRecHitsSoA = cms.InputTag('hgCalRecHitsSoAProducerEE'),
  detector = cms.string('EE'),
  recHits = cms.InputTag('HGCalRecHit', 'HGCEERecHits'),
)

hgCalLayerClustersFromAlpakaProducerHSi = hgCalLayerClustersFromAlpakaProducer_.clone(
  hgcalOutSoA = cms.InputTag('hgCalLayerClustersSoAProducerHSi'),
  hgcalRecHitsSoA = cms.InputTag('hgCalRecHitsSoAProducerHSi'),
  detector = cms.string('FH'),
  recHits = cms.InputTag('HGCalRecHit', 'HGCHEFRecHits'),
)