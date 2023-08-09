import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.hgcalMergeLayerClusters_cfi import hgcalMergeLayerClusters as hgcalMergeLayerClusters_

hgcalMergeLayerClusters = hgcalMergeLayerClusters_.clone(
)
# hgcalMergeLayerClusters = hgcalMergeLayerClusters_.clone(
#     layerClustersEE = cms.InputTag('hgCalLayerClustersFromAlpakaProducerEE'),
#     layerClustersHSi = cms.InputTag('hgCalLayerClustersFromAlpakaProducerHSi'),
#     layerClustersHSci = cms.InputTag('hgcalLayerClustersHSci'),
#     time_layerclustersEE = cms.InputTag('hgCalLayerClustersFromAlpakaProducerEE', 'timeLayerCluster'),
#     time_layerclustersHSi = cms.InputTag('hgCalLayerClustersFromAlpakaProducerHSi', 'timeLayerCluster'),
#     time_layerclustersHSci = cms.InputTag('hgcalLayerClustersHSci', 'timeLayerCluster'),
#     timeClname = cms.string('timeLayerCluster'),
#     mightGet = cms.optional.untracked.vstring
# )