import FWCore.ParameterSet.Config as cms

RecoLocalTrackerFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_clusterSummaryProducer_*_*',
        'keep DetIds_siStripDigis_*_*',
        'keep DetIdedmEDCollection_siPixelDigis_*_*',
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*',
        'keep *_siPixelClusters_*_*',
        'keep *_siStripClusters_*_*',
        'keep ClusterSummary_clusterSummaryProducer_*_*',
        'keep *_siPhase2Clusters_*_*',
        'keep *_siPhase2Clusters_*_*'
    )
)