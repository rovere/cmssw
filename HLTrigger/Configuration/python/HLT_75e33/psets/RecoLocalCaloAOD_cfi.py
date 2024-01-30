import FWCore.ParameterSet.Config as cms

RecoLocalCaloAOD = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_castorreco_*_*',
        'keep *_reducedHcalRecHits_*_*',
        'keep HcalUnpackerReport_castorDigis_*_*',
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*',
        'keep HcalUnpackerReport_hcalDigis_*_*',
        'keep *_HGCalRecHit_*_*',
        'keep recoCaloClusters_hgcalMergeLayerClusters_*_*',
        'keep *_hgcalMergeLayerClusters_timeLayerCluster_*',
        'keep *_hgcalMergeLayerClusters_InitialLayerClustersMask_*'
    )
)