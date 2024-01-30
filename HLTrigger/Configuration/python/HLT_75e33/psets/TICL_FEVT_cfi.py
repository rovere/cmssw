import FWCore.ParameterSet.Config as cms

TICL_FEVT = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_ticlSimTracksters_*_*',
        'keep *_ticlSimTICLCandidates_*_*',
        'keep *_ticlSimTrackstersFromCP_*_*',
        'keep *_ticlTrackstersCLUE3DHigh_*_*',
        'keep *_ticlTrackstersMerge_*_*',
        'keep *_ticlTrackstersHFNoseTrkEM_*_*',
        'keep *_ticlTrackstersHFNoseEM_*_*',
        'keep *_ticlTrackstersHFNoseTrk_*_*',
        'keep *_ticlTrackstersHFNoseMIP_*_*',
        'keep *_ticlTrackstersHFNoseHAD_*_*',
        'keep *_ticlTrackstersHFNoseMerge_*_*',
        'keep *_pfTICL_*_*',
        'keep CaloParticles_mix_*_*',
        'keep SimClusters_mix_*_*',
        'keep *_layerClusterSimClusterAssociationProducer_*_*',
        'keep *_layerClusterCaloParticleAssociationProducer_*_*',
        'keep *_layerClusterSimTracksterAssociationProducer_*_*',
        'keep *_tracksterSimTracksterAssociationLinking_*_*',
        'keep *_tracksterSimTracksterAssociationPR_*_*',
        'keep *_tracksterSimTracksterAssociationLinkingPU_*_*',
        'keep *_tracksterSimTracksterAssociationPRPU_*_*',
        'keep *_tracksterSimTracksterAssociationLinkingbyCLUE3D_*_*',
        'keep *_tracksterSimTracksterAssociationPRbyCLUE3D_*_*'
    )
)