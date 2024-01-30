import FWCore.ParameterSet.Config as cms

RecoMTDRECO = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep recoTrack*_trackExtenderWithMTD_*_*',
        'keep TrackingRecHitsOwned_trackExtenderWithMTD_*_*',
        'keep intedmValueMap_trackExtenderWithMTD_*_*',
        'keep floatedmValueMap_trackExtenderWithMTD_*_*',
        'keep *_mtdTrackQualityMVA_*_*'
    )
)