import FWCore.ParameterSet.Config as cms

RecoMTDAOD = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep intedmValueMap_trackExtenderWithMTD_*_*',
        'keep floatedmValueMap_trackExtenderWithMTD_*_*',
        'keep *_mtdTrackQualityMVA_*_*'
    )
)