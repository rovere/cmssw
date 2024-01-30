import FWCore.ParameterSet.Config as cms

RecoMuonIsolationRECO = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_muIsoDepositTk_*_*',
        'keep *_muIsoDepositCalByAssociatorTowers_*_*',
        'keep *_muIsoDepositCalByAssociatorHits_*_*',
        'keep *_muIsoDepositJets_*_*',
        'keep *_muIsoDepositTkDisplaced_*_*',
        'keep *_muIsoDepositCalByAssociatorTowersDisplaced_*_*',
        'keep *_muIsoDepositCalByAssociatorHitsDisplaced_*_*',
        'keep *_muIsoDepositJetsDisplaced_*_*',
        'keep *_muGlobalIsoDepositCtfTk_*_*',
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*',
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*',
        'keep *_muGlobalIsoDepositJets_*_*'
    )
)