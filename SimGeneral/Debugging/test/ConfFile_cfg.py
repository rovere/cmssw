import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("patatrackAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1)
                                        )
process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
                                'file:/data/agruber/patatrack/CMSSW_13_2_0_pre3/src/24036.0_ZTT_14TeV+2026D96/step3.root'
                                #'file:/eos/user/a/agruber/tauRecoTICL/patatrack14/CMSSW_13_2_0_pre3/src/24088.0_SinglePiPt25Eta1p7_2p7+2026D96/step3.root'
                                #'file:/eos/user/a/agruber/samples/TICLDumperTest/new/zpos/n10000/Eta_29/singletau_flatEGun_hgcalCenter/step3/step3clue3D_singletau_e100GeV_eta29_zpos_events10000_nopu.root'
                                )                
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('patatrack_testing.root')
                                   )

process.patatrackAnalyzer = cms.EDAnalyzer('patatrackAnalyzer',
CaloParticle      = cms.InputTag('mix', 'MergedCaloTruth'),
SimClusters       = cms.InputTag('mix', 'MergedCaloTruth'),
GenParticles      = cms.InputTag('genParticles'),
genParticles      = cms.InputTag('genParticles')
                              )

process.p = cms.Path(process.patatrackAnalyzer)
