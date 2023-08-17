import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('RERECO',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.EventContent.EventContent_cff')
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/a7f64f48-4ea9-458d-b192-62c533e20a08.root',
      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/5449dbca-2000-40ba-b5cb-03d27d4e580f.root',
      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/72b85563-e1a0-47de-83d2-78ce8470bc08.root',
      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/d9ffaa46-a7e7-4f38-806a-e7eb5a27caa5.root',
      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/de502344-3b81-4a8e-bcf0-1b95312db64d.root',
      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/d314306f-376e-4ead-b0b0-ddcd686bdbd1.root',
##      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/59c166d4-48ec-4f01-be56-bf068cf0fd12.root',
##      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/74f50b4b-698b-4102-8d6d-9a69cb3d1321.root',
##      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/b1e0988e-bddf-4061-972c-96cdb251a088.root',
##      '/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/f604ecf3-52da-4790-b9e9-834d9dd0abf9.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_0.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_1.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_2.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_3.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_4.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_5.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_6.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_7.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_8.root',
#      'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/wredjeb/PatatrackHackathon2023/CMSSW_13_X/CLUE3DTiming200PUFix/step3/step3_1868855_9.root',
      ),
    secondaryFileNames = cms.untracked.vstring(),
    inputCommands = cms.untracked.vstring("keep *",
                                          "drop *",
                                          "keep *_HGCalRecHit_*_*"),
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    skipEvents = cms.untracked.uint32(900),
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string('GEN-SIM-RECO'),
#        filterName = cms.untracked.string('')
#    ),
    fileName = cms.untracked.string('file:step3_HGCAL_RecHits_Uncompressed10.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
    fastCloning = cms.untracked.bool(False),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
    compressionLevel = cms.untracked.int32(0)
)
process.FEVTDEBUGHLToutput.outputCommands = []
process.FEVTDEBUGHLToutput.outputCommands.extend(['drop *', 'keep *_HGCalRecHit_*_*'])

# Path and EndPath definitions
process.end = cms.EndPath(process.FEVTDEBUGHLToutput)

# # Schedule definition
process.schedule = cms.Schedule(process.end)
#
