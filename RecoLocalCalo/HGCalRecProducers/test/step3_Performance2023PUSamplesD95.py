# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step3 -s RAW2DIGI,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM --conditions auto:phase2_realistic_T21 --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO -n 10 --eventcontent FEVTDEBUGHLT,MINIAODSIM,DQM --geometry Extended2026D88 --era Phase2C17I13M9 --no_exec --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('RERECO',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed1.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed2.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed3.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed4.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed5.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed6.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed7.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed8.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed9.root',
        'file:/data/user/rovere/CMSSW_13_2_0/src/RecoLocalCalo/HGCalRecProducers/test/step3_HGCAL_RecHits_Uncompressed10.root',
        ),
    secondaryFileNames = cms.untracked.vstring(),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Additional output definition

process.hgCalRecHitsSoAProducer = process.hgCalRecHitsSoAProducer_.clone()
process.hgCalLayerClustersSoAProducer = process.hgCalLayerClustersSoAProducer_.clone()
process.hgCalLayerClustersFromAlpakaProducer = process.hgCalLayerClustersFromAlpakaProducer_.clone()
backend_used = 'cuda_async' # 'serial_sync'
process.load("Configuration.StandardSequences.Accelerators_cff")
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import HGCalRecHit
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import fC_per_ele, HGCAL_noises, hgceeDigitizer, hgchebackDigitizer, hfnoseDigitizer
process.hgCalRecHitsSoAProducer.alpaka.backend = backend_used
process.hgCalRecHitsSoAProducer.thicknessCorrection = HGCalRecHit.thicknessCorrection.value()
process.hgCalRecHitsSoAProducer.noises = HGCAL_noises.values.value() + HGCAL_noises.values.value()
process.hgCalRecHitsSoAProducer.dEdXweights = HGCalRecHit.layerWeights.value()
process.hgCalRecHitsSoAProducer.fcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value()
# process.testHGCalRecHitSoA = cms.EDProducer('HGCalRecHitsSoAProducer@alpaka',
#     detector = cms.string('EE'),
#     recHits = cms.InputTag('HGCalRecHit', 'HGCEERecHits'),
#     mightGet = cms.optional.untracked.vstring,
#     alpaka = cms.untracked.PSet(
#       backend = cms.untracked.string('serial_sync')
#       )
#     )
process.hgCalLayerClustersSoAProducer.hgcalRecHitsSoA = 'hgCalRecHitsSoAProducer'
process.hgCalLayerClustersSoAProducer.alpaka.backend = backend_used
process.load('RecoLocalCalo.HGCalRecProducers.hgCalLayerClusterHeterogeneousDumper_cfi')
process.testHGCalSoAs = cms.EndPath(
     process.hgCalRecHitsSoAProducer
    + process.hgCalLayerClustersSoAProducer
#    + process.hgCalLayerClusterHeterogeneousDumper
    + process.hgCalLayerClustersFromAlpakaProducer
    )
process.testHGCalSoAsNoConversionToLegacy = cms.EndPath(
     process.hgCalRecHitsSoAProducer
    + process.hgCalLayerClustersSoAProducer
    )
process.testHGCalSoAsRecHitsOnly = cms.EndPath(
     process.hgCalRecHitsSoAProducer
    )
#
### LEGACY CPU EE CLUSTERS
process.legacy = cms.EndPath(process.hgcalLayerClustersEE)

### CPU EE with external KALOS/CLUE library, full C++, no ALPAKA
process.legacyKalos = cms.EndPath(process.hgcalHeterogenousLayerClustersEE)

## GENERIC CONSUMER TO MEASURE PURE I/O
process.load('FWCore.Modules.genericConsumer_cfi')
process.genericConsumer.eventProducts = ["HGCalRecHit"]
process.pureIO = cms.EndPath(process.genericConsumer)


## THROUGHPUT SERVICE
process.ThroughputService = cms.Service('ThroughputService',
    eventRange = cms.untracked.uint32(1000),
    eventResolution = cms.untracked.uint32(50),
    printEventSummary = cms.untracked.bool(True),
    enableDQM = cms.untracked.bool(False),
#    dqmPathByProcesses = cms.untracked.bool(False),
#    dqmPath = cms.untracked.string('Throughput'),
#    timeRange = cms.untracked.double(1000),
#    timeResolution = cms.untracked.double(1)
    )

process.MessageLogger.cerr.ThroughputService = cms.untracked.PSet(
    limit = cms.untracked.int32(10000000)
    )
# Schedule definition
#process.schedule = cms.Schedule(process.testHGCalSoAsRecHitsOnly)
#process.schedule = cms.Schedule(process.testHGCalSoAsNoConversionToLegacy)
process.schedule = cms.Schedule(process.testHGCalSoAs)
#process.schedule = cms.Schedule(process.legacy)
#process.schedule = cms.Schedule(process.legacyKalos)
#process.schedule = cms.Schedule(process.pureIO)

