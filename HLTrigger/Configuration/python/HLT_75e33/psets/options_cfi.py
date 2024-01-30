import FWCore.ParameterSet.Config as cms

options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(
        'IntermediateHitDoublets_highPtTripletStepHitDoublets__HLTX',
        'IntermediateHitDoublets_hltElePixelHitDoubletsForTripletsL1Seeded__HLTX',
        'IntermediateHitDoublets_hltElePixelHitDoubletsForTripletsUnseeded__HLTX',
        'IntermediateHitDoublets_hltElePixelHitDoubletsL1Seeded__HLTX',
        'IntermediateHitDoublets_hltElePixelHitDoubletsUnseeded__HLTX',
        'IntermediateHitDoublets_hltIter2Phase2L3FromL1TkMuonPixelHitDoublets__HLTX',
        'IntermediateHitDoublets_hltPhase2L3FromL1TkMuonPixelTracksHitDoublets__HLTX',
        'IntermediateHitDoublets_hltPhase2L3MuonHighPtTripletStepHitDoublets__HLTX',
        'IntermediateHitDoublets_hltPhase2L3MuonPixelTracksHitDoublets__HLTX',
        'IntermediateHitDoublets_hltPhase2PixelTracksHitDoublets__HLTX',
        'RegionsSeedingHitSets_highPtTripletStepHitTriplets__HLTX',
        'RegionsSeedingHitSets_hltElePixelHitDoubletsForTripletsL1Seeded__HLTX',
        'RegionsSeedingHitSets_hltElePixelHitDoubletsForTripletsUnseeded__HLTX',
        'RegionsSeedingHitSets_hltElePixelHitDoubletsL1Seeded__HLTX',
        'RegionsSeedingHitSets_hltElePixelHitDoubletsUnseeded__HLTX',
        'RegionsSeedingHitSets_hltElePixelHitTripletsL1Seeded__HLTX',
        'RegionsSeedingHitSets_hltElePixelHitTripletsUnseeded__HLTX',
        'RegionsSeedingHitSets_hltIter2Phase2L3FromL1TkMuonPixelHitTriplets__HLTX',
        'RegionsSeedingHitSets_hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets__HLTX',
        'RegionsSeedingHitSets_hltPhase2L3MuonHighPtTripletStepHitTriplets__HLTX',
        'RegionsSeedingHitSets_hltPhase2L3MuonPixelTracksHitQuadruplets__HLTX',
        'RegionsSeedingHitSets_hltPhase2PixelTracksHitSeeds__HLTX'
    ),
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
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)