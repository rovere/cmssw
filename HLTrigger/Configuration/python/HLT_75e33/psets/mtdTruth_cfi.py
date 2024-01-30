import FWCore.ParameterSet.Config as cms

mtdTruth = cms.PSet(
    HepMCProductLabel = cms.InputTag("generatorSmeared"),
    MaxPseudoRapidity = cms.double(5.0),
    MinEnergy = cms.double(0.5),
    accumulatorType = cms.string('MtdTruthAccumulator'),
    allowDifferentSimHitProcesses = cms.bool(False),
    bunchspace = cms.uint32(25),
    genParticleCollection = cms.InputTag("genParticles"),
    maximumPreviousBunchCrossing = cms.uint32(0),
    maximumSubsequentBunchCrossing = cms.uint32(0),
    premixStage1 = cms.bool(False),
    simHitCollections = cms.PSet(
        mtdCollections = cms.VInputTag(cms.InputTag("g4SimHits","FastTimerHitsBarrel"), cms.InputTag("g4SimHits","FastTimerHitsEndcap"))
    ),
    simTrackCollection = cms.InputTag("g4SimHits"),
    simVertexCollection = cms.InputTag("g4SimHits")
)