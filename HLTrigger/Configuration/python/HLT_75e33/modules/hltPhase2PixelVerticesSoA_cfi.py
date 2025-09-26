import FWCore.ParameterSet.Config as cms

from RecoVertex.PixelVertexFinding.PixelVertexProducerAlpakaPhase2_alpaka import PixelVertexProducerAlpakaPhase2_alpaka as _PixelVertexProducerAlpakaPhase2_alpaka

hltPhase2PixelVerticesSoA = _PixelVertexProducerAlpakaPhase2_alpaka(
    pixelTrackSrc = "hltPhase2PixelTracksSoA",
    errmax = 0.015,
    chi2max = cms.double(9.0),
    maxVertices = 512,
    doSplitting = cms.bool(True),
    oneKernel = cms.bool(False),
    useDBSCAN = cms.bool(False),
    useDensity = cms.bool(False),
    useDensityClue = cms.bool(True),
    useIterative = cms.bool(False),
    maxChi2ForFirstFit = cms.double(50.),
    maxChi2ForFinalFit = cms.double(5000.),
    maxChi2ForSplit = cms.double(4.)
)
