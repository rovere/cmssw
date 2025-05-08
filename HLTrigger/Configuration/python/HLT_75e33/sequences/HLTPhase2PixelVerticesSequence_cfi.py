import FWCore.ParameterSet.Config as cms


from ..modules.hltPhase2PixelVertices_cfi import hltPhase2PixelVertices

HLTPhase2PixelVerticesSequence = cms.Sequence(hltPhase2PixelVertices)

from Configuration.ProcessModifiers.alpaka_cff import alpaka
from RecoVertex.PixelVertexFinding.PixelVertexProducerAlpakaPhase2_alpaka import PixelVertexProducerAlpakaPhase2_alpaka as _PixelVertexProducerAlpakaPhase2_alpaka
from RecoVertex.PixelVertexFinding.PixelVertexProducerFromSoAAlpaka import PixelVertexProducerFromSoAAlpaka as _PixelVertexProducerFromSoAAlpaka



hltPhase2PixelVerticesAlpakaSoA = _PixelVertexProducerAlpakaPhase2_alpaka(
        pixelTrackSrc = "hltPhase2PixelTracksSoA",
        PtMin = 1.0)

_hltPhase2PixelVertices = _PixelVertexProducerFromSoAAlpaka(
        src = "hltPhase2PixelVerticesAlpakaSoA",
        beamSpot = "hltOnlineBeamSpot",
        TrackCollection = "hltPhase2PixelTracks")

_HLTPhase2PixelVerticesSequence = cms.Sequence(hltPhase2PixelVerticesAlpakaSoA+hltPhase2PixelVertices)

alpaka.toReplaceWith(hltPhase2PixelVertices, _hltPhase2PixelVertices)
alpaka.toReplaceWith(HLTPhase2PixelVerticesSequence, _HLTPhase2PixelVerticesSequence)
