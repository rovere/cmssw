import FWCore.ParameterSet.Config as cms
from ..psets.hgcal_reco_constants_cfi import HGCAL_reco_constants as HGCAL_reco_constants

hltHgcalSoARecHitsLayerClustersProducer = cms.EDProducer("HGCalSoARecHitsLayerClustersProducer@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('cuda_async')
    ),
    hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducer"),
    deltac = cms.double(1.3),
    kappa = cms.double(9),
    outlierDeltaFactor = cms.double(2.0)
)
