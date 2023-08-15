import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import *
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import *

from RecoLocalCalo.HGCalRecProducers.hgcalRecHitMapProducer_cfi import hgcalRecHitMapProducer

# patch particle flow clusters for HGC into local reco sequence
# (for now until global reco is going with some sort of clustering)
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHGC_cfi import *
from RecoLocalCalo.HGCalRecProducers.hgcalMultiClusters_cfi import *
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import hgcalLayerClustersHFNose, hgcalLayerClustersEE, hgcalLayerClustersHSi, hgcalLayerClustersHSci
from RecoLocalCalo.HGCalRecProducers.hgcalMergeLayerClusters_cff import hgcalMergeLayerClusters
from RecoLocalCalo.HGCalRecProducers.hgcalHeterogenousLayerClusters_cff import hgcalHeterogenousLayerClustersEE, hgcalHeterogenousLayerClustersHSi
from RecoLocalCalo.HGCalRecProducers.hgCalRecHitsSoAProducer_cff import *
from RecoLocalCalo.HGCalRecProducers.hgCalLayerClustersSoAProducer_cff import *
from RecoLocalCalo.HGCalRecProducers.hgCalLayerClustersFromAlpakaProducer_cff import *
from Configuration.StandardSequences.Accelerators_cff import *
from HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi import *

hgcalLocalRecoTask = cms.Task( HGCalUncalibRecHit,
                                       HGCalRecHit,
                                       hgcalRecHitMapProducer,
                                       hgcalLayerClustersEE,
                                    #    hgCalRecHitsSoAProducerEE,
                                    #    hgCalLayerClustersSoAProducerEE,
                                    #    hgCalLayerClustersFromAlpakaProducerEE,
				                       hgcalLayerClustersHSi,
                                    #    hgCalRecHitsSoAProducerHSi,
                                    #    hgCalLayerClustersSoAProducerHSi,
                                    #    hgCalLayerClustersFromAlpakaProducerHSi,
                                       hgcalLayerClustersHSci,
				                       hgcalMergeLayerClusters,
                                       hgcalMultiClusters,
                                       particleFlowRecHitHGC,
                                       particleFlowClusterHGCal )

_hfnose_hgcalLocalRecoTask = hgcalLocalRecoTask.copy()
_hfnose_hgcalLocalRecoTask.add(hgcalLayerClustersHFNose)

from Configuration.Eras.Modifier_phase2_hfnose_cff import phase2_hfnose
phase2_hfnose.toReplaceWith(
    hgcalLocalRecoTask, _hfnose_hgcalLocalRecoTask )

hgcalLocalRecoSequence = cms.Sequence(hgcalLocalRecoTask)
