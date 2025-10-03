
import FWCore.ParameterSet.Config as cms

from ..modules.hltHgcalDigis_cfi import *
from ..modules.hltHgcalLayerClustersEE_cfi import *
from ..modules.hltHgcalLayerClustersHSci_cfi import *
from ..modules.hltHgcalLayerClustersHSi_cfi import *
from ..modules.hltMergeLayerClusters_cfi import *
from ..modules.hltHGCalRecHit_cfi import *
from ..modules.hltHGCalUncalibRecHit_cfi import *
from ..modules.hltBarrelLayerClustersHB_cfi import *
from ..modules.hltHgcalSoARecHitsProducer_cfi import *
from ..modules.hltHgcalSoARecHitsLayerClustersProducer_cfi import *
from ..modules.hltHgcalSoALayerClustersProducer_cfi import *
from ..modules.hltHgcalLayerClustersFromSoAProducer_cfi import *


HLTHGCalLocalRecoOnlySequence = cms.Sequence((hltHgcalDigis+hltHGCalUncalibRecHit+hltHGCalRecHit+hltHgcalLayerClustersEE+hltHgcalLayerClustersHSci+hltHgcalLayerClustersHSi+hltMergeLayerClusters))

from Configuration.ProcessModifiers.alpaka_cff import alpaka
alpaka.toReplaceWith(HLTHGCalLocalRecoOnlySequence,
                     cms.Sequence(
                                  hltHgcalDigis
                                  + hltHGCalUncalibRecHit
                                  + hltHGCalRecHit
                                  + hltHgcalSoARecHitsProducer
                                  + hltHgcalSoARecHitsLayerClustersProducer
                                  + hltHgcalSoALayerClustersProducer
                                  + hltHgCalLayerClustersFromSoAProducer
                                  + hltHgcalLayerClustersHSci
                                  + hltHgcalLayerClustersHSi
                                  + hltMergeLayerClusters
                     )
                     )
