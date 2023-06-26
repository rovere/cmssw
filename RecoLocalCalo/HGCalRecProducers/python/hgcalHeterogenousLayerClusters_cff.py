import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.hgcalHeterogenousLayerClusters_cfi import hgcalHeterogenousLayerClusters as hgcalHeterogenousLayerClusters_
from RecoLocalCalo.HGCalRecProducers.hgcalMergeLayerClusters_cfi import hgcalMergeLayerClusters as hgcalMergeLayerClusters_

from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import HGCalRecHit

from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import fC_per_ele, HGCAL_noises, hgceeDigitizer, hgchebackDigitizer, hfnoseDigitizer

hgcalHeterogenousLayerClustersEE = hgcalHeterogenousLayerClusters_.clone(
    detector = 'EE',
    recHits = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
    plugin = dict(
        dEdXweights = HGCalRecHit.layerWeights.value(),
        #With the introduction of 7 regional factors (6 for silicon plus 1 for scintillator),
        #we extend fcPerMip (along with noises below) so that it is guaranteed that they have 6 entries.
        fcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value(),
        thicknessCorrection = HGCalRecHit.thicknessCorrection.value(),
        sciThicknessCorrection = HGCalRecHit.sciThicknessCorrection.value(),
        deltasi_index_regemfac = HGCalRecHit.deltasi_index_regemfac.value(),
        fcPerEle = fC_per_ele,
        #Extending noises as fcPerMip, see comment above.
        noises = HGCAL_noises.values.value() + HGCAL_noises.values.value(),
        noiseMip = hgchebackDigitizer.digiCfg.noise.value(),
        type = cms.string('SiCLUE')
    )
)

hgcalHeterogenousLayerClustersHSi = hgcalHeterogenousLayerClusters_.clone(
    detector = 'FH',
    recHits = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    plugin = dict(
        dEdXweights = HGCalRecHit.layerWeights.value(),
        #With the introduction of 7 regional factors (6 for silicon plus 1 for scintillator),
        #we extend fcPerMip (along with noises below) so that it is guaranteed that they have 6 entries.
        fcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value(),
        thicknessCorrection = HGCalRecHit.thicknessCorrection.value(),
        sciThicknessCorrection = HGCalRecHit.sciThicknessCorrection.value(),
        deltasi_index_regemfac = HGCalRecHit.deltasi_index_regemfac.value(),
        fcPerEle = fC_per_ele,
        #Extending noises as fcPerMip, see comment above.
        noises = HGCAL_noises.values.value() + HGCAL_noises.values.value(),
        noiseMip = hgchebackDigitizer.digiCfg.noise.value(),
        type = cms.string('SiCLUE')
    )
)

hgcalHeterogenousMergeLayerClusters = hgcalMergeLayerClusters_.clone(
    layerClustersEE = cms.InputTag('hgcalHeterogenousLayerClustersEE'),
    layerClustersHSi = cms.InputTag('hgcalHeterogenousLayerClustersHSi'),
    layerClustersHSci = cms.InputTag('hgcalLayerClustersHSci'),
    time_layerclustersEE = cms.InputTag('hgcalHeterogenousLayerClustersEE', 'timeLayerCluster'),
    time_layerclustersHSi = cms.InputTag('hgcalHeterogenousLayerClustersHSi', 'timeLayerCluster'),
    time_layerclustersHSci = cms.InputTag('hgcalLayerClustersHSci', 'timeLayerCluster'),
    timeClname = cms.string('timeLayerCluster'),
    mightGet = cms.optional.untracked.vstring
)