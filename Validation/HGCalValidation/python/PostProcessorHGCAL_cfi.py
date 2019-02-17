import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

eff_layers = ["effic_layer%d 'LayerCluster Efficiency vs #eta Layer%d' Num_CaloParticle_Eta_perlayer%d Denom_CaloParticle_Eta_perlayer%d" % (i, i, i, i)  for i in range(1,53) ]

postProcessorHGCAL = DQMEDHarvester('DQMGenericClient',
    subDirs = cms.untracked.vstring('HGCAL/HGCalValidator/hgcalLayerClusters/'),
    efficiency = cms.vstring(eff_layers),
    resolution = cms.vstring(),
    cumulativeDists = cms.untracked.vstring(),
    noFlowDists = cms.untracked.vstring(),
    outputFileName = cms.untracked.string(""),
    verbose = cms.untracked.uint32(4))

