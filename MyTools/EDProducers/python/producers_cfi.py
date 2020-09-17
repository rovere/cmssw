import FWCore.ParameterSet.Config as cms


HoverE = cms.EDProducer(
    "ProducerHoverE",
    
    instanceName = cms.string("HoverE"),
    
    electrons = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
    layerClusters = cms.InputTag("hgcalLayerClusters"),
    
    coneDR = cms.double(0.15),
    
    minClusE = cms.double(0.0),
    minClusET = cms.double(0.0),
    
    debug = cms.bool(False),
)


trackIso = cms.EDProducer(
    "ProducerTrackIso",
    
    instanceName = cms.string("TrackIso"),
    
    electrons = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
    tracks = cms.InputTag("generalTracks"),
    
    isoConeDR = cms.double(0.3),
    vetoConeDR = cms.double(0.01),
    vetoPhiStripDeta = cms.double(0.03),
    
    minTrackPt = cms.double(1.0),
    maxTrackEleDz = cms.double(0.15),
    
    debug = cms.bool(False),
)
