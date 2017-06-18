import FWCore.ParameterSet.Config as cms

#
# sequence to make photons from clusters in ECAL
#
# photon producer
from RecoEgamma.EgammaPhotonProducers.gedPhotonCore_cfi import *
from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *

import RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi 

gedPhotonsTmp = RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi.gedPhotons.clone()
gedPhotonsTmp.photonProducer = cms.InputTag("gedPhotonCore")
gedPhotonsTmp.candidateP4type = cms.string("fromEcalEnergy")
del gedPhotonsTmp.regressionConfig
gedPhotonsTmp.outputPhotonCollection = cms.string("")
gedPhotonsTmp.reconstructionStep = cms.string("tmp")
gedPhotonSequenceTmp = cms.Sequence(gedPhotonCore+gedPhotonsTmp)


gedPhotons = RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi.gedPhotons.clone()
gedPhotons.photonProducer = cms.InputTag("gedPhotonsTmp")
gedPhotons.outputPhotonCollection = cms.string("")
gedPhotons.reconstructionStep = cms.string("final")
gedPhotons.chargedHadronIsolation = cms.InputTag("egmPhotonIsolationCITK:h+-DR030-")
gedPhotons.neutralHadronIsolation = cms.InputTag("egmPhotonIsolationCITK:h0-DR030-")
gedPhotons.photonIsolation = cms.InputTag("egmPhotonIsolationCITK:gamma-DR030-")

#from Configuration.Eras.Modifier_phase2_hgcal_cff import phase2_hgcal
#phase2_hgcal.toModify(gedPhotons, candidateP4type='fromEcalEnergy')

gedPhotonSequence    = cms.Sequence(gedPhotons)



