import FWCore.ParameterSet.Config as cms

def _addProcessModulesDevLST(process):
    process.hltESPModulesDevLST = cms.ESProducer('LSTModulesDevESProducer@alpaka',
        appendToDataLabel = cms.string(''),
        alpaka = cms.untracked.PSet(
            backend = cms.untracked.string('')
        )
    )

from Configuration.ProcessModifiers.trackingLST_cff import trackingLST
from Configuration.ProcessModifiers.ngtScouting_cff import ngtScouting
modifyConfigurationForTrackingLSTModulesDevLST_ = (trackingLST | ngtScouting).makeProcessModifier(_addProcessModulesDevLST)
