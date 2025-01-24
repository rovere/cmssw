import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
from Configuration.Eras.Modifier_fastSim_cff import fastSim

def _addNoFlow(module):
    _noflowSeen = set()
    for eff in module.efficiency.value():
        tmp = eff.split(" ")
        if "cut" in tmp[0]:
            continue
        ind = -1
        if tmp[ind] == "fake" or tmp[ind] == "simpleratio":
            ind = -2
        if not tmp[ind] in _noflowSeen:
            module.noFlowDists.append(tmp[ind])
        if not tmp[ind-1] in _noflowSeen:
            module.noFlowDists.append(tmp[ind-1])

_defaultSubdirs = ["Tracking/TrackingMCTruth/SimDoublets"]

postProcessorSimDoublets = DQMEDHarvester("DQMGenericClient",
    subDirs = cms.untracked.vstring(_defaultSubdirs),
    efficiency = cms.vstring(
        "efficiency_vs_pT 'Efficiency vs p_{T}; True transverse momentum p_{T} [GeV]; Efficiency for valid SimDoublets' numPassVsPt numTotVsPt",
        "efficiency_vs_eta 'Efficiency vs #eta; True pseudorapidity #eta; Efficiency for valid SimDoublets' numPassVsEta numTotVsEta"
    ),
    resolution = cms.vstring(),
    cumulativeDists = cms.untracked.vstring(),
    noFlowDists = cms.untracked.vstring(),
    outputFileName = cms.untracked.string("")
)

_addNoFlow(postProcessorSimDoublets)
