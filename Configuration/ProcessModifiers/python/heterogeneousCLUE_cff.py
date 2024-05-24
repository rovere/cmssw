import FWCore.ParameterSet.Config as cms

# This modifier runs the heterogeneous version of CLUE on the EE part of HGCAL
# in  the HLT Phase2 Simplified menu. All seeded HLT Paths will keep on using the regular CPU version.
heterogeneousCLUE = cms.Modifier()
