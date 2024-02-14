import FWCore.ParameterSet.Config as cms

from ..modules.hltEG187EtL1SeededFilter_cfi import *
from ..modules.hltEgammaCandidatesL1Seeded_cfi import *
from ..modules.hltEgammaCandidatesWrapperL1Seeded_cfi import *
from ..modules.hltEgammaHGCALIDVarsL1Seeded_cfi import *
from ..modules.hltEgammaHoverEL1Seeded_cfi import *
from ..modules.hltPhoton187HEL1SeededFilter_cfi import *
from ..modules.hltPhoton187HgcalHEL1SeededFilter_cfi import *
from ..modules.l1tTkEmSingle51Filter_cfi import *
from ..sequences.HLTDoFullUnpackingEgammaEcalL1SeededSequence_cfi import *
from ..sequences.HLTEGammaDoLocalHcalSequence_cfi import *
from ..sequences.HLTHgcalTiclPFClusteringForEgammaL1SeededSequence_cfi import *
from ..sequences.HLTL1Sequence_cfi import *
from ..sequences.HLTPFClusteringForEgammaL1SeededSequence_cfi import *

HLTPhoton187L1SeededSequence = cms.Sequence(HLTL1Sequence+l1tTkEmSingle51Filter+HLTDoFullUnpackingEgammaEcalL1SeededSequence+HLTEGammaDoLocalHcalSequence+HLTPFClusteringForEgammaL1SeededSequence+HLTHgcalTiclPFClusteringForEgammaL1SeededSequence+(hltEgammaCandidatesL1Seeded+hltEgammaHGCALIDVarsL1Seeded+hltEgammaHoverEL1Seeded)+hltEgammaCandidatesWrapperL1Seeded+hltEG187EtL1SeededFilter+hltPhoton187HgcalHEL1SeededFilter+hltPhoton187HEL1SeededFilter)
