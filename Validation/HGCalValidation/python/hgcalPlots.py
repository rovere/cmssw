import os
import copy
import collections

import six
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from Validation.RecoTrack.plotting.plotting import Plot, PlotGroup, PlotFolder, Plotter, PlotOnSideGroup
from Validation.RecoTrack.plotting.html import PlotPurpose
import Validation.RecoTrack.plotting.plotting as plotting
import Validation.RecoTrack.plotting.validation as validation
import Validation.RecoTrack.plotting.html as html

#Could i get this from somewhere?
lastLayerEE = 28  # last layer of EE
lastLayerFH = 40  # last layer of FH
maxlayer = 52 # last layer of BH

_common = {"stat": True, "drawStyle": "hist"}
_legend_common = {"legendDx": -0.3,
                  "legendDy": -0.05,
                  "legendDw": 0.1}

_SelectedCaloParticles = PlotGroup("SelectedCaloParticles", [
        Plot("num_caloparticle_eta", xtitle="", **_common),
        Plot("caloparticle_energy", xtitle="", **_common),
        Plot("caloparticle_pt", xtitle="", **_common),
        Plot("caloparticle_phi", xtitle="", **_common),
        Plot("Eta vs Zorigin", xtitle="", **_common),
       ])

_num_reco_cluster_eta = PlotGroup("num_reco_cluster_eta", [
        Plot("num_reco_cluster_eta", xtitle="", **_common),
        ],ncols=1)

_mixedhitscluster = PlotGroup("mixedhitscluster", [
        Plot("mixedhitscluster", xtitle="", **_common),
        ],ncols=1)

_energyclustered = PlotGroup("energyclustered", [
        Plot("energyclustered", xtitle="", **_common),
        ],ncols=1)

_longdepthbarycentre = PlotGroup("longdepthbarycentre", [
        Plot("longdepthbarycentre", xtitle="", **_common),
        ],ncols=1)

_common_layerperthickness = {}
_common_layerperthickness.update(_common)
_common_layerperthickness['xmin'] = 0.
_common_layerperthickness['xmax'] = 100
_totclusternum_thick = PlotGroup("totclusternum_thick", [
    Plot("totclusternum_thick_120", xtitle="", **_common_layerperthickness),
    Plot("totclusternum_thick_200", xtitle="", **_common_layerperthickness),
    Plot("totclusternum_thick_300", xtitle="", **_common_layerperthickness),
    Plot("totclusternum_thick_-1", xtitle="", **_common_layerperthickness),
    Plot("mixedhitscluster", xtitle="", **_common_layerperthickness),
    ])

_totclusternum_layer_EE = PlotGroup("totclusternum_layer_EE", [
        Plot("totclusternum_layer_{:02d}".format(i+1), xtitle="", **_common) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_totclusternum_layer_FH = PlotGroup("totclusternum_layer_FH", [
        Plot("totclusternum_layer_{:02d}".format(i+1), xtitle="", **_common) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_totclusternum_layer_BH = PlotGroup("totclusternum_layer_BH", [
        Plot("totclusternum_layer_{:02d}".format(i+1), xtitle="", **_common) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

_energyclustered_perlayer_EE = PlotGroup("energyclustered_perlayer_EE", [
        Plot("energyclustered_perlayer{:02d}".format(i+1), xtitle="", **_common) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_energyclustered_perlayer_FH = PlotGroup("energyclustered_perlayer_FH", [
        Plot("energyclustered_perlayer{:02d}".format(i+1), xtitle="", **_common) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_energyclustered_perlayer_BH = PlotGroup("energyclustered_perlayer_BH", [
        Plot("energyclustered_perlayer{:02d}".format(i+1), xtitle="", **_common) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )
_cellsenedens_thick =  PlotGroup("cellsenedens_thick", [
    Plot("cellsenedens_thick_120", xtitle="", **_common),
    Plot("cellsenedens_thick_200", xtitle="", **_common),
    Plot("cellsenedens_thick_300", xtitle="", **_common),
    Plot("cellsenedens_thick_-1", xtitle="", **_common),
    ])

#----------------------------------------------------------------------------------------------------------------
#120 um
_common_cells = {}
_common_cells.update(_common)
_common_cells["xmin"] = 0
_common_cells["xmax"] = 50
_common_cells["ymin"] = 0.1
_common_cells["ymax"] = 10000
_common_cells["ylog"] = True
_cellsnum_perthick_perlayer_120_EE = PlotGroup("cellsnum_perthick_perlayer_120_EE", [
  Plot("cellsnum_perthick_perlayer_120_{:02}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_cellsnum_perthick_perlayer_120_FH = PlotGroup("cellsnum_perthick_perlayer_120_FH", [
        Plot("cellsnum_perthick_perlayer_120_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_cellsnum_perthick_perlayer_120_BH = PlotGroup("cellsnum_perthick_perlayer_120_BH", [
        Plot("cellsnum_perthick_perlayer_120_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#200 um
_cellsnum_perthick_perlayer_200_EE = PlotGroup("cellsnum_perthick_perlayer_200_EE", [
        Plot("cellsnum_perthick_perlayer_200_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_cellsnum_perthick_perlayer_200_FH = PlotGroup("cellsnum_perthick_perlayer_200_FH", [
        Plot("cellsnum_perthick_perlayer_200_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_cellsnum_perthick_perlayer_200_BH = PlotGroup("cellsnum_perthick_perlayer_200_BH", [
        Plot("cellsnum_perthick_perlayer_200_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#300 um
_cellsnum_perthick_perlayer_300_EE = PlotGroup("cellsnum_perthick_perlayer_300_EE", [
        Plot("cellsnum_perthick_perlayer_300_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_cellsnum_perthick_perlayer_300_FH = PlotGroup("cellsnum_perthick_perlayer_300_FH", [
        Plot("cellsnum_perthick_perlayer_300_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_cellsnum_perthick_perlayer_300_BH = PlotGroup("cellsnum_perthick_perlayer_300_BH", [
        Plot("cellsnum_perthick_perlayer_300_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#scint um
_cellsnum_perthick_perlayer_scint_EE = PlotGroup("cellsnum_perthick_perlayer_-1_EE", [
        Plot("cellsnum_perthick_perlayer_-1_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_cellsnum_perthick_perlayer_scint_FH = PlotGroup("cellsnum_perthick_perlayer_-1_FH", [
        Plot("cellsnum_perthick_perlayer_-1_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_cellsnum_perthick_perlayer_scint_BH = PlotGroup("cellsnum_perthick_perlayer_-1_BH", [
        Plot("cellsnum_perthick_perlayer_-1_{:02d}".format(i+1), xtitle="", **_common_cells) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#----------------------------------------------------------------------------------------------------------------
#120 um
_common_distance = {}
_common_distance.update(_common)
_common_distance.update(_legend_common)
_common_distance["xmax"] = 150
_common_distance["stat"] = False
_common_distance["ymin"] = 1e-3
_common_distance["ymax"] = 10000
_common_distance["ylog"] = True


_distancetomaxcell_perthickperlayer_120_EE = PlotGroup("distancetomaxcell_perthickperlayer_120_EE", [
        Plot("distancetomaxcell_perthickperlayer_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_120_FH = PlotGroup("distancetomaxcell_perthickperlayer_120_FH", [
        Plot("distancetomaxcell_perthickperlayer_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_120_BH = PlotGroup("distancetomaxcell_perthickperlayer_120_BH", [
        Plot("distancetomaxcell_perthickperlayer_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#200 um
_distancetomaxcell_perthickperlayer_200_EE = PlotGroup("distancetomaxcell_perthickperlayer_200_EE", [
        Plot("distancetomaxcell_perthickperlayer_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_200_FH = PlotGroup("distancetomaxcell_perthickperlayer_200_FH", [
        Plot("distancetomaxcell_perthickperlayer_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_200_BH = PlotGroup("distancetomaxcell_perthickperlayer_200_BH", [
        Plot("distancetomaxcell_perthickperlayer_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#300 um
_distancetomaxcell_perthickperlayer_300_EE = PlotGroup("distancetomaxcell_perthickperlayer_300_EE", [
        Plot("distancetomaxcell_perthickperlayer_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_300_FH = PlotGroup("distancetomaxcell_perthickperlayer_300_FH", [
        Plot("distancetomaxcell_perthickperlayer_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_300_BH = PlotGroup("distancetomaxcell_perthickperlayer_300_BH", [
        Plot("distancetomaxcell_perthickperlayer_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#scint um
_distancetomaxcell_perthickperlayer_scint_EE = PlotGroup("distancetomaxcell_perthickperlayer_-1_EE", [
        Plot("distancetomaxcell_perthickperlayer_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_scint_FH = PlotGroup("distancetomaxcell_perthickperlayer_-1_FH", [
        Plot("distancetomaxcell_perthickperlayer_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_scint_BH = PlotGroup("distancetomaxcell_perthickperlayer_-1_BH", [
        Plot("distancetomaxcell_perthickperlayer_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )


#----------------------------------------------------------------------------------------------------------------
#120 um
_distancetoseedcell_perthickperlayer_120_EE = PlotGroup("distancetoseedcell_perthickperlayer_120_EE", [
        Plot("distancetoseedcell_perthickperlayer_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_120_FH = PlotGroup("distancetoseedcell_perthickperlayer_120_FH", [
        Plot("distancetoseedcell_perthickperlayer_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_120_BH = PlotGroup("distancetoseedcell_perthickperlayer_120_BH", [
        Plot("distancetoseedcell_perthickperlayer_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#200 um
_distancetoseedcell_perthickperlayer_200_EE = PlotGroup("distancetoseedcell_perthickperlayer_200_EE", [
        Plot("distancetoseedcell_perthickperlayer_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_200_FH = PlotGroup("distancetoseedcell_perthickperlayer_200_FH", [
        Plot("distancetoseedcell_perthickperlayer_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_200_BH = PlotGroup("distancetoseedcell_perthickperlayer_200_BH", [
        Plot("distancetoseedcell_perthickperlayer_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#300 um
_distancetoseedcell_perthickperlayer_300_EE = PlotGroup("distancetoseedcell_perthickperlayer_300_EE", [
        Plot("distancetoseedcell_perthickperlayer_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_300_FH = PlotGroup("distancetoseedcell_perthickperlayer_300_FH", [
        Plot("distancetoseedcell_perthickperlayer_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_300_BH = PlotGroup("distancetoseedcell_perthickperlayer_300_BH", [
        Plot("distancetoseedcell_perthickperlayer_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#scint um
_distancetoseedcell_perthickperlayer_scint_EE = PlotGroup("distancetoseedcell_perthickperlayer_-1_EE", [
        Plot("distancetoseedcell_perthickperlayer_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_scint_FH = PlotGroup("distancetoseedcell_perthickperlayer_-1_FH", [
        Plot("distancetoseedcell_perthickperlayer_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_scint_BH = PlotGroup("distancetoseedcell_perthickperlayer_-1_BH", [
        Plot("distancetoseedcell_perthickperlayer_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#=====================================================================================================================
#----------------------------------------------------------------------------------------------------------------
#120 um
_distancetomaxcell_perthickperlayer_eneweighted_120_EE = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_120_EE", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_eneweighted_120_FH = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_120_FH", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_eneweighted_120_BH = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_120_BH", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#200 um
_distancetomaxcell_perthickperlayer_eneweighted_200_EE = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_200_EE", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_eneweighted_200_FH = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_200_FH", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_eneweighted_200_BH = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_200_BH", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#300 um
_distancetomaxcell_perthickperlayer_eneweighted_300_EE = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_300_EE", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_eneweighted_300_FH = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_300_FH", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_eneweighted_300_BH = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_300_BH", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#scint um
_distancetomaxcell_perthickperlayer_eneweighted_scint_EE = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_-1_EE", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_eneweighted_scint_FH = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_-1_FH", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetomaxcell_perthickperlayer_eneweighted_scint_BH = PlotGroup("distancetomaxcell_perthickperlayer_eneweighted_-1_BH", [
        Plot("distancetomaxcell_perthickperlayer_eneweighted_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )


#----------------------------------------------------------------------------------------------------------------
#120 um
_distancetoseedcell_perthickperlayer_eneweighted_120_EE = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_120_EE", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_eneweighted_120_FH = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_120_FH", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_eneweighted_120_BH = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_120_BH", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_120_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#200 um
_distancetoseedcell_perthickperlayer_eneweighted_200_EE = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_200_EE", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_eneweighted_200_FH = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_200_FH", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_eneweighted_200_BH = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_200_BH", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_200_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#300 um
_distancetoseedcell_perthickperlayer_eneweighted_300_EE = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_300_EE", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_eneweighted_300_FH = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_300_FH", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_eneweighted_300_BH = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_300_BH", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_300_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )

#scint um
_distancetoseedcell_perthickperlayer_eneweighted_scint_EE = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_-1_EE", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_eneweighted_scint_FH = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_-1_FH", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerEE,lastLayerFH)
        ],
                                    ncols=4
                                    )

_distancetoseedcell_perthickperlayer_eneweighted_scint_BH = PlotGroup("distancetoseedcell_perthickperlayer_eneweighted_-1_BH", [
        Plot("distancetoseedcell_perthickperlayer_eneweighted_-1_{:02d}".format(i+1), xtitle="", **_common_distance) for i in range(lastLayerFH,maxlayer)
        ],
                                    ncols=4
                                    )


_common_score = {"title": "Score CaloParticle to LayerClusters",
                 "stat": False,
                 "ymin": 0.1,
                 "ymax": 1000,
                 "xmin": 0,
                 "xmax": 1,
                 "drawStyle": "hist",
                 "lineWidth": 1,
                 "ylog": True
                }
_common_score.update(_legend_common)
_score_caloparticle_to_layerclusters = PlotGroup("score_caloparticle_to_layercluster", [
        Plot("Score_caloparticle2layercl_perlayer{:02d}".format(i+1), xtitle="Layer {:02d}".format(i+1), **_common_score) for i in range(0,maxlayer)
        ],
                                    ncols=8
                                    )

_common_score = {"title": "Score LayerCluster to CaloParticles",
                 "stat": False,
                 "ymin": 0.1,
                 "ymax": 1000,
                 "xmin": 0,
                 "xmax": 1,
                 "drawStyle": "hist",
                 "lineWidth": 1,
                 "ylog": True
                }
_common_score.update(_legend_common)
_score_layercluster_to_caloparticles = PlotGroup("score_layercluster_to_caloparticle", [
        Plot("Score_layercl2caloparticle_perlayer{:02d}".format(i+1), xtitle="Layer {:02d}".format(i+1), **_common_score) for i in range(0,maxlayer)
        ],
                                    ncols=8
                                    )

_common_shared= {"title": "Shared Energy CaloParticle To Layer Cluster",
                 "stat": False,
                 "legend": False,
                }
_common_shared.update(_legend_common)
_shared_plots = [Plot("SharedEnergy_caloparticle2layercl_perlayer{:02d}".format(i+1), xtitle="", **_common_shared) for i in range(0,maxlayer)]
_shared_plots.extend([Plot("SharedEnergy_caloparticle2layercl_vs_eta_perlayer{:02d}".format(i+1), xtitle="Layer {:02d}".format(i+1), **_common_shared) for i in range(0,maxlayer)])
_shared_plots.extend([Plot("SharedEnergy_caloparticle2layercl_vs_phi_perlayer{:02d}".format(i+1), xtitle="Layer {:02d}".format(i+1), **_common_shared) for i in range(0,maxlayer)])
_sharedEnergy_caloparticle_to_layercluster = PlotGroup("sharedEnergy_caloparticle_to_layercluster", _shared_plots, ncols=8)

_common_shared= {"title": "Shared Energy Layer Cluster To CaloParticle",
                 "stat": False,
                 "legend": False,
                }
_common_shared.update(_legend_common)
_shared_plots2 = [Plot("SharedEnergy_layercluster2caloparticle_perlayer{:02d}".format(i+1), xtitle="", **_common_shared) for i in range(0,maxlayer)]
_shared_plots2.extend([Plot("SharedEnergy_layercl2caloparticle_vs_eta_perlayer{:02d}".format(i+1), xtitle="Layer {:02d}".format(i+1), **_common_shared) for i in range(0,maxlayer)])
_shared_plots2.extend([Plot("SharedEnergy_layercl2caloparticle_vs_phi_perlayer{:02d}".format(i+1), xtitle="Layer {:02d}".format(i+1), **_common_shared) for i in range(0,maxlayer)])
_sharedEnergy_layercluster_to_caloparticle = PlotGroup("sharedEnergy_layercluster_to_caloparticle", _shared_plots2, ncols=8)


_common_assoc = {#"title": "Cell Association Table",
                 "stat": False,
                 "legend": False,
                 "xbinlabels": ["", "TN(pur)", "FN(ineff.)", "FP(fake)", "TP(eff)"],
                 "drawStyle": "hist",
                 "ymin": 0.1,
                 "ymax": 10000,
                 "ylog": True}
_common_assoc.update(_legend_common)
_cell_association_table = PlotGroup("cellAssociation_table", [
        Plot("cellAssociation_perlayer{:02d}".format(i+1), xtitle="Layer {:02d}".format(i+1), **_common_assoc) for i in range(0,maxlayer)
        ],
                                    ncols=8
                                    )

_bin_count = 0
_label_period = 4
_xbinlabels = ["%d"%(i+1) if (i+1)%_label_period==0 else "" for i in range(0,maxlayer)]
_xbinlabels.extend(["%d"%(i+1) if (i+1)%_label_period == 0 else "" for i in range(0,maxlayer)])
_common_eff = {"stat": False, "legend": False, "xbinlabels": _xbinlabels, "xbinlabelsize": 12, "xbinlabeloptions": "v"}
_effplots = [Plot("effic_eta_layer{:02d}".format(i+1), xtitle="", **_common_eff) for i in range(0,maxlayer)]
_effplots.extend([Plot("effic_phi_layer{:02d}".format(i+1), xtitle="", **_common_eff) for i in range(0,maxlayer)])
_common_eff["xmin"] = 0.
_bin_count += 52*2.
_common_eff["xmax"] =_bin_count
_effplots.extend([Plot("globalEfficiencies", xtitle="Global Efficiencies", **_common_eff)])
_efficiencies = PlotGroup("Efficiencies", _effplots, ncols=8)


_common_dup = {"stat": False, "legend": False, "title": "Global Duplicates", "xbinlabels": _xbinlabels, "xbinlabelsize": 12, "xbinlabeloptions": "v"}
_dupplots = [Plot("duplicate_eta_layer{:02d}".format(i+1), xtitle="", **_common_dup) for i in range(0,maxlayer)]
_dupplots.extend([Plot("duplicate_phi_layer{:02d}".format(i+1), xtitle="", **_common_dup) for i in range(0,maxlayer)])
_common_dup["xmin"] = _bin_count+1
_bin_count += 52*2.
_common_dup["xmax"] = _bin_count
_dupplots.extend([Plot("globalEfficiencies", xtitle="Global Duplicates", **_common_dup)])
_duplicates = PlotGroup("Duplicates", _dupplots, ncols=8)

_common_fake = {"stat": False, "legend": False, "title": "Global Fake Rates", "xbinlabels": _xbinlabels, "xbinlabelsize": 12, "xbinlabeloptions": "v"}
_fakeplots = [Plot("fake_eta_layer{:02d}".format(i+1), xtitle="", **_common_fake) for i in range(0,maxlayer)]
_fakeplots.extend([Plot("fake_phi_layer{:02d}".format(i+1), xtitle="", **_common_fake) for i in range(0,maxlayer)])
_common_fake["xmin"] = _bin_count+1
_bin_count += 52*2.
_common_fake["xmax"] = _bin_count
_fakeplots.extend([Plot("globalEfficiencies", xtitle="Global Fake Rate", **_common_fake)])
_fakes = PlotGroup("FakeRate", _fakeplots, ncols=8)

_common_merge = {"stat": False, "legend": False, "title": "Global Merge Rates", "xbinlabels": _xbinlabels, "xbinlabelsize": 12, "xbinlabeloptions": "v"}
_mergeplots = [Plot("merge_eta_layer{:02d}".format(i+1), xtitle="", **_common_merge) for i in range(0,maxlayer)]
_mergeplots.extend([Plot("merge_phi_layer{:02d}".format(i+1), xtitle="", **_common_merge) for i in range(0,maxlayer)])
_common_merge["xmin"] = _bin_count+1
_bin_count += 52*2.
_common_merge["xmax"] = _bin_count
_mergeplots.extend([Plot("globalEfficiencies", xtitle="Global merge Rate", **_common_merge)])
_merges = PlotGroup("MergeRate", _mergeplots, ncols=8)


_common_energy_score = dict(removeEmptyBins=False, xbinlabelsize=10, xbinlabeloption="d", ncols=4)
_energyscore_cp2lc = []
for i in range(0, maxlayer):
  _energyscore_cp2lc.append(PlotOnSideGroup("Energy_vs_Score_Layer{:02d}".format(i+1), Plot("Energy_vs_Score_caloparticle2layer_perlayer{:02d}".format(i+1), drawStyle="COLZ", adjustMarginLeft=0.1, adjustMarginRight=0.1, **_common_energy_score)))

_energyscore_lc2cp = []
for i in range(0, maxlayer):
  _energyscore_lc2cp.append(PlotOnSideGroup("Energy_vs_Score_Layer{:02d}".format(i+1), Plot("Energy_vs_Score_layer2caloparticle_perlayer{:02d}".format(i+1), drawStyle="COLZ", adjustMarginLeft=0.1, adjustMarginRight=0.1, **_common_energy_score)))
#_energyclustered =



hgcalLayerClustersPlotter = Plotter()
#We follow Chris categories in folders
# [A] calculated "energy density" for cells in a) 120um, b) 200um, c) 300um, d) scint
# (one entry per rechit, in the appropriate histo)
hgcalLayerClustersPlotter.append("CellsEnergyDensityPerThickness", [
        "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
        ], PlotFolder(
        _cellsenedens_thick,
        loopSubFolders=False,
        purpose=PlotPurpose.Timing, page="CellsEnergyDensityPerThickness"
        ))

# [B] number of layer clusters per event in a) 120um, b) 200um, c) 300um, d) scint
# (one entry per event in each of the four histos)
hgcalLayerClustersPlotter.append("TotalNumberofLayerClustersPerThickness", [
        "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
        ], PlotFolder(
        _totclusternum_thick,
        loopSubFolders=False,
        purpose=PlotPurpose.Timing, page="TotalNumberofLayerClustersPerThickness"
        ))

# [C] number of layer clusters per layer (one entry per event in each histo)
hgcalLayerClustersPlotter.append("NumberofLayerClustersPerLayer", [
        "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
        ], PlotFolder(
        _totclusternum_layer_EE,
        _totclusternum_layer_FH,
        _totclusternum_layer_BH,
        loopSubFolders=False,
        purpose=PlotPurpose.Timing, page="NumberofLayerClustersPerLayer"
        ))

# [D] For each layer cluster:
# number of cells in layer cluster, by layer - separate histos in each layer for 120um Si, 200/300um Si, Scint
# NB: not all combinations exist; e.g. no 120um Si in layers with scint.
# (One entry in the appropriate histo per layer cluster).
hgcalLayerClustersPlotter.append("CellsNumberPerLayerPerThickness", [
        "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
        ], PlotFolder(
        _cellsnum_perthick_perlayer_120_EE,
        _cellsnum_perthick_perlayer_120_FH,
        _cellsnum_perthick_perlayer_120_BH,
        _cellsnum_perthick_perlayer_200_EE,
        _cellsnum_perthick_perlayer_200_FH,
        _cellsnum_perthick_perlayer_200_BH,
        _cellsnum_perthick_perlayer_300_EE,
        _cellsnum_perthick_perlayer_300_FH,
        _cellsnum_perthick_perlayer_300_BH,
        _cellsnum_perthick_perlayer_scint_EE,
        _cellsnum_perthick_perlayer_scint_FH,
        _cellsnum_perthick_perlayer_scint_BH,
        loopSubFolders=False,
        purpose=PlotPurpose.Timing, page="CellsNumberPerLayerPerThickness"
        ))

# [E] For each layer cluster:
# distance of cells from a) seed cell, b) max cell; and c), d): same with entries weighted by cell energy
# separate histos in each layer for 120um Si, 200/300um Si, Scint
# NB: not all combinations exist; e.g. no 120um Si in layers with scint.
# (One entry in each of the four appropriate histos per cell in a layer cluster)
hgcalLayerClustersPlotter.append("CellsDistanceToSeedAndMaxCellPerLayerPerThickness", [
        "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
        ], PlotFolder(
        _distancetomaxcell_perthickperlayer_120_EE,
        _distancetomaxcell_perthickperlayer_120_FH,
        _distancetomaxcell_perthickperlayer_120_BH,
        _distancetomaxcell_perthickperlayer_200_EE,
        _distancetomaxcell_perthickperlayer_200_FH,
        _distancetomaxcell_perthickperlayer_200_BH,
        _distancetomaxcell_perthickperlayer_300_EE,
        _distancetomaxcell_perthickperlayer_300_FH,
        _distancetomaxcell_perthickperlayer_300_BH,
        _distancetomaxcell_perthickperlayer_scint_EE,
        _distancetomaxcell_perthickperlayer_scint_FH,
        _distancetomaxcell_perthickperlayer_scint_BH,
        _distancetoseedcell_perthickperlayer_120_EE,
        _distancetoseedcell_perthickperlayer_120_FH,
        _distancetoseedcell_perthickperlayer_120_BH,
        _distancetoseedcell_perthickperlayer_200_EE,
        _distancetoseedcell_perthickperlayer_200_FH,
        _distancetoseedcell_perthickperlayer_200_BH,
        _distancetoseedcell_perthickperlayer_300_EE,
        _distancetoseedcell_perthickperlayer_300_FH,
        _distancetoseedcell_perthickperlayer_300_BH,
        _distancetoseedcell_perthickperlayer_scint_EE,
        _distancetoseedcell_perthickperlayer_scint_FH,
        _distancetoseedcell_perthickperlayer_scint_BH,
        _distancetomaxcell_perthickperlayer_eneweighted_120_EE,
        _distancetomaxcell_perthickperlayer_eneweighted_120_FH,
        _distancetomaxcell_perthickperlayer_eneweighted_120_BH,
        _distancetomaxcell_perthickperlayer_eneweighted_200_EE,
        _distancetomaxcell_perthickperlayer_eneweighted_200_FH,
        _distancetomaxcell_perthickperlayer_eneweighted_200_BH,
        _distancetomaxcell_perthickperlayer_eneweighted_300_EE,
        _distancetomaxcell_perthickperlayer_eneweighted_300_FH,
        _distancetomaxcell_perthickperlayer_eneweighted_300_BH,
        _distancetomaxcell_perthickperlayer_eneweighted_scint_EE,
        _distancetomaxcell_perthickperlayer_eneweighted_scint_FH,
        _distancetomaxcell_perthickperlayer_eneweighted_scint_BH,
        _distancetoseedcell_perthickperlayer_eneweighted_120_EE,
        _distancetoseedcell_perthickperlayer_eneweighted_120_FH,
        _distancetoseedcell_perthickperlayer_eneweighted_120_BH,
        _distancetoseedcell_perthickperlayer_eneweighted_200_EE,
        _distancetoseedcell_perthickperlayer_eneweighted_200_FH,
        _distancetoseedcell_perthickperlayer_eneweighted_200_BH,
        _distancetoseedcell_perthickperlayer_eneweighted_300_EE,
        _distancetoseedcell_perthickperlayer_eneweighted_300_FH,
        _distancetoseedcell_perthickperlayer_eneweighted_300_BH,
        _distancetoseedcell_perthickperlayer_eneweighted_scint_EE,
        _distancetoseedcell_perthickperlayer_eneweighted_scint_FH,
        _distancetoseedcell_perthickperlayer_eneweighted_scint_BH,
        loopSubFolders=False,
        purpose=PlotPurpose.Timing, page="CellsDistanceToSeedAndMaxCellPerLayerPerThickness"
        ))

# [F] Looking at the fraction of true energy that has been clustered; by layer and overall
hgcalLayerClustersPlotter.append("EnergyClusteredByLayerAndOverall", [
        "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
        ], PlotFolder(
        _energyclustered,
        _energyclustered_perlayer_EE,
        _energyclustered_perlayer_FH,
        _energyclustered_perlayer_BH,
        loopSubFolders=False,
        purpose=PlotPurpose.Timing, page="EnergyClusteredByLayerAndOverall"
        ))

# [G] Miscellaneous plots:
# longdepthbarycentre: The longitudinal depth barycentre. One entry per event.
# mixedhitscluster: Number of clusters per event with hits in different thicknesses.
# num_reco_cluster_eta: Number of reco clusters vs eta

hgcalLayerClustersPlotter.append("Miscellaneous", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _num_reco_cluster_eta,
            _mixedhitscluster,
            _energyclustered,
            _longdepthbarycentre,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="Miscellaneous"
            ))

# [H] SelectedCaloParticles plots
hgcalLayerClustersPlotter.append("SelectedCaloParticles_Photons", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/SelectedCaloParticles/22",
            ], PlotFolder(
            _SelectedCaloParticles,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="SelectedCaloParticles_Photons"
            ))

# [I] Score of CaloParticles wrt Layer Clusters
hgcalLayerClustersPlotter.append("ScoreCaloParticlesToLayerClusters", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _score_caloparticle_to_layerclusters,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="ScoreCaloParticlesToLayerClusters"))

# [J] Score of LayerClusters wrt CaloParticles
hgcalLayerClustersPlotter.append("ScoreLayerClustersToCaloParticles", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _score_layercluster_to_caloparticles,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="ScoreLayerClustersToCaloParticles"))

# [K] Shared Energy between CaloParticle and LayerClusters
hgcalLayerClustersPlotter.append("SharedEnergy", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _sharedEnergy_caloparticle_to_layercluster,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="SharedEnergyCaloParticleToLayerCluster"))

# [K2] Shared Energy between LayerClusters and CaloParticle
hgcalLayerClustersPlotter.append("SharedEnergy", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _sharedEnergy_layercluster_to_caloparticle,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="SharedEnergyLayerClusterToCaloParticle"))

# [L] Cell Association per Layer
hgcalLayerClustersPlotter.append("CellAssociation", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _cell_association_table,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="CellAssociation"))

# [M] Efficiency Plots
hgcalLayerClustersPlotter.append("Efficiencies", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _efficiencies,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="Efficiencies"))
# [L] Duplicate Plots
hgcalLayerClustersPlotter.append("Duplicates", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _duplicates,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="Duplicates"))
# [M] Fake Rate Plots
hgcalLayerClustersPlotter.append("FakeRate", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _fakes,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="Fakes"))
# [N] Merge Rate Plots
hgcalLayerClustersPlotter.append("MergeRate", [
            "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
            ], PlotFolder(
            _merges,
            loopSubFolders=False,
            purpose=PlotPurpose.Timing, page="Merges"))
# [O] Energy vs Score 2D plots CP to LC
for i,item in enumerate(_energyscore_cp2lc, start=1):
  hgcalLayerClustersPlotter.append("Energy_vs_Score_CP2LC", [
              "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
              ], PlotFolder(
              item,
              loopSubFolders=False,
              purpose=PlotPurpose.Timing, page="Energy_vs_Score_CP2LC"))
# [P] Energy vs Score 2D plots LC to CP
for i,item in enumerate(_energyscore_lc2cp, start=1):
  hgcalLayerClustersPlotter.append("Energy_vs_Score_LC2CP", [
              "DQMData/Run 1/HGCAL/Run summary/HGCalValidator/hgcalLayerClusters",
              ], PlotFolder(
              item,
              loopSubFolders=False,
              purpose=PlotPurpose.Timing, page="Energy_vs_Score_LC2CP"))

