import FWCore.ParameterSet.Config as cms

# import CMS geometry
from Geometry.CMSCommonData.cmsIdealGeometry2015XML_cfi import XMLIdealGeometryESSource

# add our custom detector grouping to DDD
XMLIdealGeometryESSource.geomXMLFiles.extend(['SimTracker/TrackerMaterialAnalysis/data/trackingMaterialGroups.xml'])
