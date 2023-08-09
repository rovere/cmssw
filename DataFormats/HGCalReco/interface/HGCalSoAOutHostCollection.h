#ifndef DataFormats_HGCalReco_interface_HGCalSoAOutHostCollection_h
#define DataFormats_HGCalReco_interface_HGCalSoAOutHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAOut.h"


// SoA with delta, rho, nearestHigher, clusterIndex, isSeed, and numberOfClusters fields in host memory
using HGCalSoAOutHostCollection = PortableHostCollection<HGCalCellsOutSoA>;

#endif  // DataFormats_HGCalReco_interface_HGCalSoAOutHostCollection_h
