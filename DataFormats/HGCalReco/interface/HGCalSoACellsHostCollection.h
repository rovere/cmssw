#ifndef DataFormats_HGCalReco_interface_HGCalSoACellsHostCollection_h
#define DataFormats_HGCalReco_interface_HGCalSoACellsHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoACells.h"


  // SoA with x, y, z, id fields in host memory
using HGCalSoACellsHostCollection = PortableHostCollection<HGCalCellsSoA>;


#endif  // DataFormats_HGCalReco_interface_HGCalSoACellsHostCollection_h