#ifndef DataFormats_PortableTestObjects_interface_alpaka_HGCalSoAOutDeviceCollection_h
#define DataFormats_PortableTestObjects_interface_alpaka_HGCalSoAOutDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAOut.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

    // SoA with delta, rho, weight, nearestHigher, clusterIndex, layer, isSeed, and cellsCount fields in device global memory
    using HGCalSoAOutDeviceCollection = PortableCollection<HGCalCellsOutSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // DataFormats_PortableTestObjects_interface_alpaka_HGCalSoAOutDeviceCollection_h
