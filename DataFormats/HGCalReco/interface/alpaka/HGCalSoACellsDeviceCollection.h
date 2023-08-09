#ifndef DataFormats_PortableTestObjects_interface_alpaka_HGCalSoACellsDeviceCollection_h
#define DataFormats_PortableTestObjects_interface_alpaka_HGCalSoACellsDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoACells.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

    // SoA with x, y, weight, sigmaNoise, recHitsIndex layer fields in device global memory
    using HGCalSoACellsDeviceCollection = PortableCollection<HGCalCellsSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // DataFormats_PortableTestObjects_interface_alpaka_HGCalSoACellsDeviceCollection_h
