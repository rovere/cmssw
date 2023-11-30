#ifndef DataFormats_PortableTestObjects_interface_alpaka_HGCalSoAClustersServiceDeviceCollection_h
#define DataFormats_PortableTestObjects_interface_alpaka_HGCalSoAClustersServiceDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAClustersService.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

    using HGCalSoAClustersServiceDeviceCollection = PortableCollection<HGCalClustersServiceSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // DataFormats_PortableTestObjects_interface_alpaka_HGCalSoAClustersServiceDeviceCollection_h
