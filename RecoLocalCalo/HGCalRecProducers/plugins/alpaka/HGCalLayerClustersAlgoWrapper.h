#ifndef RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersAlgoWrapper_h
#define RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersAlgoWrapper_h

#include "DataFormats/HGCalReco/interface/HGCalSoACellsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoACellsDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAOutDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalLayerClustersAlgoWrapper {
    public:
      void run(Queue& queue,
          const unsigned int size,
          const HGCalSoACellsDeviceCollection::ConstView inputs,
          HGCalSoAOutDeviceCollection::View outputs
          ) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersAlgoWrapper_h

