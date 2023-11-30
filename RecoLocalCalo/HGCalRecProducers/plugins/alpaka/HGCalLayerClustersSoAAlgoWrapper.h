#ifndef RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersSoAAlgoWrapper_h
#define RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersSoAAlgoWrapper_h

#include "DataFormats/HGCalReco/interface/HGCalSoACellsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoACellsDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAOutDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersServiceDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalLayerClustersSoAAlgoWrapper {
    public:
      void run(Queue& queue,
          const unsigned int numer_of_clusters,
          float thresholdW0,
          float positionDeltaRho2,
          const HGCalSoACellsDeviceCollection::ConstView input_rechits_soa,
          const HGCalSoAOutDeviceCollection::ConstView input_clusters_soa,
          HGCalSoAClustersDeviceCollection::View outputs,
          HGCalSoAClustersServiceDeviceCollection::View outputs_service
          ) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersSoAAlgoWrapper_h

