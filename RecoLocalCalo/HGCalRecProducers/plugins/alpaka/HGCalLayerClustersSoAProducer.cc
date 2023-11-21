#include "DataFormats/PortableTestObjects/interface/alpaka/TestDeviceCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaTest/interface/AlpakaESTestRecords.h"
#include "HeterogeneousCore/AlpakaTest/interface/alpaka/AlpakaESTestData.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/HGCalReco/interface/HGCalSoACellsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoACellsDeviceCollection.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalLayerClustersSoAProducer : public stream::EDProducer<> {
    public:
      HGCalLayerClustersSoAProducer(edm::ParameterSet const& config)
        : getTokenDevice_{consumes(config.getParameter<edm::InputTag>("hgcalRecHitsSoA"))},
      deviceToken_{produces()}{}

      ~HGCalLayerClustersSoAProducer() override = default;

      void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {

        auto const& deviceInput = iEvent.get(getTokenDevice_);
        std::cout << "Size of device collection: " << deviceInput->metadata().size() << std::endl;

        /*
        if constexpr (! std::is_same_v<ALPAKA_ACCELERATOR_NAMESPACE::Device, alpaka_common::DevHost>) {
          // Trigger copy async to GPU
          HGCalSoACellsDeviceCollection deviceProduct{cells->metadata().size(), iEvent.queue()}; // QUEUE TO BE VERIFIED
          alpaka::memcpy(iEvent.queue(), deviceProduct.buffer(), cells.const_buffer());
          iEvent.emplace(deviceToken_, std::move(deviceProduct));
        } else {
          iEvent.emplace(deviceToken_, std::move(cells));
        }
        */
      }

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
        edm::ParameterSetDescription desc;
        desc.add<edm::InputTag>("hgcalRecHitsSoA", edm::InputTag("TO BE DEFINED"));
        descriptions.addWithDefaultLabel(desc);
      }

    private:
      // use device::EDGetToken<T> to read from device memory space
      device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::PortableCollection<HGCalCellsSoA>> const getTokenDevice_;
      device::EDPutToken<ALPAKA_ACCELERATOR_NAMESPACE::PortableCollection<HGCalCellsSoA>> const deviceToken_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalLayerClustersSoAProducer);
