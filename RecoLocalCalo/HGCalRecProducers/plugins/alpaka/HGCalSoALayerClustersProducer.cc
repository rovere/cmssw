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
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalTilesConstants.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/HGCalReco/interface/HGCalSoACellsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAOutDeviceCollection.h"
#include "HGCalLayerClustersSoAAlgoWrapper.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalSoALayerClustersProducer : public stream::EDProducer<> {
    public:
      HGCalSoALayerClustersProducer(edm::ParameterSet const& config)
        : getTokenDeviceRecHits_{consumes(config.getParameter<edm::InputTag>("hgcalRecHitsSoA"))},
          getTokenDeviceClusters_{consumes(config.getParameter<edm::InputTag>("hgcalLayerClustersSoA"))},
          deviceTokenSoAClusters_{produces()}{}

      ~HGCalSoALayerClustersProducer() override = default;

      void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {

        // Get RecHitsSoA on the device
        auto const & deviceInputRecHits = iEvent.get(getTokenDeviceRecHits_);
        auto const inputRechits_v = deviceInputRecHits.view();

        // Get LayerClusters almost-SoA on device: this has still the same
        // cardinality as the RecHitsSoA, but has all the required information
        // to assemble the clusters, i.e., it has the cluster index assigned to
        // each rechit.
        auto const & deviceInputClusters = iEvent.get(getTokenDeviceClusters_);
        auto const inputClusters_v = deviceInputClusters.view();
        //
        // Allocate output SoA for the clusters, one entry for each cluster
        auto numclusters = cms::alpakatools::make_device_view<const unsigned int>(alpaka::getDev(iEvent.queue()), inputClusters_v.numberOfClustersScalar());
        unsigned int p;
        auto host_p = cms::alpakatools::make_host_view<unsigned int>(p);
        alpaka::memcpy(iEvent.queue(), host_p, numclusters);
        alpaka::wait(iEvent.queue());

        ALPAKA_ACCELERATOR_NAMESPACE::PortableCollection<HGCalClustersSoA> output(p, iEvent.queue());
        auto output_v = output.view();

        algo_.run(iEvent.queue(),
            p,
            inputRechits_v, inputClusters_v, output_v);
        iEvent.emplace(deviceTokenSoAClusters_, std::move(output));
      }


      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
        edm::ParameterSetDescription desc;
        desc.add<edm::InputTag>("hgcalLayerClustersSoA", edm::InputTag("TO BE DEFINED"));
        desc.add<edm::InputTag>("hgcalRecHitsSoA", edm::InputTag("TO BE DEFINED"));
        descriptions.addWithDefaultLabel(desc);
      }

    private:
      device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::PortableCollection<HGCalCellsSoA>> const getTokenDeviceRecHits_;
      device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::PortableCollection<HGCalCellsOutSoA>> const getTokenDeviceClusters_;
      device::EDPutToken<ALPAKA_ACCELERATOR_NAMESPACE::PortableCollection<HGCalClustersSoA>> const deviceTokenSoAClusters_;
      HGCalLayerClustersSoAAlgoWrapper algo_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalSoALayerClustersProducer);
