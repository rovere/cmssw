// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <alpaka/alpaka.hpp>

#include "DataFormats/HGCalReco/interface/HGCalSoACellsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoACellsDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAOutDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalTilesConstants.h"

#include "HGCalLayerClustersSoAAlgoWrapper.h"

#include "CLUEAlgoAlpaka.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  class HGCalLayerClustersSoAAlgoKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
        const unsigned int numer_of_clusters,
        const HGCalSoACellsDeviceCollection::ConstView input_rechits_soa,
        const HGCalSoAOutDeviceCollection::ConstView input_clusters_soa,
        HGCalSoAClustersDeviceCollection::View outputs) const {

      // make a strided loop over the kernel grid, covering up to "size" elements
      for (int32_t i : elements_with_stride(acc, input_rechits_soa.metadata().size())) {
        // Skip unassigned rechits
        if (input_clusters_soa[i].clusterIndex() == -1) {
          continue;
        }
        auto clIdx = input_clusters_soa[i].clusterIndex();
        alpaka::atomicAdd(acc, &outputs[clIdx].energy(), input_rechits_soa[i].weight());
      }
    }
  };


  void HGCalLayerClustersSoAAlgoWrapper::run(Queue& queue,
      const unsigned int size,
      const HGCalSoACellsDeviceCollection::ConstView input_rechits_soa,
      const HGCalSoAOutDeviceCollection::ConstView input_clusters_soa,
      HGCalSoAClustersDeviceCollection::View outputs) const {

    auto energy = cms::alpakatools::make_device_view<float>(alpaka::getDev(queue), outputs.energy(), size);
    alpaka::memset(queue, energy, 0x0);

    // use 64 items per group (this value is arbitrary, but it's a reasonable starting point)
    uint32_t items = 64;

    // use as many groups as needed to cover the whole problem
    uint32_t groups = divide_up_by(input_rechits_soa.metadata().size(), items);

    // map items to
    //   - threads with a single element per thread on a GPU backend
    //   - elements within a single thread on a CPU backend
    auto workDiv = make_workdiv<Acc1D>(groups, items);

    alpaka::exec<Acc1D>(queue, workDiv, HGCalLayerClustersSoAAlgoKernel{}, size, input_rechits_soa, input_clusters_soa, outputs);
  }
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
