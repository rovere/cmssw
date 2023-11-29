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
      /*
      // global index of the thread within the grid
      const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];
      const portabletest::Matrix matrix{{1, 2, 3, 4, 5, 6}, {2, 4, 6, 8, 10, 12}, {3, 6, 9, 12, 15, 18}};

      // set this only once in the whole kernel grid
      if (thread == 0) {
        view.r() = 1.;
      }

      // make a strided loop over the kernel grid, covering up to "size" elements
      for (int32_t i : elements_with_stride(acc, size)) {
        view[i] = {xvalue, 0., 0., i, matrix * i};
      }
      */
    }
  };


  void HGCalLayerClustersSoAAlgoWrapper::run(Queue& queue,
      const unsigned int size,
      const HGCalSoACellsDeviceCollection::ConstView input_rechits_soa,
      const HGCalSoAOutDeviceCollection::ConstView input_clusters_soa,
      HGCalSoAClustersDeviceCollection::View outputs) const {

    /*
    // use 64 items per group (this value is arbitrary, but it's a reasonable starting point)
    uint32_t items = 64;

    // use as many groups as needed to cover the whole problem
    uint32_t groups = divide_up_by(collection->metadata().size(), items);

    // map items to
    //   - threads with a single element per thread on a GPU backend
    //   - elements within a single thread on a CPU backend
    auto workDiv = make_workdiv<Acc1D>(groups, items);

    alpaka::exec<Acc1D>(queue, workDiv, TestAlgoKernel{}, collection.view(), collection->metadata().size(), xvalue);
    */
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
