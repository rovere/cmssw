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

#include "HGCalLayerClustersAlgoWrapper.h"

#include "CLUEAlgoAlpaka.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  /*
  class TestAlgoKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  portabletest::TestDeviceCollection::View view,
                                  int32_t size,
                                  double xvalue) const {
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
    }
  };
  */

  void HGCalLayerClustersAlgoWrapper::run(Queue& queue,
      const unsigned int size,
      const HGCalSoACellsDeviceCollection::ConstView inputs,
      HGCalSoAOutDeviceCollection::View outputs) const {

    CLUEAlgoAlpaka<ALPAKA_ACCELERATOR_NAMESPACE::Acc1D, Queue,
    HGCalSiliconTilesConstants, 96> algoStandalone(queue, 1.3f,9.f,2.f,false);

    /*
    if constexpr (std::is_same_v<ALPAKA_ACCELERATOR_NAMESPACE::Device, alpaka_common::DevHost>) {
      std::cout << "Collection from HGCalLayerClustersAlgoWrapper@CPU " << inputs.metadata().size() << std::endl;
    }
    */
    algoStandalone.makeClustersCMSSW(size,
        inputs.dim1(),
        inputs.dim2(),
        inputs.layer(),
        inputs.weight(),
        inputs.sigmaNoise(),
        inputs.detid(),
        outputs.rho(),
        outputs.delta(),
        outputs.nearestHigher(),
        outputs.clusterIndex(),
        outputs.isSeed(),
        &outputs.numberOfClustersScalar());
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
