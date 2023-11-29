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

  class HGCalLayerClustersSoAAlgoKernelEnergy {
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

  class HGCalLayerClustersSoAAlgoKernelPosition {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
        const unsigned int numer_of_clusters,
        float thresholdW0,
        float positionDeltaRho2,
        const HGCalSoACellsDeviceCollection::ConstView input_rechits_soa,
        const HGCalSoAOutDeviceCollection::ConstView input_clusters_soa,
        HGCalSoAClustersDeviceCollection::View outputs) const {

      // global index of the thread within the grid
      const int32_t thread_idx = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];

      outputs[thread_idx].x() = 0.f;
      outputs[thread_idx].y() = 0.f;
      outputs[thread_idx].z() = 0.f;
      float maxEnergyValue = 0.f;
      float total_weight = 0.f;
      unsigned int maxEnergyIndex = 0;

      for (auto hit = 0; hit < input_rechits_soa.metadata().size(); ++hit) {
        if (input_clusters_soa[hit].clusterIndex() != thread_idx) {
          continue;
        }
        total_weight += input_rechits_soa[hit].weight();
        if (input_rechits_soa[hit].weight() > maxEnergyValue) {
          maxEnergyValue = input_rechits_soa[hit].weight();
          maxEnergyIndex = hit;
        }
      }

      float total_weight_log = 0.f;
      float reference_x = input_rechits_soa[maxEnergyIndex].dim1();
      float reference_y = input_rechits_soa[maxEnergyIndex].dim2();
      float reference_z = input_rechits_soa[maxEnergyIndex].dim3();
      for (auto hit = 0; hit < input_rechits_soa.metadata().size(); ++hit) {
        if (input_clusters_soa[hit].clusterIndex() != thread_idx) {
          continue;
        }
        //for silicon only just use 1+6 cells = 1.3cm for all thicknesses
        const float d1 = input_rechits_soa[hit].dim1() - reference_x;
        const float d2 = input_rechits_soa[hit].dim2() - reference_y;
        if ((d1 * d1 + d2 * d2) > positionDeltaRho2) {
          continue;
        }
        float Wi = std::max(thresholdW0 + std::log(input_rechits_soa[hit].weight() / total_weight), 0.f);
        outputs[thread_idx].x() += input_rechits_soa[hit].dim1() * Wi;
        outputs[thread_idx].y() += input_rechits_soa[hit].dim2() * Wi;
        total_weight_log += Wi;
      }
      total_weight = total_weight_log;
      float inv_tot_weight = 1.f / total_weight;
      outputs[thread_idx].x() *= inv_tot_weight;
      outputs[thread_idx].y() *= inv_tot_weight;
      outputs[thread_idx].z() = reference_z;
    }
  };


  void HGCalLayerClustersSoAAlgoWrapper::run(Queue& queue,
      const unsigned int size,
      float thresholdW0,
      float positionDeltaRho2,
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

    alpaka::exec<Acc1D>(queue, workDiv, HGCalLayerClustersSoAAlgoKernelEnergy{}, size, input_rechits_soa, input_clusters_soa, outputs);

    // From now on, we divide work by cluster...

    // use as many groups as needed to cover the whole problem
    uint32_t group_clusters = divide_up_by(size, items);
    auto workDivClusters = make_workdiv<Acc1D>(group_clusters, items);

    alpaka::exec<Acc1D>(queue, workDivClusters, HGCalLayerClustersSoAAlgoKernelPosition{},
        size, thresholdW0, positionDeltaRho2, input_rechits_soa, input_clusters_soa, outputs);
  }
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
