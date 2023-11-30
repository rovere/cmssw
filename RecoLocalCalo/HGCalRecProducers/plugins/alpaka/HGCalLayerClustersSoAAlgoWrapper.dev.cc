// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <alpaka/alpaka.hpp>

#include "DataFormats/HGCalReco/interface/HGCalSoACellsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoACellsDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAOutDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersServiceDeviceCollection.h"
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

  // Kernel to find the max for every cluster
  class HGCalLayerClustersSoAAlgoKernelPositionByHits {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
        ALPAKA_FN_ACC void operator()(TAcc const& acc,
            const unsigned int numer_of_clusters,
            float thresholdW0,
            float positionDeltaRho2,
            const HGCalSoACellsDeviceCollection::ConstView input_rechits_soa,
            const HGCalSoAOutDeviceCollection::ConstView input_clusters_soa,
            HGCalSoAClustersDeviceCollection::View outputs,
            HGCalSoAClustersServiceDeviceCollection::View outputs_service) const {

          // make a strided loop over the kernel grid, covering up to "size" elements
          for (int32_t hit_index : elements_with_stride(acc, input_rechits_soa.metadata().size())) {
            const int cluster_index = input_clusters_soa[hit_index].clusterIndex();

            // Bail out if you are not part of any cluster
            if (cluster_index == -1) {
              continue;
            }

            alpaka::atomicAdd(acc, &outputs_service[cluster_index].total_weight(), input_rechits_soa[hit_index].weight());
            // Read the current seed index, and the associated energy.
            int clusterSeed = outputs_service[cluster_index].maxEnergyIndex();
            float clusterEnergy = (clusterSeed == -1) ? 0. : input_rechits_soa[clusterSeed].weight();

            while (input_rechits_soa[hit_index].weight() > clusterEnergy) {
              // If output_service[cluster_index].maxEnergyIndex() did not change,
              // store the new value and exit the loop.  Otherwise return the value
              // that has been updated, and decide again if the maximum needs to be
              // updated.
              int seed = alpaka::atomicCas(acc, & outputs_service[cluster_index].maxEnergyIndex(), clusterSeed, hit_index);
              if (seed == hit_index) {
                // atomicCas has stored the new value.
                break;
              } else {
                // Update the seed index and re-read the associated energy.
                clusterSeed = seed;
                clusterEnergy = input_rechits_soa[clusterSeed].weight();
              }
            }  // CAS
          }    // elements_with_stride
        }      // operator()
  };

    // Real Kernel position
    class HGCalLayerClustersSoAAlgoKernelPositionByHits2 {
      public:
        template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
          ALPAKA_FN_ACC void operator()(TAcc const& acc,
              const unsigned int numer_of_clusters,
              float thresholdW0,
              float positionDeltaRho2,
              const HGCalSoACellsDeviceCollection::ConstView input_rechits_soa,
              const HGCalSoAOutDeviceCollection::ConstView input_clusters_soa,
              HGCalSoAClustersDeviceCollection::View outputs,
              HGCalSoAClustersServiceDeviceCollection::View outputs_service) const {

            // make a strided loop over the kernel grid, covering up to "size" elements
            for (int32_t hit_index : elements_with_stride(acc, input_rechits_soa.metadata().size())) {
              const int cluster_index = input_clusters_soa[hit_index].clusterIndex();
              const int max_energy_index = outputs_service[cluster_index].maxEnergyIndex();

              // Bail out if you are not part of any cluster
              if (cluster_index == -1) {
                continue;
              }

              /*
                 float total_weight_log = 0.f;
                 float reference_x = input_rechits_soa[maxEnergyIndex].dim1();
                 float reference_y = input_rechits_soa[maxEnergyIndex].dim2();
                 float reference_z = input_rechits_soa[maxEnergyIndex].dim3();
                 */
              //for silicon only just use 1+6 cells = 1.3cm for all thicknesses
              const float d1 = input_rechits_soa[hit_index].dim1() - input_rechits_soa[max_energy_index].dim1();
              const float d2 = input_rechits_soa[hit_index].dim2() - input_rechits_soa[max_energy_index].dim2();
              if ((d1 * d1 + d2 * d2) > positionDeltaRho2) {
                continue;
              }
              float Wi = std::max(thresholdW0 + std::log(input_rechits_soa[hit_index].weight() / outputs_service[cluster_index].total_weight()), 0.f);
              alpaka::atomicAdd(acc, &outputs[cluster_index].x(), input_rechits_soa[hit_index].dim1() * Wi);
              alpaka::atomicAdd(acc, &outputs[cluster_index].y(), input_rechits_soa[hit_index].dim2() * Wi);
              alpaka::atomicAdd(acc, &outputs_service[cluster_index].total_weight_log(), Wi);
            }  // elements_with_stride
          }    // operator()
    };

    class HGCalLayerClustersSoAAlgoKernelPositionByHits3 {
      public:
        template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
          ALPAKA_FN_ACC void operator()(TAcc const& acc,
              const unsigned int numer_of_clusters,
              float thresholdW0,
              float positionDeltaRho2,
              const HGCalSoACellsDeviceCollection::ConstView input_rechits_soa,
              const HGCalSoAOutDeviceCollection::ConstView input_clusters_soa,
              HGCalSoAClustersDeviceCollection::View outputs,
              HGCalSoAClustersServiceDeviceCollection::View outputs_service) const {

            // make a strided loop over the kernel grid, covering up to "size" elements
            for (int32_t cluster_index : elements_with_stride(acc, outputs.metadata().size())) {
              const int max_energy_index = outputs_service[cluster_index].maxEnergyIndex();

              float inv_tot_weight = 1.f / outputs_service[cluster_index].total_weight();
              outputs[cluster_index].x() *= inv_tot_weight;
              outputs[cluster_index].y() *= inv_tot_weight;
              outputs[cluster_index].z() = input_rechits_soa[max_energy_index].dim3();
            }  // elements_with_stride
          }    // operator()
    };

  void HGCalLayerClustersSoAAlgoWrapper::run(Queue& queue,
      const unsigned int size,
      float thresholdW0,
      float positionDeltaRho2,
      const HGCalSoACellsDeviceCollection::ConstView input_rechits_soa,
      const HGCalSoAOutDeviceCollection::ConstView input_clusters_soa,
      HGCalSoAClustersDeviceCollection::View outputs,
      HGCalSoAClustersServiceDeviceCollection::View outputs_service
      ) const {

    auto energy = cms::alpakatools::make_device_view<float>(alpaka::getDev(queue), outputs.energy(), size);
    alpaka::memset(queue, energy, 0x0);
    auto total_weight = cms::alpakatools::make_device_view<float>(alpaka::getDev(queue), outputs_service.total_weight(), size);
    alpaka::memset(queue, total_weight, 0x0);
    auto total_weight_log = cms::alpakatools::make_device_view<float>(alpaka::getDev(queue), outputs_service.total_weight_log(), size);
    alpaka::memset(queue, total_weight_log, 0x0);
    auto maxEnergyValue = cms::alpakatools::make_device_view<float>(alpaka::getDev(queue), outputs_service.maxEnergyValue(), size);
    alpaka::memset(queue, maxEnergyValue, 0x0);
    auto maxEnergyIndex = cms::alpakatools::make_device_view<int>(alpaka::getDev(queue), outputs_service.maxEnergyIndex(), size);
    alpaka::memset(queue, maxEnergyIndex, 0xff);

    // use 64 items per group (this value is arbitrary, but it's a reasonable starting point)
    uint32_t items = 64;

    // use as many groups as needed to cover the whole problem
    uint32_t groups = divide_up_by(input_rechits_soa.metadata().size(), items);

    // map items to
    //   - threads with a single element per thread on a GPU backend
    //   - elements within a single thread on a CPU backend
    auto workDiv = make_workdiv<Acc1D>(groups, items);

    alpaka::exec<Acc1D>(queue, workDiv, HGCalLayerClustersSoAAlgoKernelEnergy{}, size, input_rechits_soa, input_clusters_soa, outputs);

#if 0
    // From now on, we divide work by cluster...

    // use as many groups as needed to cover the whole problem
    uint32_t group_clusters = divide_up_by(size, items);
    auto workDivClusters = make_workdiv<Acc1D>(group_clusters, items);

    alpaka::exec<Acc1D>(queue, workDivClusters, HGCalLayerClustersSoAAlgoKernelPosition{},
        size, thresholdW0, positionDeltaRho2, input_rechits_soa, input_clusters_soa, outputs);
#endif
    alpaka::exec<Acc1D>(queue, workDiv, HGCalLayerClustersSoAAlgoKernelPositionByHits{},
        size, thresholdW0, positionDeltaRho2, input_rechits_soa, input_clusters_soa, outputs, outputs_service);
    alpaka::exec<Acc1D>(queue, workDiv, HGCalLayerClustersSoAAlgoKernelPositionByHits2{},
        size, thresholdW0, positionDeltaRho2, input_rechits_soa, input_clusters_soa, outputs, outputs_service);
    uint32_t group_clusters = divide_up_by(size, items);
    auto workDivClusters = make_workdiv<Acc1D>(group_clusters, items);
    alpaka::exec<Acc1D>(queue, workDivClusters, HGCalLayerClustersSoAAlgoKernelPositionByHits3{},
        size, thresholdW0, positionDeltaRho2, input_rechits_soa, input_clusters_soa, outputs, outputs_service);
  }
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
