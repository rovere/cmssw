#ifndef RecoLocalCalo_HGCalRecProducers_DumpClustersDetails_h
#define RecoLocalCalo_HGCalRecProducers_DumpClustersDetails_h

#include <vector>
#include <iostream>
#include <fmt/format.h>

namespace hgcalUtils {

  template<int N>
    float roundToN(float value){
      static float scale = std::pow(10., N);
      return std::round(value * scale) / scale;
    };

  static bool sortByDetId(const std::pair<int, float>& pair1, const std::pair<int, float>& pair2) {
        return pair1.first < pair2.first;
  }


  template<int N>
  class DumpClustersDetails {
    public:
      DumpClustersDetails() {};

      template<typename T>
        void dumpInfos(const T& clusters, bool dumpCellsDetId) const {
          int count=0;
          for (auto &i: clusters) {
            std::cout << fmt::format("Seed: {}, Idx: {}, x: {}, y: {}, z: {}, eta: {}, phi: {}",
                i.seed().rawId(),
                count++,
                roundToN<N>(i.x()),
                roundToN<N>(i.y()),
                roundToN<N>(i.z()),
                roundToN<N>(i.eta()),
                roundToN<N>(i.phi()));
            if (dumpCellsDetId) {
              auto sorted = i.hitsAndFractions(); // copy...
              std::stable_sort(std::begin(sorted), std::end(sorted), sortByDetId);
              for (auto const & c : sorted) {
                std::cout << fmt::format(" ({}, {})", c.first, c.second) ;
              } // loop on hits and fractions
            }
            std::cout << std::endl;
          } // loop on clusters
        }
  };

  template<int N>
    class DumpClustersSoADetails {
      public:
        DumpClustersSoADetails() {};

        template<typename T>
          void dumpInfos(const T& clustersSoA) const {
            for (int i = 0; i < clustersSoA->metadata().size(); ++i) {
              auto clusterSoAV = clustersSoA.view()[i];
              std::cout << fmt::format("Idx: {}, delta: {:.{}f}, rho: {:.{}f}, nearest: {}, clsIdx: {}, isSeed: {}",
                  i,
                  clusterSoAV.delta(), std::numeric_limits<float>::max_digits10,
                  clusterSoAV.rho(), std::numeric_limits<float>::max_digits10,
                  clusterSoAV.nearestHigher(),
                  clusterSoAV.clusterIndex(),
                  clusterSoAV.isSeed()
                  ) << std::endl;
              /*
              std::cout << fmt::format("Idx: {}, delta: {}, rho: {}, nearest: {}, clsIdx: {}, isSeed: {}",
                  i,
                  roundToN<N>(clusterSoAV.delta()),
                  roundToN<N>(clusterSoAV.rho()),
                  roundToN<N>(clusterSoAV.nearestHigher()),
                  clusterSoAV.clusterIndex(),
                  clusterSoAV.isSeed()
                  ) << std::endl;
              */
            }
          }
    };

  template<int N>
    class DumpCellsSoADetails {
      public:
        DumpCellsSoADetails() {};

        template<typename T>
          void dumpInfos(const T& cells) const {
            for (int i = 0; i < cells->metadata().size(); ++i) {
              auto cellSoAV = cells.view()[i];
              std::cout << fmt::format("Idx Cell: {}, x: {}, y: {}, layer: {}, weight: {}, sigmaNoise: {}, detid: {}",
                  i,
                  roundToN<N>(cellSoAV.dim1()),
                  roundToN<N>(cellSoAV.dim2()),
                  cellSoAV.layer(),
                  roundToN<N>(cellSoAV.weight()),
                  roundToN<N>(cellSoAV.sigmaNoise()),
                  cellSoAV.detid()
                  ) << std::endl;
            }
          }
    };

  template<int N>
    class DumpLegacySoADetails {
      public:
        DumpLegacySoADetails() {};

        template<typename T>
          void dumpInfos(T & cells) const {
            for (unsigned int l = 0; l < cells.size(); l++) {
              for (unsigned int i = 0; i < cells.at(l).dim1.size(); ++i) {
                std::cout << fmt::format("Idx Cell: {}, x: {}, y: {}, layer: {}, weight: {}, sigmaNoise: {}, delta: {:.{}f}, rho: {:.{}f}, nearestHigher: {}, clsIdx: {}, isSeed: {}, detid: {}",
                    i,
                    roundToN<N>(cells.at(l).dim1.at(i)),
                    roundToN<N>(cells.at(l).dim2.at(i)),
                    l,
                    roundToN<N>(cells.at(l).weight.at(i)),
                    roundToN<N>(cells.at(l).sigmaNoise.at(i)),
                    cells.at(l).delta.at(i), std::numeric_limits<float>::max_digits10,
                    cells.at(l).rho.at(i), std::numeric_limits<float>::max_digits10,
                    cells.at(l).nearestHigher.at(i),
                    cells.at(l).clusterIndex.at(i),
                    cells.at(l).isSeed.at(i),
                    cells.at(l).detid.at(i)
                    ) << std::endl;
              }
            }
          }
    };

  using DumpClusters = DumpClustersDetails<6>;
  using DumpClustersSoA = DumpClustersSoADetails<6>;
  using DumpCellsSoA = DumpCellsSoADetails<6>;
  using DumpLegacySoA = DumpLegacySoADetails<6>;
}  // namespace hgcalUtils


#endif
