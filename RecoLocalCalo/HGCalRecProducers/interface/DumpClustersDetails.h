#ifndef RecoLocalCalo_HGCalRecProducers_DumpClustersDetails_h
#define RecoLocalCalo_HGCalRecProducers_DumpClustersDetails_h

#include <vector>
#include <iostream>
#include <fmt/format.h>

namespace hgcalUtils {

  static bool sortByDetId(const std::pair<int, float>& pair1, const std::pair<int, float>& pair2) {
    return pair1.first < pair2.first;
  }


  class DumpClustersDetails {
    public:
      DumpClustersDetails() {};

      template<typename T>
        void dumpInfos(const T& clusters, bool dumpCellsDetId) const {
          int count=0;
          for (auto &i: clusters) {
            std::cout << fmt::format("Seed: {}, Idx: {}, x: {:a}, y: {:a}, z: {:a}, eta: {:a}, phi: {:a}",
                i.seed().rawId(),
                count++,
                i.x(),
                i.y(),
                i.z(),
                i.eta(),
                i.phi());
            if (dumpCellsDetId) {
              auto sorted = i.hitsAndFractions(); // copy...
              std::stable_sort(std::begin(sorted), std::end(sorted), sortByDetId);
              for (auto const & c : sorted) {
                std::cout << fmt::format(" ({}, {:a})", c.first, c.second) ;
              } // loop on hits and fractions
            }
            std::cout << std::endl;
          } // loop on clusters
        }
  };

    class DumpClustersSoADetails {
      public:
        DumpClustersSoADetails() {};

        template<typename T>
          void dumpInfos(const T& clustersSoA) const {
            for (int i = 0; i < clustersSoA->metadata().size(); ++i) {
              auto clusterSoAV = clustersSoA.view()[i];
              std::cout << fmt::format("Idx: {}, delta: {:.h}, rho: {:.h}, nearest: {}, clsIdx: {}, isSeed: {}",
                  i,
                  clusterSoAV.delta(),
                  clusterSoAV.rho(),
                  clusterSoAV.nearestHigher(),
                  clusterSoAV.clusterIndex(),
                  clusterSoAV.isSeed()
                  ) << std::endl;
            }
          }
    };

    class DumpCellsSoADetails {
      public:
        DumpCellsSoADetails() {};

        template<typename T>
          void dumpInfos(const T& cells) const {
            for (int i = 0; i < cells->metadata().size(); ++i) {
              auto cellSoAV = cells.view()[i];
              std::cout << fmt::format("Idx Cell: {}, x: {:a}, y: {:a}, layer: {}, weight: {:a}, sigmaNoise: {:a}, detid: {}",
                  i,
                  cellSoAV.dim1(),
                  cellSoAV.dim2(),
                  cellSoAV.layer(),
                  cellSoAV.weight(),
                  cellSoAV.sigmaNoise(),
                  cellSoAV.detid()
                  ) << std::endl;
            }
          }
    };

    class DumpLegacySoADetails {
      public:
        DumpLegacySoADetails() {};

        template<typename T>
          void dumpInfos(T & cells) const {
            for (unsigned int l = 0; l < cells.size(); l++) {
              for (unsigned int i = 0; i < cells.at(l).dim1.size(); ++i) {
                std::cout << fmt::format("Idx Cell: {}, x: {:a}, y: {:a}, layer: {}, weight: {:a}, sigmaNoise: {:a}, delta: {:a}, rho: {:a}, nearestHigher: {}, clsIdx: {}, isSeed: {}, detid: {}",
                    i,
                    cells.at(l).dim1.at(i),
                    cells.at(l).dim2.at(i),
                    l,
                    cells.at(l).weight.at(i),
                    cells.at(l).sigmaNoise.at(i),
                    cells.at(l).delta.at(i),
                    cells.at(l).rho.at(i),
                    cells.at(l).nearestHigher.at(i),
                    cells.at(l).clusterIndex.at(i),
                    cells.at(l).isSeed.at(i),
                    cells.at(l).detid.at(i)
                    ) << std::endl;
              }
            }
          }
    };

  using DumpClusters = DumpClustersDetails;
  using DumpClustersSoA = DumpClustersSoADetails;
  using DumpCellsSoA = DumpCellsSoADetails;
  using DumpLegacySoA = DumpLegacySoADetails;
}  // namespace hgcalUtils


#endif
