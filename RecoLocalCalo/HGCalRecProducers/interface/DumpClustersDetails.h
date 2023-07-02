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
          for (auto &i: clusters) {
            std::cout << fmt::format("Seed: {}, x: {}, y: {}, z: {}, eta: {}, phi: {}",
                i.seed().rawId(),
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

  using DumpClusters = DumpClustersDetails<6>;
}  // namespace hgcalUtils


#endif
