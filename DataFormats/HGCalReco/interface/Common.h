#ifndef DataFormats_HGCalReco_Common_h
#define DataFormats_HGCalReco_Common_h

#include <vector>
#include <array>
#include <cstdint>

namespace ticl::constants {
  constexpr float minEta = 1.5f;
  constexpr float maxEta = 3.2f;
  constexpr int nEtaBins = 34;
  constexpr int nPhiBins = 126;
  constexpr int nLayers = 104;
  constexpr int iterations = 4;
}  // namespace ticl::constants

namespace ticl {
  typedef std::vector<std::pair<unsigned int, float> > HgcalClusterFilterMask;
  typedef std::array<std::vector<unsigned int>, constants::nEtaBins * constants::nPhiBins> Tile;
}  // namespace ticl

#endif  // DataFormats_HGCalReco_Common_h
