// Authors: Marco Rovere, Felice Pantaleo - marco.rovere@cern.ch, felice.pantaleo@cern.ch
// Date: 05/2019

#ifndef RecoHGCal_TICL_TICLLayerTile_h
#define RecoHGCal_TICL_TICLLayerTile_h

#include "DataFormats/TICL/interface/Common.h"
#include "DataFormats/Math/interface/normalizedPhi.h"

namespace ticl {
  class TICLLayerTile {
  public:
    void fill(double eta, double phi, unsigned int layerClusterId) {
      tile_[getGlobalBin(eta, phi)].push_back(layerClusterId);
    }

    int getEtaBin(float eta) const {
      constexpr float etaRange = ticl::constants::maxEta - ticl::constants::minEta;
      static_assert(etaRange >= 0.f);
      float r = ticl::constants::nEtaBins / etaRange;
      int etaBin = (std::abs(eta) - ticl::constants::minEta) * r;
      etaBin = std::clamp(etaBin, 0, ticl::constants::nEtaBins);
      return etaBin;
    }

    int getPhiBin(float phi) const {
      auto normPhi = normalizedPhi(phi);
      float r = ticl::constants::nPhiBins * M_1_PI * 0.5f;
      int phiBin = (normPhi + M_PI) * r;

      return phiBin;
    }

    int getGlobalBin(int etaBin, int phiBin) const { return phiBin + etaBin * ticl::constants::nPhiBins; }

    int getGlobalBin(double eta, double phi) const {
      return getPhiBin(phi) + getEtaBin(eta) * ticl::constants::nPhiBins;
    }

    void clear() {
      auto nBins = ticl::constants::nEtaBins * ticl::constants::nPhiBins;
      for (int j = 0; j < nBins; ++j)
        tile_[j].clear();
    }

    const std::vector<unsigned int>& operator[](int globalBinId) const { return tile_[globalBinId]; }

  private:
    Tile tile_;
  };
  typedef std::array<TICLLayerTile, ticl::constants::nLayers> Tiles;

  class TICLLayerTiles {
  public:
    const TICLLayerTile& operator[](int layer) const { return tiles_[layer]; }
    void fill(int layer, double eta, double phi, unsigned int layerClusterId) {
      tiles_[layer].fill(eta, phi, layerClusterId);
    }

  private:
    Tiles tiles_;
  };
}  // namespace ticl

#endif
