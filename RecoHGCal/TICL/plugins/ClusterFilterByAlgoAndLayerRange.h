// Authors: Marco Rovere - marco.rovere@cern.ch, Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 09/2020

#ifndef RecoHGCal_TICL_ClusterFilterByAlgoAndLayerRange_H__
#define RecoHGCal_TICL_ClusterFilterByAlgoAndLayerRange_H__

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "ClusterFilterBase.h"

#include <memory>
#include <utility>

// Filter clusters that belong to a specific algorithm
namespace ticl {
  class ClusterFilterByAlgoAndLayerRange final : public ClusterFilterBase {
  public:
    ClusterFilterByAlgoAndLayerRange(const edm::ParameterSet& ps)
        : ClusterFilterBase(ps),
          algo_number_(ps.getParameter<int>("algo_number")),
          min_layerId_(ps.getParameter<int>("min_layerId")),
          max_layerId_(ps.getParameter<int>("max_layerId")) {}
    ~ClusterFilterByAlgoAndLayerRange() override{};

    void filter(const std::vector<reco::CaloCluster>& layerClusters,
                const HgcalClusterFilterMask& availableLayerClusters,
                std::vector<float>& layerClustersMask,
                hgcal::RecHitTools& rhtools) const override {
      auto filteredLayerClusters = std::make_unique<HgcalClusterFilterMask>();
      for (auto const& cl : availableLayerClusters) {
        auto const& layerCluster = layerClusters[cl.first];
        auto layerId = rhtools.getLayerWithOffset(layerCluster.hitsAndFractions()[0].first);
        if (layerCluster.algo() == algo_number_ and layerId <= max_layerId_ and
            layerId >= min_layerId_ ) {
          filteredLayerClusters->emplace_back(cl);
        } else {
          layerClustersMask[cl.first] = 0.;
        }
      }
    }

  private:
    int algo_number_;
    unsigned int min_layerId_;
    unsigned int max_layerId_;
  };
}  // namespace ticl

#endif
