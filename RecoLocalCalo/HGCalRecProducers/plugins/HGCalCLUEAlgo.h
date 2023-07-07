#ifndef RecoLocalCalo_HGCalRecProducers_HGCalCLUEAlgo_h
#define RecoLocalCalo_HGCalRecProducers_HGCalCLUEAlgo_h

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalLayerTiles.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalCLUEStrategy.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

// C/C++ headers
#include <set>
#include <string>
#include <vector>
#include <cmath>

template <typename TILE, typename STRATEGY>
class HGCalCLUEAlgoT : public HGCalClusteringAlgoBase {
public:
  HGCalCLUEAlgoT(const edm::ParameterSet& ps)
      : HGCalClusteringAlgoBase(
            (HGCalClusteringAlgoBase::VerbosityLevel)ps.getUntrackedParameter<unsigned int>("verbosity", 3),
            reco::CaloCluster::undefined),
        vecDeltas_(ps.getParameter<std::vector<double>>("deltac")),
        kappa_(ps.getParameter<double>("kappa")),
        use2x2_(ps.getParameter<bool>("use2x2")),
        initialized_(false) {}

  ~HGCalCLUEAlgoT() override {}

  void getEventSetupPerAlgorithm(const edm::EventSetup& es) override;

  void populate(const HGCRecHitCollection& hits)override{};
  void setCellsOnLayer(std::shared_ptr<std::vector<CellsOnLayer>>  cells)override{
    cells_ = cells;
  }

  // this is the method that will start the clusterisation (it is possible to invoke this method
  // more than once - but make sure it is with different hit collections (or else use reset)

  void makeClusters() override;

  // this is the method to get the cluster collection out
  std::vector<reco::BasicCluster> getClusters(bool) override;

  void reset() override {
    clusters_v_.clear();
    clusters_v_.shrink_to_fit();
    for (auto& cl : numberOfClustersPerLayer_) {
      cl = 0;
    }

    for (auto& cells : *cells_) {
      cells.clear();
      cells.shrink_to_fit();
    }
  }

  static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    iDesc.add<std::vector<double>>("thresholdW0", {2.9, 2.9, 2.9});
    iDesc.add<double>("positionDeltaRho2", 1.69);
    iDesc.add<std::vector<double>>("deltac",
                                   {
                                       1.3,
                                       1.3,
                                       1.3,
                                       0.0315,  // for scintillator
                                   });
    iDesc.add<bool>("dependSensor", true);
    iDesc.add<double>("ecut", 3.0);
    iDesc.add<double>("kappa", 9.0);
    iDesc.addUntracked<unsigned int>("verbosity", 3);
    iDesc.add<std::vector<double>>("dEdXweights", {});
    iDesc.add<std::vector<double>>("thicknessCorrection", {});
    iDesc.add<double>("sciThicknessCorrection", 0.9);
    iDesc.add<int>("deltasi_index_regemfac", 3);
    iDesc.add<unsigned>("maxNumberOfThickIndices", 6);
    iDesc.add<std::vector<double>>("fcPerMip", {});
    iDesc.add<double>("fcPerEle", 0.0);
    iDesc.add<std::vector<double>>("noises", {});
    edm::ParameterSetDescription descNestedNoiseMIP;
    descNestedNoiseMIP.add<bool>("scaleByDose", false);
    descNestedNoiseMIP.add<unsigned int>("scaleByDoseAlgo", 0);
    descNestedNoiseMIP.add<double>("scaleByDoseFactor", 1.);
    descNestedNoiseMIP.add<std::string>("doseMap", "");
    descNestedNoiseMIP.add<std::string>("sipmMap", "");
    descNestedNoiseMIP.add<double>("referenceIdark", -1);
    descNestedNoiseMIP.add<double>("referenceXtalk", -1);
    descNestedNoiseMIP.add<double>("noise_MIP", 1. / 100.);
    iDesc.add<edm::ParameterSetDescription>("noiseMip", descNestedNoiseMIP);
    iDesc.add<bool>("use2x2", true);  // use 2x2 or 3x3 scenario for scint density calculation
  }

  /// point in the space
  typedef math::XYZPoint Point;

private:
  // The two parameters used to identify clusters
  std::vector<double> vecDeltas_;
  double kappa_;

  // various parameters used for calculating the noise levels for a given sensor (and whether to use
  // them)
  std::vector<double> dEdXweights_;
  double sciThicknessCorrection_;

  bool use2x2_;

  // initialization bool
  bool initialized_;

  float outlierDeltaFactor_ = 2.f;


  std::shared_ptr<std::vector<CellsOnLayer>> cells_;

  std::vector<int> numberOfClustersPerLayer_;

  inline float distance2(const TILE& lt, int cell1, int cell2, int layerId) const {  // 2-d distance on the layer (x-y)
    return (lt.distance2(cells_->at(layerId).dim1[cell1],
                                  cells_->at(layerId).dim2[cell1],
                                  cells_->at(layerId).dim1[cell2],
                                  cells_->at(layerId).dim2[cell2]));
  }

  inline float distance(const TILE& lt, int cell1, int cell2, int layerId) const {  // 2-d distance on the layer (x-y)
    return std::sqrt(lt.distance2(cells_->at(layerId).dim1[cell1],
                                  cells_->at(layerId).dim2[cell1],
                                  cells_->at(layerId).dim1[cell2],
                                  cells_->at(layerId).dim2[cell2]));
  }

  void prepareDataStructures(const unsigned int layerId);
  void calculateLocalDensity(const TILE& lt, const unsigned int layerId,
                             float delta);  // return max density
  void calculateLocalDensity(const TILE& lt, const unsigned int layerId, float delta, HGCalSiliconStrategy strategy);
  void calculateLocalDensity(const TILE& lt,
                             const unsigned int layerId,
                             float delta,
                             HGCalScintillatorStrategy strategy);
  void calculateDistanceToHigher(const TILE& lt, const unsigned int layerId, float delta);
  int findAndAssignClusters(const unsigned int layerId, float delta);
};

// explicit template instantiation
extern template class HGCalCLUEAlgoT<HGCalSiliconLayerTiles, HGCalSiliconStrategy>;
extern template class HGCalCLUEAlgoT<HGCalScintillatorLayerTiles, HGCalScintillatorStrategy>;
extern template class HGCalCLUEAlgoT<HFNoseLayerTiles, HGCalSiliconStrategy>;

using HGCalSiCLUEAlgo = HGCalCLUEAlgoT<HGCalSiliconLayerTiles, HGCalSiliconStrategy>;
using HGCalSciCLUEAlgo = HGCalCLUEAlgoT<HGCalScintillatorLayerTiles, HGCalScintillatorStrategy>;
using HFNoseCLUEAlgo = HGCalCLUEAlgoT<HFNoseLayerTiles, HGCalSiliconStrategy>;

#endif
