#include "RecoLocalCalo/HGCalRecProducers/plugins/HGCalCLUEAlgo.h"

// Geometry
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/DumpClustersDetails.h"
//
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "oneapi/tbb/task_arena.h"
#include "oneapi/tbb.h"
#include <limits>
#include "DataFormats/DetId/interface/DetId.h"

using namespace hgcal_clustering;

template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::getEventSetupPerAlgorithm(const edm::EventSetup& es) {
  numberOfClustersPerLayer_.clear();
  numberOfClustersPerLayer_.resize(2 * (maxlayer_ + 1), 0);
}

template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::prepareDataStructures(unsigned int l) {
  auto cellsSize = cells_->at(l).detid.size();
  cells_->at(l).rho.resize(cellsSize, 0.f);
  cells_->at(l).delta.resize(cellsSize, 9999999);
  cells_->at(l).nearestHigher.resize(cellsSize, -1);
  cells_->at(l).clusterIndex.resize(cellsSize, -1);
  cells_->at(l).followers.resize(cellsSize);
  cells_->at(l).isSeed.resize(cellsSize, false);
}

// Create a vector of Hexels associated to one cluster from a collection of
// HGCalRecHits - this can be used directly to make the final cluster list -
// this method can be invoked multiple times for the same event with different
// input (reset should be called between events)
template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::makeClusters() {
  // assign all hits in each layer to a cluster core
  tbb::this_task_arena::isolate([&] {
    tbb::parallel_for(size_t(0), size_t(2 * maxlayer_ + 2), [&](size_t i) {
      prepareDataStructures(i);
      T lt;
      lt.clear();
      lt.fill(cells_->at(i).dim1, cells_->at(i).dim2);

      float delta;
      if constexpr (std::is_same_v<STRATEGY, HGCalSiliconStrategy>) {
        // maximum search distance (critical distance) for local density calculation
        float delta_c;
        if (i % maxlayer_ < lastLayerEE_)
          delta_c = vecDeltas_[0];
        else if (i % maxlayer_ < (firstLayerBH_ - 1))
          delta_c = vecDeltas_[1];
        else
          delta_c = vecDeltas_[2];
        delta = delta_c;
      } else {
        float delta_r = vecDeltas_[3];
        delta = delta_r;
      }
      LogDebug("HGCalCLUEAlgo") << "maxlayer: " << maxlayer_ << " lastLayerEE: " << lastLayerEE_
                                << " firstLayerBH: " << firstLayerBH_ << "\n";

      calculateLocalDensity(lt, i, delta);
      calculateDistanceToHigher(lt, i, delta);
      numberOfClustersPerLayer_[i] = findAndAssignClusters(i, delta);
    });
  });
}

template <typename T, typename STRATEGY>
std::vector<reco::BasicCluster> HGCalCLUEAlgoT<T, STRATEGY>::getClusters(bool) {
  std::vector<int> offsets(numberOfClustersPerLayer_.size(), 0);

  int maxClustersOnLayer = numberOfClustersPerLayer_[0];

  for (unsigned layerId = 1; layerId < offsets.size(); ++layerId) {
    offsets[layerId] = offsets[layerId - 1] + numberOfClustersPerLayer_[layerId - 1];
    maxClustersOnLayer = std::max(maxClustersOnLayer, numberOfClustersPerLayer_[layerId]);
  }

  auto totalNumberOfClusters = offsets.back() + numberOfClustersPerLayer_.back();
  clusters_v_.resize(totalNumberOfClusters);
  std::vector<std::vector<int>> cellsIdInCluster;
  cellsIdInCluster.reserve(maxClustersOnLayer);

  for (unsigned int layerId = 0; layerId < 2 * maxlayer_ + 2; ++layerId) {
    cellsIdInCluster.resize(numberOfClustersPerLayer_[layerId]);
    auto& cellsOnLayer = cells_->at(layerId);
    unsigned int numberOfCells = cellsOnLayer.detid.size();
    auto firstClusterIdx = offsets[layerId];

    for (unsigned int i = 0; i < numberOfCells; ++i) {
      auto clusterIndex = cellsOnLayer.clusterIndex[i];
      if (clusterIndex != -1)
        cellsIdInCluster[clusterIndex].push_back(i);
    }

    std::vector<std::pair<DetId, float>> thisCluster;

    for (auto& cl : cellsIdInCluster) {
      math::XYZPoint position = math::XYZPoint(0.f, 0.f, 0.f);
      float energy = 0.f;
      int seedDetId = -1;

      for (auto cellIdx : cl) {
        energy += cellsOnLayer.weight[cellIdx];
        thisCluster.emplace_back(cellsOnLayer.detid[cellIdx], 1.f);
        if (cellsOnLayer.isSeed[cellIdx]) {
          seedDetId = cellsOnLayer.detid[cellIdx];
        }
      }
      auto globalClusterIndex = cellsOnLayer.clusterIndex[cl[0]] + firstClusterIdx;

      clusters_v_[globalClusterIndex] =
          reco::BasicCluster(energy, position, reco::CaloID::DET_HGCAL_ENDCAP, thisCluster, algoId_);
      clusters_v_[globalClusterIndex].setSeed(seedDetId);
      thisCluster.clear();
    }

    cellsIdInCluster.clear();
  }
  return clusters_v_;
}
template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::calculateLocalDensity(const T& lt,
                                                        const unsigned int layerId,
                                                        float delta,
                                                        HGCalSiliconStrategy strategy) {
  auto& cellsOnLayer = cells_->at(layerId);
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  for (unsigned int i = 0; i < numberOfCells; i++) {
    std::array<int, 4> search_box = lt.searchBox(cellsOnLayer.dim1[i] - delta,
                                                 cellsOnLayer.dim1[i] + delta,
                                                 cellsOnLayer.dim2[i] - delta,
                                                 cellsOnLayer.dim2[i] + delta);

    for (int xBin = search_box[0]; xBin < search_box[1] + 1; ++xBin) {
      for (int yBin = search_box[2]; yBin < search_box[3] + 1; ++yBin) {
        int binId = lt.getGlobalBinByBin(xBin, yBin);
        size_t binSize = lt[binId].size();

        for (unsigned int j = 0; j < binSize; j++) {
          unsigned int otherId = lt[binId][j];
          if (distance2(lt, i, otherId, layerId) < delta*delta) {
            cellsOnLayer.rho[i] += (i == otherId ? 1.f : 0.5f) * cellsOnLayer.weight[otherId];
          }
        }
      }
    }
    LogDebug("HGCalCLUEAlgo") << "Debugging calculateLocalDensity: \n"
                              << "  cell: " << i << " eta: " << cellsOnLayer.dim1[i] << " phi: " << cellsOnLayer.dim2[i]
                              << " energy: " << cellsOnLayer.weight[i] << " density: " << cellsOnLayer.rho[i] << "\n";
  }
}
template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::calculateLocalDensity(const T& lt,
                                                        const unsigned int layerId,
                                                        float delta,
                                                        HGCalScintillatorStrategy strategy) {
  auto& cellsOnLayer = cells_->at(layerId);
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  for (unsigned int i = 0; i < numberOfCells; i++) {
    std::array<int, 4> search_box = lt.searchBox(cellsOnLayer.dim1[i] - delta,
                                                 cellsOnLayer.dim1[i] + delta,
                                                 cellsOnLayer.dim2[i] - delta,
                                                 cellsOnLayer.dim2[i] + delta);
    cellsOnLayer.rho[i] += cellsOnLayer.weight[i];
    float northeast(0), northwest(0), southeast(0), southwest(0), all(0);
    for (int etaBin = search_box[0]; etaBin < search_box[1] + 1; ++etaBin) {
      for (int phiBin = search_box[2]; phiBin < search_box[3] + 1; ++phiBin) {
        int phi = (phiBin % T::type::nRows);
        int binId = lt.getGlobalBinByBin(etaBin, phi);
        size_t binSize = lt[binId].size();
        for (unsigned int j = 0; j < binSize; j++) {
          unsigned int otherId = lt[binId][j];
          if (distance(lt, i, otherId, layerId) < delta) {
            int iPhi = HGCScintillatorDetId(cellsOnLayer.detid[i]).iphi();
            int otherIPhi = HGCScintillatorDetId(cellsOnLayer.detid[otherId]).iphi();
            int iEta = HGCScintillatorDetId(cellsOnLayer.detid[i]).ieta();
            int otherIEta = HGCScintillatorDetId(cellsOnLayer.detid[otherId]).ieta();
            int dIPhi = otherIPhi - iPhi;
            dIPhi += abs(dIPhi) < 2 ? 0
                     : dIPhi < 0    ? scintMaxIphi_
                                    : -scintMaxIphi_;  // cells with iPhi=288 and iPhi=1 should be neiboring cells
            int dIEta = otherIEta - iEta;
            LogDebug("HGCalCLUEAlgo") << "  Debugging calculateLocalDensity for Scintillator: \n"
                                      << "    cell: " << otherId << " energy: " << cellsOnLayer.weight[otherId]
                                      << " otherIPhi: " << otherIPhi << " iPhi: " << iPhi << " otherIEta: " << otherIEta
                                      << " iEta: " << iEta << "\n";

            if (otherId != i) {
              auto neighborCellContribution = 0.5f * cellsOnLayer.weight[otherId];
              all += neighborCellContribution;
              if (dIPhi >= 0 && dIEta >= 0)
                northeast += neighborCellContribution;
              if (dIPhi <= 0 && dIEta >= 0)
                southeast += neighborCellContribution;
              if (dIPhi >= 0 && dIEta <= 0)
                northwest += neighborCellContribution;
              if (dIPhi <= 0 && dIEta <= 0)
                southwest += neighborCellContribution;
            }
            LogDebug("HGCalCLUEAlgo") << "  Debugging calculateLocalDensity for Scintillator: \n"
                                      << "    northeast: " << northeast << " southeast: " << southeast
                                      << " northwest: " << northwest << " southwest: " << southwest << "\n";
          }
        }
      }
    }
    float neighborsval = (std::max(northeast, northwest) > std::max(southeast, southwest))
                             ? std::max(northeast, northwest)
                             : std::max(southeast, southwest);
    if (use2x2_)
      cellsOnLayer.rho[i] += neighborsval;
    else
      cellsOnLayer.rho[i] += all;
    LogDebug("HGCalCLUEAlgo") << "Debugging calculateLocalDensity: \n"
                              << "  cell: " << i << " eta: " << cellsOnLayer.dim1[i] << " phi: " << cellsOnLayer.dim2[i]
                              << " energy: " << cellsOnLayer.weight[i] << " density: " << cellsOnLayer.rho[i] << "\n";
  }
}
template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::calculateLocalDensity(const T& lt, const unsigned int layerId, float delta) {
  if constexpr (std::is_same_v<STRATEGY, HGCalSiliconStrategy>) {
    calculateLocalDensity(lt, layerId, delta, HGCalSiliconStrategy());
  } else {
    calculateLocalDensity(lt, layerId, delta, HGCalScintillatorStrategy());
  }
}

template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::calculateDistanceToHigher(const T& lt, const unsigned int layerId, float delta) {
  auto& cellsOnLayer = cells_->at(layerId);
  unsigned int numberOfCells = cellsOnLayer.detid.size();

  for (unsigned int i = 0; i < numberOfCells; i++) {
    // initialize delta and nearest higher for i
    float maxDelta = std::numeric_limits<float>::max();
    float i_delta = maxDelta;
    float rho_max = 0.f;
    int i_nearestHigher = -1;
    auto range = outlierDeltaFactor_ * delta;
    std::array<int, 4> search_box = lt.searchBox(cellsOnLayer.dim1[i] - range,
                                                 cellsOnLayer.dim1[i] + range,
                                                 cellsOnLayer.dim2[i] - range,
                                                 cellsOnLayer.dim2[i] + range);
    // loop over all bins in the search box
    for (int dim1Bin = search_box[0]; dim1Bin < search_box[1] + 1; ++dim1Bin) {
      for (int dim2Bin = search_box[2]; dim2Bin < search_box[3] + 1; ++dim2Bin) {
        // get the id of this bin
        size_t binId = lt.getGlobalBinByBin(dim1Bin, dim2Bin);
        if constexpr (std::is_same_v<STRATEGY, HGCalScintillatorStrategy>)
          binId = lt.getGlobalBinByBin(dim1Bin, (dim2Bin % T::type::nRows));
        // get the size of this bin
        size_t binSize = lt[binId].size();

        // loop over all hits in this bin
        for (unsigned int j = 0; j < binSize; j++) {
          unsigned int otherId = lt[binId][j];
          float dist = distance2(lt, i, otherId, layerId);
          bool foundHigher =
              (cellsOnLayer.rho[otherId] > cellsOnLayer.rho[i]) ||
              (cellsOnLayer.rho[otherId] == cellsOnLayer.rho[i] && cellsOnLayer.detid[otherId] > cellsOnLayer.detid[i]);
          
          if (foundHigher && dist < i_delta) {
            rho_max = cellsOnLayer.rho[otherId];
            i_delta = dist;
            i_nearestHigher = otherId;
          } else if (foundHigher && dist == i_delta && cellsOnLayer.rho[otherId] > rho_max) {
            rho_max = cellsOnLayer.rho[otherId];
            i_delta = dist;
            i_nearestHigher = otherId;
          } else if (foundHigher && dist == i_delta && cellsOnLayer.rho[otherId] == rho_max && cellsOnLayer.detid[otherId] > cellsOnLayer.detid[i]) {
            rho_max = cellsOnLayer.rho[otherId];
            i_delta = dist;
            i_nearestHigher = otherId;
          }
        }
      }
    }
    bool foundNearestHigherInSearchBox = (i_delta != maxDelta);
    if (foundNearestHigherInSearchBox) {
      cellsOnLayer.delta[i] = std::sqrt(i_delta);
      cellsOnLayer.nearestHigher[i] = i_nearestHigher;
    } else {
      // otherwise delta is guaranteed to be larger outlierDeltaFactor_*delta_c
      // we can safely maximize delta to be maxDelta
      cellsOnLayer.delta[i] = maxDelta;
      cellsOnLayer.nearestHigher[i] = -1;
    }

    LogDebug("HGCalCLUEAlgo") << "Debugging calculateDistanceToHigher: \n"
                              << "  cell: " << i << " eta: " << cellsOnLayer.dim1[i] << " phi: " << cellsOnLayer.dim2[i]
                              << " energy: " << cellsOnLayer.weight[i] << " density: " << cellsOnLayer.rho[i]
                              << " nearest higher: " << cellsOnLayer.nearestHigher[i]
                              << " distance: " << cellsOnLayer.delta[i] << "\n";
  }
}

template <typename T, typename STRATEGY>
int HGCalCLUEAlgoT<T, STRATEGY>::findAndAssignClusters(const unsigned int layerId, float delta) {
  // this is called once per layer and endcap...
  // so when filling the cluster temporary vector of Hexels we resize each time
  // by the number  of clusters found. This is always equal to the number of
  // cluster centers...
  unsigned int nClustersOnLayer = 0;
  auto& cellsOnLayer = cells_->at(layerId);
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  std::vector<int> localStack;
  // find cluster seeds and outlier
  for (unsigned int i = 0; i < numberOfCells; i++) {
    float rho_c = kappa_ * cellsOnLayer.sigmaNoise[i];
    // initialize clusterIndex
    cellsOnLayer.clusterIndex[i] = -1;
    bool isSeed = (cellsOnLayer.delta[i] > delta) && (cellsOnLayer.rho[i] >= rho_c);
    bool isOutlier = (cellsOnLayer.delta[i] > outlierDeltaFactor_ * delta) && (cellsOnLayer.rho[i] < rho_c);
    if (isSeed) {
      cellsOnLayer.clusterIndex[i] = nClustersOnLayer;
      cellsOnLayer.isSeed[i] = true;
      nClustersOnLayer++;
      localStack.push_back(i);

    } else if (!isOutlier) {
      cellsOnLayer.followers[cellsOnLayer.nearestHigher[i]].push_back(i);
    }
  }

  // need to pass clusterIndex to their followers
  while (!localStack.empty()) {
    int endStack = localStack.back();
    auto& thisSeed = cellsOnLayer.followers[endStack];
    localStack.pop_back();

    // loop over followers
    for (int j : thisSeed) {
      // pass id to a follower
      cellsOnLayer.clusterIndex[j] = cellsOnLayer.clusterIndex[endStack];
      // push this follower to localStack
      localStack.push_back(j);
    }
  }
  return nClustersOnLayer;
}

// explicit template instantiation
template class HGCalCLUEAlgoT<HGCalSiliconLayerTiles, HGCalSiliconStrategy>;
template class HGCalCLUEAlgoT<HGCalScintillatorLayerTiles, HGCalScintillatorStrategy>;
template class HGCalCLUEAlgoT<HFNoseLayerTiles, HGCalSiliconStrategy>;
