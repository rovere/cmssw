
#include "CellularAutomaton.h"

template <unsigned int numberOfLayers>
void CellularAutomaton<numberOfLayers>::create_and_connect_cells(
    const std::map<std::pair<int, int>, HitDoublets> &layersMap,
    const std::unordered_map<int, int>& hits_on_layer,
    const TrackingRegion& region, const float thetaCut, const float phiCut) {
  float ptmin = region.ptMin();
  float region_origin_x = region.origin().x();
  float region_origin_y = region.origin().y();
  float region_origin_radius = region.originRBound();
  unsigned int layerPairId = 0;
  int innerLayerId = 0;
  int outerLayerId = 0;
  // for (auto const & lh : hits_on_layer)
  //   std::cout << "Layer: " << lh.first << " has hits: " << lh.second << std::endl;
  for (auto const & aDoublet : layersMap) {
    std::tie(innerLayerId, outerLayerId) = aDoublet.first;
    auto numberOfDoublets = aDoublet.second.size();
    auto const & doublets = aDoublet.second;
    // std::cout << "Resizing isOuterHitOfCell on layer: " << outerLayerId
    //           << " to: " << hits_on_layer.at(outerLayerId) << std::endl;
    isOuterHitOfCell[outerLayerId].resize(hits_on_layer.at(outerLayerId));
    // std::cout << "Filling layer: " << layerPairId << " out of " << layersMap.size() << std::endl;
    theFoundCellsPerLayer[layerPairId].reserve(numberOfDoublets);
    for (unsigned int i = 0; i < numberOfDoublets; ++i) {
      theFoundCellsPerLayer[layerPairId].emplace_back(
          &doublets, i, /* cellId,*/ doublets.innerHitId(i),
          doublets.outerHitId(i));
      // std::cout << "Filling isOuterHitOfCell: " << outerLayerId
      //           << ", " << doublets.outerHitId(i) << std::endl;
      isOuterHitOfCell[outerLayerId][doublets.outerHitId(i)]
          .push_back(&(theFoundCellsPerLayer[layerPairId][i]));
      if (innerLayerId > 0)
        for (auto neigCell : isOuterHitOfCell[innerLayerId][doublets.innerHitId(i)]) {
          theFoundCellsPerLayer[layerPairId][i].check_alignment_and_tag(
              neigCell, ptmin, region_origin_x, region_origin_y,
              region_origin_radius, thetaCut, phiCut);
        }
    }
    layerPairId++;
  }
  dump("create_and_connect_cells");
}

template <unsigned int numberOfLayers>
void CellularAutomaton<numberOfLayers>::dump(const char * label) {
  int num_layer = 0;
  for (auto layer : theFoundCellsPerLayer) {
    std::cout << label << ": Layer " << num_layer << " has size: " << layer.size() << std::endl;
    for (auto const & cell : layer) {
      std::cout << label << ": Layer: " << num_layer << " cell: " << cell.get_cell_id()
                << "( " << cell.get_inner_hit_id() << ", " << cell.get_outer_hit_id()
                << " ) ";
      cell.dump();
    }
    num_layer++;
  }
}

template <unsigned int numberOfLayers>
void CellularAutomaton<numberOfLayers>::evolve() {
  constexpr unsigned int numberOfIterations = numberOfLayers - 2;
  unsigned int numberOfCellsFound;
  for (unsigned int iteration = 0; iteration < numberOfIterations - 1;
       ++iteration) {
    for (unsigned int innerLayerId = 0;
         innerLayerId < numberOfIterations - iteration; ++innerLayerId) {
      for (auto& cell : theFoundCellsPerLayer[innerLayerId]) {
        cell.evolve();
      }
    }

    for (unsigned int innerLayerId = 0;
         innerLayerId < numberOfLayers - iteration - 2; ++innerLayerId) {
      for (auto& cell : theFoundCellsPerLayer[innerLayerId]) {
        cell.update_state();
      }
    }
  }

  //last iteration
  numberOfCellsFound = theFoundCellsPerLayer[0].size();

  for (unsigned int cellId = 0; cellId < numberOfCellsFound; ++cellId) {
    theFoundCellsPerLayer[0][cellId].evolve();
  }

  for (auto& cell : theFoundCellsPerLayer[0]) {
    cell.update_state();
    if (cell.is_root_cell(numberOfLayers - 2)) {
      theRootCells.push_back(&cell);
    }
  }
  dump("evolve");
}

template <unsigned int numberOfLayers>
void CellularAutomaton<numberOfLayers>::find_ntuplets(
    std::vector<CACell::CAntuplet>& foundNtuplets,
    const unsigned int minHitsPerNtuplet) {
  std::vector<CACell*> tmpNtuplet;
  tmpNtuplet.reserve(numberOfLayers);

  for (CACell* root_cell : theRootCells) {
    tmpNtuplet.clear();
    tmpNtuplet.push_back(root_cell);
    root_cell->find_ntuplets(foundNtuplets, tmpNtuplet, minHitsPerNtuplet);
  }

  int v_num = 0;
  for (auto const & v: foundNtuplets) {
    int vv_num = 0;
    for (auto const & vv : v) {
      std::cout << "find_ntuplets: " << v_num << ", " << vv_num;
      vv->dump();
      vv_num++;
    }
    v_num++;
  }
}

template class CellularAutomaton<13>;
