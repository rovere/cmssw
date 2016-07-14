
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
//    std::cout << "CA_layerPairId: " << layerPairId
//              << " assigned to: " << aDoublet.first.first << " , "
//              << aDoublet.first.second << std::endl;
    ca_connetions_.push_back(CAConnections(layerPairId, aDoublet.first.first, aDoublet.first.second));
    std::tie(innerLayerId, outerLayerId) = aDoublet.first;
    auto numberOfDoublets = aDoublet.second.size();
    auto const & doublets = aDoublet.second;
//    std::cout << "Resizing isOuterHitOfCell on layer: " << outerLayerId
//              << " to: " << hits_on_layer.at(outerLayerId) << std::endl;
    isOuterHitOfCell[outerLayerId].resize(hits_on_layer.at(outerLayerId));
//    std::cout << "Filling layer: " << layerPairId << " out of " << layersMap.size() << std::endl;
    theFoundCellsPerLayer[layerPairId].reserve(numberOfDoublets);
    for (unsigned int i = 0; i < numberOfDoublets; ++i) {
//      std::cout << "theFoundCellsPerLayer for "
//                << layerPairId << " isOuterHitOfCell for " << outerLayerId << std::endl;
      theFoundCellsPerLayer[layerPairId].emplace_back(
          &doublets, i, /* cellId,*/ doublets.innerHitId(i),
          doublets.outerHitId(i));
//      std::cout << "Filling isOuterHitOfCell: " << outerLayerId
//                << ", " << doublets.outerHitId(i) << std::endl;
      isOuterHitOfCell[outerLayerId][doublets.outerHitId(i)]
          .push_back(&(theFoundCellsPerLayer[layerPairId][i]));
      if (innerLayerId > 0)
        for (auto neigCell : isOuterHitOfCell[innerLayerId][doublets.innerHitId(i)]) {
//          std::cout << "neigCell : theFoundCellsPerLayer for "
//                    << layerPairId << std::endl;
          theFoundCellsPerLayer[layerPairId][i].check_alignment_and_tag(
              neigCell, ptmin, region_origin_x, region_origin_y,
              region_origin_radius, thetaCut, phiCut);
        }
    }
    layerPairId++;
  }
//  dump("create_and_connect_cells");
}

template <unsigned int numberOfLayers>
void CellularAutomaton<numberOfLayers>::dump(const char * label) {
//  int num_layer = 0;
//  for (auto layer : theFoundCellsPerLayer) {
//    std::cout << label << ": Layer " << num_layer << " has size: " << layer.size() << std::endl;
//    for (auto const & cell : layer) {
//      std::cout << label << ": Layer: " << num_layer << " cell: " << cell.get_cell_id()
//                << "( " << cell.get_inner_hit_id() << ", " << cell.get_outer_hit_id()
//                << " ) ";
//      cell.dump();
//    }
//    num_layer++;
//  }
}

template <unsigned int numberOfLayers>
void CellularAutomaton<numberOfLayers>::evolve() {
//  std::vector<int> layersToEvolve = {0, 3, 4, 5, 1, 2, 9, 11};
  std::vector<int> layersToEvolve = {0, 1, 2, 3, 4, 5, 9, 10};
//  std::vector<int> layersToEvolve = {0, 3}; //Works OK
//  constexpr unsigned int numberOfIterations = numberOfLayers - 1;
  unsigned int the_numberOfIterations = 1;
  unsigned int numberOfCellsFound;
  for (unsigned int iteration = 0; iteration < the_numberOfIterations;
       ++iteration) {
//    std::cout << "Evolve at iteration: " << iteration << " of " << the_numberOfIterations << std::endl;
    for (auto const innerLayerId : layersToEvolve) {
//      std::cout << "Evolving cells of CA_layerPairId: " << innerLayerId << std::endl;
      for (auto& cell : theFoundCellsPerLayer[innerLayerId]) {
        cell.evolve();
      }
    }

//   std::cout << "Updating status of cells of CA_layerPairId: " << innerLayerId << std::endl; for (unsigned int innerLayerId = 0;
//         innerLayerId < theFoundCellsPerLayer.size() - iteration - 1; ++innerLayerId) {
    for (auto const innerLayerId : layersToEvolve) {
//      std::cout << "Updating status of cells of CA_layerPairId: " << innerLayerId << std::endl;
      for (auto& cell : theFoundCellsPerLayer[innerLayerId]) {
        cell.update_state();
      }
    }
  }

  //last iteration
  std::vector<int> initialLayers = {0, 1, 2};
  for (auto const layer : initialLayers) {
    numberOfCellsFound = theFoundCellsPerLayer[layer].size();

//  std::cout << "Number of found cells before last iteration: " << numberOfCellsFound << std::endl;
    for (unsigned int cellId = 0; cellId < numberOfCellsFound; ++cellId) {
      theFoundCellsPerLayer[layer][cellId].evolve();
    }

    for (auto& cell : theFoundCellsPerLayer[layer]) {
      cell.update_state();
//    std::cout << "Target[2]; comparing cell: " << cell.get_cell_id()
//              << " with state:  " << cell.get_CA_state() << std::endl;

      if (cell.is_root_cell(2)) {
        theRootCells.push_back(&cell);
      }
    }
  }
//  std::cout << "I have found " << theRootCells.size()
//            << " final cells!" << std::endl;
//  for (auto const & c : theRootCells) {
//    std::cout << "Final Root Cell: ";
//    c->dump();
//  }

//  dump("evolve");
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

//  int v_num = 0;
//  for (auto const & v: foundNtuplets) {
//    int vv_num = 0;
//    for (auto const & vv : v) {
//      std::cout << "find_ntuplets: " << v_num << ", " << vv_num;
//      vv->dump();
//      vv_num++;
//    }
//    v_num++;
//  }
}

template class CellularAutomaton<13>;
