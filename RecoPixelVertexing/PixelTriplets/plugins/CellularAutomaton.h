#ifndef RECOPIXELVERTEXING_PIXELTRIPLETS_PLUGINS_CELLULARAUTOMATON_H_
#define RECOPIXELVERTEXING_PIXELTRIPLETS_PLUGINS_CELLULARAUTOMATON_H_
#include <array>
#include <unordered_map>
#include "CACell.h"
#include "TrackingTools/TransientTrackingRecHit/interface/SeedingLayerSetsHits.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"

typedef struct logic {
  logic(int layer, int inner, int outer): layer_index(layer), inner_index(inner), outer_index(outer){}
  int layer_index;
  int inner_index;
  int outer_index;
} CAConnections;

template<unsigned int theNumberOfLayers>
class CellularAutomaton {
 public:
  CellularAutomaton() {}

  void create_and_connect_cells(const std::map<std::pair<int, int>, HitDoublets> &,
                                const std::unordered_map<int, int>&,
                                const TrackingRegion&,
                                const float, const float);
  void evolve();
  void find_ntuplets(std::vector<CACell::CAntuplet>&, const unsigned int);
  void dump(const char *);

 private:
  //for each hit in each layer, store the pointers of the Cells of which it is outerHit
  std::array<std::vector<std::vector<CACell*> >, theNumberOfLayers> isOuterHitOfCell;
  std::array<std::vector<CACell>, theNumberOfLayers> theFoundCellsPerLayer;
  std::vector<CACell*> theRootCells;
  std::vector<std::vector<CACell*> > theNtuplets;
  std::array<std::string, theNumberOfLayers> layers_;
  std::vector<CAConnections> ca_connetions_;
};


#endif 






