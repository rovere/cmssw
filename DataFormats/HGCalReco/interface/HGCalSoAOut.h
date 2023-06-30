#ifndef DataFormats_HGCalReco_interface_HGCalSoAOut_h
#define DataFormats_HGCalReco_interface_HGCalSoAOut_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

  // SoA layout with delta, rho, weight, nearestHigher, clusterIndex, layer, isSeed, and cellsCount fields
  GENERATE_SOA_LAYOUT(HGCalSoAOut,
                      // columns: one value per element
                      SOA_COLUMN(float, delta),
                      SOA_COLUMN(float, rho),
                      SOA_COLUMN(unsigned int, nearestHigher),
                      SOA_COLUMN(int, clusterIndex),
                      SOA_COLUMN(uint8_t, isSeed),
                      SOA_SCALAR(unsigned int, numberOfClustersScalar)
                      )

  using HGCalCellsOutSoA = HGCalSoAOut<>;

#endif
