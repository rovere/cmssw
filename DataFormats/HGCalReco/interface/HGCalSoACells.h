#ifndef DataFormats_HGCalReco_interface_HGCalSoACells_h
#define DataFormats_HGCalReco_interface_HGCalSoACells_h

#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

  // SoA layout with x, y, z, id fields
  GENERATE_SOA_LAYOUT(HGCalSoACells,
                      // columns: one value per element
                      SOA_COLUMN(float, dim1),
                      SOA_COLUMN(float, dim2),
                      SOA_COLUMN(int, layer),
                      SOA_COLUMN(float, weight),
                      SOA_COLUMN(float, sigmaNoise))

  using HGCalCellsSoA = HGCalSoACells<>;

#endif  // DataFormats_PortableTestObjects_interface_TestSoA_h