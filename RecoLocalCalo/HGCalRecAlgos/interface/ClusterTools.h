#ifndef __RecoLocalCalo_HGCalRecAlgos_ClusterTools_h__
#define __RecoLocalCalo_HGCalRecAlgos_ClusterTools_h__

#include <array>
#include <cmath>
#include <numeric>

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

class HGCalGeometry;
class HGCalDDDConstants;
class DetId;

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

namespace hgcal {
  class ClusterTools {
  public:
    ClusterTools();
    ClusterTools(const edm::ParameterSet &, edm::ConsumesCollector &);
    ~ClusterTools() {}

    void getEvent(const edm::Event &);
    void getEventSetup(const edm::EventSetup &);

    float getClusterHadronFraction(const reco::CaloCluster &) const;

    int getLayer(const DetId) const;

    // only for EE
    bool getWidths(const reco::CaloCluster &clus,
                   double &sigmaetaeta,
                   double &sigmaphiphi,
                   double &sigmaetaetalog,
                   double &sigmaphiphilog) const;

  private:
    RecHitTools rhtools_;
    const edm::EDGetTokenT<HGCRecHitCollection> eetok, fhtok, bhtok;
    const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;

    const HGCRecHitCollection *eerh_, *fhrh_, *bhrh_;
  };
}  // namespace hgcal

#endif
