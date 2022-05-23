#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "vdt/vdtMath.h"

#include <iostream>

using namespace hgcal;
ClusterTools::ClusterTools() {}

ClusterTools::ClusterTools(const edm::ParameterSet& conf, edm::ConsumesCollector& sumes)
    : eetok(sumes.consumes<HGCRecHitCollection>(conf.getParameter<edm::InputTag>("HGCEEInput"))),
      fhtok(sumes.consumes<HGCRecHitCollection>(conf.getParameter<edm::InputTag>("HGCFHInput"))),
      bhtok(sumes.consumes<HGCRecHitCollection>(conf.getParameter<edm::InputTag>("HGCBHInput"))),
      caloGeometryToken_{sumes.esConsumes()} {}

void ClusterTools::getEvent(const edm::Event& ev) {
  eerh_ = &ev.get(eetok);
  fhrh_ = &ev.get(fhtok);
  bhrh_ = &ev.get(bhtok);
}

void ClusterTools::getEventSetup(const edm::EventSetup& es) { rhtools_.setGeometry(es.getData(caloGeometryToken_)); }

float ClusterTools::getClusterHadronFraction(const reco::CaloCluster& clus) const {
  float energy = 0.f, energyHad = 0.f;
  const auto& hits = clus.hitsAndFractions();
  for (const auto& hit : hits) {
    const auto& id = hit.first;
    const float fraction = hit.second;
    if (id.det() == DetId::HGCalEE) {
      energy += eerh_->find(id)->energy() * fraction;
    } else if (id.det() == DetId::HGCalHSi) {
      const float temp = fhrh_->find(id)->energy();
      energy += temp * fraction;
      energyHad += temp * fraction;
    } else if (id.det() == DetId::HGCalHSc) {
      const float temp = bhrh_->find(id)->energy();
      energy += temp * fraction;
      energyHad += temp * fraction;
    } else if (id.det() == DetId::Forward) {
      switch (id.subdetId()) {
        case HGCEE:
          energy += eerh_->find(id)->energy() * fraction;
          break;
        case HGCHEF: {
          const float temp = fhrh_->find(id)->energy();
          energy += temp * fraction;
          energyHad += temp * fraction;
        } break;
        default:
          throw cms::Exception("HGCalClusterTools") << " Cluster contains hits that are not from HGCal! " << std::endl;
      }
    } else if (id.det() == DetId::Hcal && id.subdetId() == HcalEndcap) {
      const float temp = bhrh_->find(id)->energy();
      energy += temp * fraction;
      energyHad += temp * fraction;
    } else {
      throw cms::Exception("HGCalClusterTools") << " Cluster contains hits that are not from HGCal! " << std::endl;
    }
  }
  float fraction = -1.f;
  if (energy > 0.f) {
    fraction = energyHad / energy;
  }
  return fraction;
}

int ClusterTools::getLayer(const DetId detid) const { return rhtools_.getLayerWithOffset(detid); }

bool ClusterTools::getWidths(const reco::CaloCluster& clus,
                             double& sigmaetaeta,
                             double& sigmaphiphi,
                             double& sigmaetaetal,
                             double& sigmaphiphil) const {
  if (getLayer(clus.hitsAndFractions()[0].first) > (int)rhtools_.lastLayerEE())
    return false;
  const math::XYZPoint& position(clus.position());
  unsigned nhit = clus.hitsAndFractions().size();

  sigmaetaeta = 0.;
  sigmaphiphi = 0.;
  sigmaetaetal = 0.;
  sigmaphiphil = 0.;

  double sumw = 0.;
  double sumlogw = 0.;

  for (unsigned int ih = 0; ih < nhit; ++ih) {
    const DetId& id = (clus.hitsAndFractions())[ih].first;
    if ((clus.hitsAndFractions())[ih].second == 0.)
      continue;

    if ((id.det() == DetId::HGCalEE) || (id.det() == DetId::Forward && id.subdetId() == HGCEE)) {
      const HGCRecHit* theHit = &(*eerh_->find(id));

      GlobalPoint cellPos = rhtools_.getPosition(id);
      double weight = theHit->energy();
      // take w0=2 To be optimized
      double logweight = 0;
      if (clus.energy() != 0) {
        logweight = std::max(0., 2 + log(theHit->energy() / clus.energy()));
      }
      double deltaetaeta2 = (cellPos.eta() - position.eta()) * (cellPos.eta() - position.eta());
      double deltaphiphi2 = (cellPos.phi() - position.phi()) * (cellPos.phi() - position.phi());
      sigmaetaeta += deltaetaeta2 * weight;
      sigmaphiphi += deltaphiphi2 * weight;
      sigmaetaetal += deltaetaeta2 * logweight;
      sigmaphiphil += deltaphiphi2 * logweight;
      sumw += weight;
      sumlogw += logweight;
    }
  }

  if (sumw <= 0.)
    return false;

  sigmaetaeta /= sumw;
  sigmaetaeta = std::sqrt(sigmaetaeta);
  sigmaphiphi /= sumw;
  sigmaphiphi = std::sqrt(sigmaphiphi);

  if (sumlogw != 0) {
    sigmaetaetal /= sumlogw;
    sigmaetaetal = std::sqrt(sigmaetaetal);
    sigmaphiphil /= sumlogw;
    sigmaphiphil = std::sqrt(sigmaphiphil);
  }

  return true;
}
