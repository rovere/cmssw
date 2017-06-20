// user include files
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <map>
#include <string>

class HGCalHitCalibration : public DQMEDAnalyzer {
 public:
  explicit HGCalHitCalibration(const edm::ParameterSet&);
  ~HGCalHitCalibration();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&,
                      edm::EventSetup const&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  void fillWithRecHits(std::map<DetId, const HGCRecHit*>&, DetId, unsigned int,
                       float, int&, float&);

  edm::EDGetTokenT<HGCRecHitCollection> recHitsEE_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsFH_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsBH_;
  edm::EDGetTokenT<std::vector<CaloParticle> > caloParticles_;
  edm::EDGetTokenT<std::vector<reco::PFCluster> > hgcalMultiClusters_;

  std::string detector_;
  int algo_;
  HGCalDepthPreClusterer pre_;
  bool rawRecHits_;
  hgcal::RecHitTools recHitTools_;

  std::map<int, MonitorElement*> h_EoP_CPene_calib_fraction_;
  std::map<int, MonitorElement*> hgcal_EoP_CPene_calib_fraction_;
  MonitorElement* LayerOccupancy_;

  std::vector<float> Energy_layer_calib_;
  std::vector<float> Energy_layer_calib_fraction_;
};

HGCalHitCalibration::HGCalHitCalibration(const edm::ParameterSet& iConfig)
    : detector_(iConfig.getParameter<std::string>("detector")),
      rawRecHits_(iConfig.getParameter<bool>("rawRecHits")) {
  if (detector_ == "all") {
    recHitsEE_ = consumes<HGCRecHitCollection>(
        edm::InputTag("HGCalRecHit", "HGCEERecHits"));
    recHitsFH_ = consumes<HGCRecHitCollection>(
        edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
    recHitsBH_ = consumes<HGCRecHitCollection>(
        edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
    algo_ = 1;
  } else if (detector_ == "EM") {
    recHitsEE_ = consumes<HGCRecHitCollection>(
        edm::InputTag("HGCalRecHit", "HGCEERecHits"));
    algo_ = 2;
  } else if (detector_ == "HAD") {
    recHitsFH_ = consumes<HGCRecHitCollection>(
        edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
    recHitsBH_ = consumes<HGCRecHitCollection>(
        edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
    algo_ = 3;
  }
  caloParticles_ = consumes<std::vector<CaloParticle> >(
      edm::InputTag("mix", "MergedCaloTruth"));
  hgcalMultiClusters_ = consumes<std::vector<reco::PFCluster> >(
      edm::InputTag("particleFlowClusterHGCalFromMC"));
}

HGCalHitCalibration::~HGCalHitCalibration() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

void HGCalHitCalibration::bookHistograms(DQMStore::IBooker& ibooker,
                                         edm::Run const& iRun,
                                         edm::EventSetup const& /* iSetup */) {
  ibooker.cd();
  ibooker.setCurrentFolder("HGCalHitCalibration");
  h_EoP_CPene_calib_fraction_[100] =
      ibooker.book1D("h_EoP_CPene_100_calib_fraction", "", 1000, -0.5, 2.5);
  h_EoP_CPene_calib_fraction_[200] =
      ibooker.book1D("h_EoP_CPene_200_calib_fraction", "", 1000, -0.5, 2.5);
  h_EoP_CPene_calib_fraction_[300] =
      ibooker.book1D("h_EoP_CPene_300_calib_fraction", "", 1000, -0.5, 2.5);
  hgcal_EoP_CPene_calib_fraction_[100] =
      ibooker.book1D("hgcal_EoP_CPene_100_calib_fraction", "", 1000, -0.5, 2.5);
  hgcal_EoP_CPene_calib_fraction_[200] =
      ibooker.book1D("hgcal_EoP_CPene_200_calib_fraction", "", 1000, -0.5, 2.5);
  hgcal_EoP_CPene_calib_fraction_[300] =
      ibooker.book1D("hgcal_EoP_CPene_300_calib_fraction", "", 1000, -0.5, 2.5);
  LayerOccupancy_ = ibooker.book1D("LayerOccupancy", "", 60, 0., 60.);
}

void HGCalHitCalibration::fillWithRecHits(
    std::map<DetId, const HGCRecHit*>& hitmap, const DetId hitid,
    const unsigned int hitlayer, const float fraction, int& seedDet,
    float& seedEnergy) {
  if (hitmap.find(hitid) == hitmap.end()) {
    // Hit was not reconstructed
    return;
  }
  unsigned int layer = recHitTools_.getLayerWithOffset(hitid);
  assert(hitlayer == layer);
  Energy_layer_calib_fraction_[layer] += hitmap[hitid]->energy() * fraction;
  LayerOccupancy_->Fill(layer);
  if (seedEnergy < hitmap[hitid]->energy()) {
    seedEnergy = hitmap[hitid]->energy();
    seedDet = recHitTools_.getSiThickness(hitid);
  }
}

void HGCalHitCalibration::analyze(const edm::Event& iEvent,
                                  const edm::EventSetup& iSetup) {
  using namespace edm;

  // size should be HGC layers 52 is enough
  Energy_layer_calib_.clear();
  Energy_layer_calib_fraction_.clear();
  for (unsigned int ij = 0; ij < 60; ++ij) {
    Energy_layer_calib_.push_back(0.);
    Energy_layer_calib_fraction_.push_back(0.);
  }

  recHitTools_.getEventSetup(iSetup);

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;

  Handle<std::vector<CaloParticle> > caloParticleHandle;
  iEvent.getByToken(caloParticles_, caloParticleHandle);
  const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;

  Handle<std::vector<reco::PFCluster> > hgcalMultiClustersHandle;
  iEvent.getByToken(hgcalMultiClusters_, hgcalMultiClustersHandle);

  // make a map detid-rechit
  std::map<DetId, const HGCRecHit*> hitmap;
  switch (algo_) {
    case 1: {
      iEvent.getByToken(recHitsEE_, recHitHandleEE);
      iEvent.getByToken(recHitsFH_, recHitHandleFH);
      iEvent.getByToken(recHitsBH_, recHitHandleBH);
      const auto& rechitsEE = *recHitHandleEE;
      const auto& rechitsFH = *recHitHandleFH;
      const auto& rechitsBH = *recHitHandleBH;
      for (unsigned int i = 0; i < rechitsEE.size(); ++i) {
        hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      for (unsigned int i = 0; i < rechitsFH.size(); ++i) {
        hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for (unsigned int i = 0; i < rechitsBH.size(); ++i) {
        hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
    case 2: {
      iEvent.getByToken(recHitsEE_, recHitHandleEE);
      const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
      for (unsigned int i = 0; i < rechitsEE.size(); i++) {
        hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      break;
    }
    case 3: {
      iEvent.getByToken(recHitsFH_, recHitHandleFH);
      iEvent.getByToken(recHitsBH_, recHitHandleBH);
      const auto& rechitsFH = *recHitHandleFH;
      const auto& rechitsBH = *recHitHandleBH;
      for (unsigned int i = 0; i < rechitsFH.size(); i++) {
        hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for (unsigned int i = 0; i < rechitsBH.size(); i++) {
        hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
    default:
      break;
  }

  // loop over caloParticles
  int seedDet = 0;
  float seedEnergy = 0.;
  std::cout << "Number of caloParticles: " << caloParticles.size() << std::endl;
  for (const auto& it_caloPart : caloParticles) {
    const SimClusterRefVector& simClusterRefVector = it_caloPart.simClusters();
    std::cout << "Simclusters linked to a single CaloParticle have size: "
              << simClusterRefVector.size() << std::endl;
    // Bail out if the number of caloparticles is not equal to the
    // number of simClusters: that's the case, for example, of a
    // photon that converted before HGCal.  TODO (rovere): prepare a
    // 0-material scenarion for the PhaseII, so that these kind of
    // tricks won't be needed.
    if (caloParticles.size() != simClusterRefVector.size()) return;
    Energy_layer_calib_.clear();
    Energy_layer_calib_fraction_.clear();
    for (unsigned int ij = 0; ij < 60; ++ij) {
      Energy_layer_calib_.push_back(0.);
      Energy_layer_calib_fraction_.push_back(0.);
    }

    seedDet = 0;
    seedEnergy = 0.;
    for (const auto& it_sc : simClusterRefVector) {
      const SimCluster& simCluster = (*(it_sc));
      std::cout << ">>> simCluster.energy() = " << simCluster.energy()
                << std::endl;
      const std::vector<std::pair<uint32_t, float> >& hits_and_fractions =
          simCluster.hits_and_fractions();

      // loop over hits
      for (const auto& it_haf : hits_and_fractions) {
        unsigned int hitlayer = recHitTools_.getLayerWithOffset(it_haf.first);
        DetId hitid = (it_haf.first);
        // std::cout << "looking for " << hitid.rawId() << std::endl;
        // dump raw RecHits and match
        if (rawRecHits_) {
          if (hitid.det() == DetId::Forward &&
              (hitid.subdetId() == HGCEE or hitid.subdetId() == HGCHEF or
               hitid.subdetId() == HGCHEB))
            fillWithRecHits(hitmap, hitid, hitlayer, it_haf.second, seedDet,
                            seedEnergy);
        }
      }  // end simHit
    }    // end simCluster

    float sumCalibRecHitCalib_fraction = 0;
    for (unsigned int iL = 0; iL < Energy_layer_calib_fraction_.size(); ++iL) {
      sumCalibRecHitCalib_fraction += Energy_layer_calib_fraction_[iL];
    }

    if (h_EoP_CPene_calib_fraction_.find(seedDet) !=
        h_EoP_CPene_calib_fraction_.end())
      h_EoP_CPene_calib_fraction_[seedDet]->Fill(sumCalibRecHitCalib_fraction /
                                                 it_caloPart.energy());
    // TODO
    // * Loop over HGCalMultiClusters
    // * Select only the most energetic one
    // * Loop over its recHits_fractions
    // * Make the same plots as before.
    const auto& clusters = *hgcalMultiClustersHandle;
    float total_energy = 0.;
    for (const auto& c : clusters) {
      // no need to go through the hits_and_fractions gymnastic, since
      // the energy should be already correct. Also, in case of multiple
      // MultiClusters, consider only the most energetic one.
      if (c.correctedEnergy() > total_energy) {
        total_energy = c.correctedEnergy();
        seedDet = recHitTools_.getSiThickness(c.seed());
      }
    }
    if (hgcal_EoP_CPene_calib_fraction_.find(seedDet) !=
        hgcal_EoP_CPene_calib_fraction_.end())
      hgcal_EoP_CPene_calib_fraction_[seedDet]->Fill(total_energy /
                                                     it_caloPart.energy());
  }  // end caloparticle
}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void HGCalHitCalibration::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(HGCalHitCalibration);
