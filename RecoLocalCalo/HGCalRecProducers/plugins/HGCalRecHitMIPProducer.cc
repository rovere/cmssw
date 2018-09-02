/** \class HGCalRecHitMIPProducer
 *   produce HGCAL rechits (MIP-like)from calibrated rechits
 *
 *  \author Marco Rovere
 *
 **/
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include <iostream>

class HGCalRecHitMIPProducer : public edm::stream::EDProducer<> {

 public:
  explicit HGCalRecHitMIPProducer(const edm::ParameterSet& ps);
  ~HGCalRecHitMIPProducer() override;
  void produce(edm::Event& evt, const edm::EventSetup& es) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:

  void computeMipEnergies();
  const edm::EDGetTokenT<HGCRecHitCollection> eeRecHitCollection_;
  const edm::EDGetTokenT<HGCRecHitCollection> hefRecHitCollection_;
  const edm::EDGetTokenT<HGCRecHitCollection> hebRecHitCollection_;
  const std::string ee_mip_RechitCollection_;
  const std::string hef_mip_RechitCollection_;
  const std::string heb_mip_RechitCollection_;
  std::vector<double> weights_;
  std::vector<double> thickness_corrections_;
  std::vector<double> cce_;
  unsigned maxlayers_;
  unsigned thicknesses_;
  double mip_cut_;
  std::vector<std::vector<double> > mip_energy_gev_; // first index is layer number, second is thickness
  hgcal::RecHitTools rhtools_;
};

HGCalRecHitMIPProducer::HGCalRecHitMIPProducer(const edm::ParameterSet& ps) :
  eeRecHitCollection_( consumes<HGCRecHitCollection>( ps.getParameter<edm::InputTag>("HGCEERecHitCollection") ) ),
  hefRecHitCollection_( consumes<HGCRecHitCollection>( ps.getParameter<edm::InputTag>("HGCHEFRecHitCollection") ) ),
  hebRecHitCollection_( consumes<HGCRecHitCollection>( ps.getParameter<edm::InputTag>("HGCHEBRecHitCollection") ) ),
  ee_mip_RechitCollection_( ps.getUntrackedParameter<std::string>("HGCEEMIPrechitCollection") ),
  hef_mip_RechitCollection_( ps.getUntrackedParameter<std::string>("HGCHEFMIPrechitCollection") ),
  heb_mip_RechitCollection_( ps.getUntrackedParameter<std::string>("HGCHEBMIPrechitCollection") ),
  weights_( ps.getParameter<std::vector<double> >("weights")),
  thickness_corrections_( ps.getParameter<std::vector<double> >("thickness_corrections")),
  cce_ (ps.getParameter<std::vector<double> >("cce")),
  maxlayers_(ps.getUntrackedParameter<unsigned>("maxlayers")),
  thicknesses_(ps.getUntrackedParameter<unsigned>("thicknesses")),
  mip_cut_(ps.getUntrackedParameter<double>("mip_cut")) {
  produces< HGCRecHitCollection >(ee_mip_RechitCollection_);
  produces< HGCRecHitCollection >(hef_mip_RechitCollection_);
  produces< HGCRecHitCollection >(heb_mip_RechitCollection_);
  computeMipEnergies();
}

HGCalRecHitMIPProducer::~HGCalRecHitMIPProducer() {
}

void
HGCalRecHitMIPProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  using namespace edm;

  Handle< HGCRecHitCollection > HGCeeRecHits_handle;
  Handle< HGCRecHitCollection > HGChefRecHits_handle;
  Handle< HGCRecHitCollection > HGChebRecHits_handle;

  const HGCRecHitCollection*  eeRecHits = nullptr;
  const HGCRecHitCollection*  hefRecHits = nullptr;
  const HGCRecHitCollection*  hebRecHits = nullptr;

  // Initialize the RecHitTool
  rhtools_.getEventSetup(es);

  // get the HGC uncalib rechit collection
  evt.getByToken( eeRecHitCollection_, HGCeeRecHits_handle);
  eeRecHits = HGCeeRecHits_handle.product();

  evt.getByToken( hefRecHitCollection_, HGChefRecHits_handle);
  hefRecHits = HGChefRecHits_handle.product();

  evt.getByToken( hebRecHitCollection_, HGChebRecHits_handle);
  hebRecHits = HGChebRecHits_handle.product();

  auto mip_selection = [&](const HGCRecHit & hit)->bool {
    DetId detid = hit.detid();
    int layer = rhtools_.getLayerWithOffset(detid);
    int thickness_idx = rhtools_.getSiThickIndex(detid);
    return (hit.signalOverSigmaNoise() > 3.0f
        && hit.energy() > mip_cut_*mip_energy_gev_[layer][thickness_idx]);
  };
  // Loop over RecHits and filter them

  // collection of rechits to put in the event
  auto ee_mipRecHits = std::make_unique<HGCRecHitCollection>();
  ee_mipRecHits->reserve(eeRecHits->size());
  std::copy_if(eeRecHits->begin(), eeRecHits->end(), ee_mipRecHits->begin(), mip_selection);
//  auto last_element = std::copy_if(eeRecHits->begin(), eeRecHits->end(), ee_mipRecHits->begin(), mip_selection);
//  ee_mipRecHits->resize(std::distance(ee_mipRecHits->begin(), last_element));

  auto hef_mipRecHits = std::make_unique<HGChefRecHitCollection>();
  auto heb_mipRecHits = std::make_unique<HGChebRecHitCollection>();


  // put the collection of recunstructed hits in the event
  std::cout << "total # HGCee calibrated rechits: " << eeRecHits->size() << std::endl;
  std::cout << "total # HGChef calibrated rechits: " << hefRecHits->size() << std::endl;
  std::cout << "total # HGCheb calibrated rechits: " << hebRecHits->size() << std::endl;
  std::cout << "total # HGCeeMip calibrated rechits: " << ee_mipRecHits->size() << std::endl;
  std::cout << "total # HGChefMip calibrated rechits: " << hef_mipRecHits->size() << std::endl;
  std::cout << "total # HGChebMip calibrated rechits: " << heb_mipRecHits->size() << std::endl;

  evt.put(std::move(ee_mipRecHits), ee_mip_RechitCollection_);
  evt.put(std::move(hef_mipRecHits), hef_mip_RechitCollection_);
  evt.put(std::move(heb_mipRecHits), heb_mip_RechitCollection_);
}

void HGCalRecHitMIPProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("HGCEERecHitCollection", edm::InputTag("HGCalRecHit:HGCEERecHits"));
  desc.add<edm::InputTag>("HGCHEFRecHitCollection", edm::InputTag("HGCalRecHit:HGCHEFRecHits"));
  desc.add<edm::InputTag>("HGCHEBRecHitCollection", edm::InputTag("HGCalRecHit:HGCHEBRecHits"));
  desc.addUntracked<std::string>("HGCEEMIPrechitCollection", "HGCEEMipRecHits");
  desc.addUntracked<std::string>("HGCHEFMIPrechitCollection", "HGCHEFMipRecHits");
  desc.addUntracked<std::string>("HGCHEBMIPrechitCollection", "HGCHEBMipRecHits");
  desc.add<std::vector<double> >("weights", {0.});
  desc.add<std::vector<double> >("thickness_corrections", {0.});
  desc.add<std::vector<double> >("cce", {1.0, 1.0, 1.0});
  desc.addUntracked<unsigned>("maxlayers", 52);
  desc.addUntracked<unsigned>("thicknesses", 3);
  desc.addUntracked<double>("mip_cut", 3.0);
  descriptions.add("HGCalMipLikeRecHit",desc);
}


void HGCalRecHitMIPProducer::computeMipEnergies() {

  assert(weights_.size() >= maxlayers_);
  assert(thickness_corrections_.size() >= thicknesses_);
  assert(cce_.size() >= thicknesses_);
  std::vector<double> dummy;
  dummy.resize(thicknesses_, 0);
  mip_energy_gev_.resize(maxlayers_, dummy);


  for (unsigned l = 0; l < maxlayers_; ++l) {
    for (unsigned t = 0; t < thicknesses_; ++t) {
      mip_energy_gev_[l][t] = weights_[1+l] * thickness_corrections_[t] * 0.001f/cce_[t];
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( HGCalRecHitMIPProducer );
