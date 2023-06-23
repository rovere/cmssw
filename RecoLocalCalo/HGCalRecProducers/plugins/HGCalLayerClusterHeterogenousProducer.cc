// Authors: Olivie Franklova - olivie.abigail.franklova@cern.ch
// Date: 03/2023
// @file create layer clusters
#ifndef __RecoLocalCalo_HGCRecProducers_HGCalLayerClusterHeterogenousProducer_H__
#define __RecoLocalCalo_HGCRecProducers_HGCalLayerClusterHeterogenousProducer_H__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "RecoParticleFlow/PFClusterProducer/interface/RecHitTopologicalCleanerBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/SeedFinderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterBuilderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalLayerClusterAlgoFactory.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "CLUEAlgo.h"
#include "TilesConstants.h"

class HGCalLayerClusterHeterogenousProducer : public edm::stream::EDProducer<> {
public:
  /**
   * @brief Constructor with parameter settings - which can be changed in hgcalLayerCluster_cff.py.
   * Constructor will set all variables by input param ps. 
   * algoID variables will be set accordingly to the detector type.
   * 
   * @param[in] ps parametr set to set variables
  */
  HGCalLayerClusterHeterogenousProducer(const edm::ParameterSet&);
  ~HGCalLayerClusterHeterogenousProducer() override {}
  /**
   * @brief Method fill description which will be used in pyhton file.
   * 
   * @param[out] description to be fill
  */
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  /**
   * @brief Method run the algoritm to get clusters.
   * 
   * @param[in, out] evt from get info and put result
   * @param[in] es to get event setup info
  */
  void produce(edm::Event&, const edm::EventSetup&) override;



private:
  edm::EDGetTokenT<HGCRecHitCollection> hits_token_;

  reco::CaloCluster::AlgoId algoId_;

  std::string detector_;

  std::string timeClname_;
  unsigned int hitsTime_;

  // for calculate position
  std::vector<double> thresholdW0_;
  double positionDeltaRho2_;

  bool dependSensor_;
  bool initialized_;
  unsigned maxNumberOfThickIndices_;
  std::vector<std::vector<double>> thresholds_;
  std::vector<std::vector<double>> v_sigmaNoise_;
  bool isNose_;
  unsigned int maxlayer_;
  double sciThicknessCorrection_;
  std::vector<double> fcPerMip_;
  double fcPerEle_;
  std::vector<double> nonAgedNoises_;
  double noiseMip_;
  double ecut_;
  std::vector<double> dEdXweights_;
  std::vector<double> thicknessCorrection_;
  int deltasi_index_regemfac_;
  CellsOnLayers cells_;
  // double dc_;

  hgcal::RecHitTools rhtools_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  float roundTo5(float value){
    return std::round(value * 100000.0) / 100000.0;
    // return std::round(value / 100000.0) * 100000.0;
  }

  std::vector<reco::BasicCluster> createClusters(int numberOfClusters){
    
    int numbeOfseeds = 0;
    int numberOfCells = cells_.detid.size();

    std::vector<reco::BasicCluster> clusters;
    clusters.resize(numberOfClusters);
    std::set<int> usedIdx;
    for (int i = 0; i < numberOfCells; i++){
        if (cells_.clusterIndex[i] == -1){
          continue;
        }
        int globalClusterIdx = cells_.clusterIndex[i];
        if (usedIdx.find(globalClusterIdx) != usedIdx.end()){
          clusters[globalClusterIdx].setEnergy(clusters[globalClusterIdx].energy() + cells_.weight[i]);
          clusters[globalClusterIdx].addHitAndFraction(cells_.detid[i], 1.f);
        }
        else{
          std::vector<std::pair<DetId, float>> thisCluster;
          thisCluster.emplace_back(cells_.detid[i], 1.f);
          math::XYZPoint position = math::XYZPoint(0.f, 0.f, 0.f);
          clusters[globalClusterIdx] = reco::BasicCluster(cells_.weight[i], position, reco::CaloID::DET_HGCAL_ENDCAP, std::move(thisCluster), algoId_ );
          usedIdx.emplace(globalClusterIdx);
        }
        if (cells_.isSeed[i]){
          clusters[globalClusterIdx].setSeed(cells_.detid[i]);
          numbeOfseeds ++;
        }
    }
    // for (unsigned i = 0; i < clusters.size(); ++i) {
    //   const reco::CaloCluster& sCl = clusters[i];
    //   std::cout << clusters[i].hitsAndFractions().size() << std::endl;
    //   std::cout << sCl.hitsAndFractions().size() <<  " " << clusters.size() << std::endl; 
    // }
    return clusters;
  }
void computeThreshold() {

  // To support the TDR geometry and also the post-TDR one (v9 onwards), we
  // need to change the logic of the vectors containing signal to noise and
  // thresholds. The first 3 indices will keep on addressing the different
  // thicknesses of the Silicon detectors in CE_E , the next 3 indices will address
  // the thicknesses of the Silicon detectors in CE_H, while the last one, number 6 (the
  // seventh) will address the Scintillators. This change will support both
  // geometries at the same time.

  if (initialized_)
    return;  // only need to calculate thresholds once

  initialized_ = true;

  std::vector<double> dummy;

  dummy.resize(maxNumberOfThickIndices_ + !isNose_, 0);  // +1 to accomodate for the Scintillators
  thresholds_.resize(maxlayer_, dummy);
  v_sigmaNoise_.resize(maxlayer_, dummy);

  for (unsigned ilayer = 1; ilayer <= maxlayer_; ++ilayer) {
    for (unsigned ithick = 0; ithick < maxNumberOfThickIndices_; ++ithick) {
      float sigmaNoise = 0.001f * fcPerEle_ * nonAgedNoises_[ithick] * dEdXweights_[ilayer] /
                         (fcPerMip_[ithick] * thicknessCorrection_[ithick]);
      thresholds_[ilayer - 1][ithick] = sigmaNoise * ecut_;
      v_sigmaNoise_[ilayer - 1][ithick] = sigmaNoise;
      LogDebug("HGCalCLUEAlgo") << "ilayer: " << ilayer << " nonAgedNoises: " << nonAgedNoises_[ithick]
                                << " fcPerEle: " << fcPerEle_ << " fcPerMip: " << fcPerMip_[ithick]
                                << " noiseMip: " << fcPerEle_ * nonAgedNoises_[ithick] / fcPerMip_[ithick]
                                << " sigmaNoise: " << sigmaNoise << "\n";
    }

    if (!isNose_) {
      float scintillators_sigmaNoise = 0.001f * noiseMip_ * dEdXweights_[ilayer] / sciThicknessCorrection_;
      thresholds_[ilayer - 1][maxNumberOfThickIndices_] = ecut_ * scintillators_sigmaNoise;
      v_sigmaNoise_[ilayer - 1][maxNumberOfThickIndices_] = scintillators_sigmaNoise;
      LogDebug("HGCalCLUEAlgo") << "ilayer: " << ilayer << " noiseMip: " << noiseMip_
                                << " scintillators_sigmaNoise: " << scintillators_sigmaNoise << "\n";
    }
  }
}

void populate(const HGCRecHitCollection& hits) {
  // loop over all hits and create the Hexel structure, skip energies below ecut
  if (dependSensor_) {
    // for each layer and wafer calculate the thresholds (sigmaNoise and energy)
    // once
    computeThreshold(); 
  }
  int index = 0;
  for (unsigned int i = 0; i < hits.size(); ++i) {
    const HGCRecHit& hgrh = hits[i];
    DetId detid = hgrh.detid();
    unsigned int layerOnSide = (rhtools_.getLayerWithOffset(detid) - 1);

    // set sigmaNoise default value 1 to use kappa value directly in case of
    // sensor-independent thresholds
    float sigmaNoise = 1.f;
    if (dependSensor_) {
      int thickness_index = rhtools_.getSiThickIndex(detid);
      if (thickness_index == -1){
        thickness_index = maxNumberOfThickIndices_;
      }
      double storedThreshold = thresholds_[layerOnSide][thickness_index];
      if (detid.det() == DetId::HGCalHSi || detid.subdetId() == HGCHEF) {
        storedThreshold = thresholds_[layerOnSide][thickness_index + deltasi_index_regemfac_];
      }
      sigmaNoise = v_sigmaNoise_[layerOnSide][thickness_index];

      if (hgrh.energy() < storedThreshold)
        continue;  // this sets the ZS threshold at ecut times the sigma noise
                   // for the sensor
    }

    if (!dependSensor_ && hgrh.energy() < ecut_)
      continue;

    const GlobalPoint position(rhtools_.getPosition(detid));
    int offset = ((rhtools_.zside(detid) + 1) >> 1) * maxlayer_;
    int layer = layerOnSide + offset;
    cells_.detid[index] = detid;
    // cells_.detid.emplace_back(detid);
    if  (detector_ == "BH") {
      // cells_.dim1.emplace_back(position.eta());
      cells_.dim1[index] = position.eta();
      // cells_.dim2.emplace_back(position.phi());
      cells_.dim2[index] = position.phi();
    }  // else, isSilicon == true and eta phi values will not be used
    else {
      // cells_.dim1.emplace_back(position.x());
      cells_.dim1[index] = position.x();
      // cells_.dim2.emplace_back(position.y());
      cells_.dim2[index] = position.y();
    }
    cells_.weight[index] = hgrh.energy();
    cells_.sigmaNoise[index] = sigmaNoise;
    cells_.layer[index] = layer;
    index ++;
  }

  cells_.resize(index);
}


  /**
   * @brief Sets algoId accordingly to the detector type
  */
  void setAlgoId();

  /**
   * @brief Counts position for all points in the cluster
   * 
   * @param[in] hitmap hitmap to find correct RecHit
   * @param[in] hitsAndFraction all hits in the cluster
   * @return counted position
  */
  math::XYZPoint calculatePosition(std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
                                   const std::vector<std::pair<DetId, float>>& hitsAndFractions);

  /**
   * @brief Counts time for all points in the cluster
   * 
   * @param[in] hitmap hitmap to find correct RecHit only for silicon (not for BH-HSci)
   * @param[in] hitsAndFraction all hits in the cluster
   * @return counted time
  */
  std::pair<float, float> calculateTime(std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
                                        const std::vector<std::pair<DetId, float>>& hitsAndFractions,
                                        size_t sizeCluster);
};

DEFINE_FWK_MODULE(HGCalLayerClusterHeterogenousProducer);

HGCalLayerClusterHeterogenousProducer::HGCalLayerClusterHeterogenousProducer(const edm::ParameterSet& ps)
    : algoId_(reco::CaloCluster::undefined),
      detector_(ps.getParameter<std::string>("detector")),  // one of EE, FH, BH, HFNose
      timeClname_(ps.getParameter<std::string>("timeClname")),
      hitsTime_(ps.getParameter<unsigned int>("nHitsTime")),
      caloGeomToken_(consumesCollector().esConsumes<CaloGeometry, CaloGeometryRecord>())
      {
  initialized_ = false;
  setAlgoId();  //sets algo id according to detector type
  hits_token_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("recHits"));
  isNose_ = false;

  auto pluginPSet = ps.getParameter<edm::ParameterSet>("plugin");

  thresholdW0_ = pluginPSet.getParameter<std::vector<double>>("thresholdW0");
  positionDeltaRho2_ = pluginPSet.getParameter<double>("positionDeltaRho2");
  dependSensor_ = pluginPSet.getParameter<bool>("dependSensor");
  maxNumberOfThickIndices_ = pluginPSet.getParameter<unsigned>("maxNumberOfThickIndices");
  sciThicknessCorrection_ = pluginPSet.getParameter<double>("sciThicknessCorrection"),
  fcPerEle_ = pluginPSet.getParameter<double>("fcPerEle");
  fcPerMip_ = pluginPSet.getParameter<std::vector<double>>("fcPerMip");
  nonAgedNoises_ = pluginPSet.getParameter<std::vector<double>>("noises");
  noiseMip_ = pluginPSet.getParameter<edm::ParameterSet>("noiseMip").getParameter<double>("noise_MIP");
  ecut_ = pluginPSet.getParameter<double>("ecut");
  dEdXweights_ = pluginPSet.getParameter<std::vector<double>>("dEdXweights");
  thicknessCorrection_ = pluginPSet.getParameter<std::vector<double>>("thicknessCorrection");
  deltasi_index_regemfac_ = pluginPSet.getParameter<int>("deltasi_index_regemfac");
  // dc_ = pluginPSet.getParameter<std::vector<double>>("deltac");
        
        

  produces<std::vector<float>>("InitialLayerClustersMask");
  produces<std::vector<reco::BasicCluster>>();
  //time for layer clusters
  produces<edm::ValueMap<std::pair<float, float>>>(timeClname_);
//   std::cout << "Hterogenous" << detector_ << " :" << " end constructor" << std::endl;
}

//todo
void HGCalLayerClusterHeterogenousProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalLayerClusters
  edm::ParameterSetDescription desc;
  edm::ParameterSetDescription pluginDesc;
  pluginDesc.addNode(edm::PluginDescription<HGCalLayerClusterAlgoFactory>("type", "SiCLUE", true));

  desc.add<edm::ParameterSetDescription>("plugin", pluginDesc);
  desc.add<std::string>("detector", "EE")->setComment("options EE, FH, BH,  HFNose; other value defaults to EE");
  desc.add<edm::InputTag>("recHits", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  desc.add<std::string>("timeClname", "timeLayerCluster");
  desc.add<unsigned int>("nHitsTime", 3);
  descriptions.add("hgcalHeterogenousLayerClusters", desc);
}

math::XYZPoint HGCalLayerClusterHeterogenousProducer::calculatePosition(
    std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
    const std::vector<std::pair<DetId, float>>& hitsAndFractions) {
  float total_weight = 0.f;
  float maxEnergyValue = 0.f;
  DetId maxEnergyIndex;
  float x = 0.f;
  float y = 0.f;
  // std::cout << "Hterogenous" << detector_ << " :" << " start calculate position" << std::endl;
  for (auto const& hit : hitsAndFractions) {
    //time is computed wrt  0-25ns + offset and set to -1 if no time
    const HGCRecHit* rechit = hitmap[hit.first];
    total_weight += rechit->energy();
    if (rechit->energy() > maxEnergyValue) {
      maxEnergyValue = rechit->energy();
      maxEnergyIndex = rechit->detid();
    }
  }
  float total_weight_log = 0.f;
  auto thick = rhtools_.getSiThickIndex(maxEnergyIndex);
  const GlobalPoint positionMaxEnergy(rhtools_.getPosition(maxEnergyIndex));
  for (auto const& hit : hitsAndFractions) {
    //time is computed wrt  0-25ns + offset and set to -1 if no time
    const HGCRecHit* rechit = hitmap[hit.first];

    const GlobalPoint position(rhtools_.getPosition(rechit->detid()));

    if (thick != -1) {  //silicon
      //for silicon only just use 1+6 cells = 1.3cm for all thicknesses
      const float d1 = position.x() - positionMaxEnergy.x();
      const float d2 = position.y() - positionMaxEnergy.y();
      if ((d1 * d1 + d2 * d2) > positionDeltaRho2_)
        continue;

      float Wi = std::max(thresholdW0_[thick] + std::log(rechit->energy() / total_weight), 0.);
      x += position.x() * Wi;
      y += position.y() * Wi;
      total_weight_log += Wi;
    } else {  //scintillator
      x += position.x() * rechit->energy();
      y += position.y() * rechit->energy();
    }
  }
  if (thick != -1) {
    total_weight = total_weight_log;
  }
  if (total_weight != 0.) {
    float inv_tot_weight = 1.f / total_weight;
    return math::XYZPoint(x * inv_tot_weight, y * inv_tot_weight, positionMaxEnergy.z());
  } else {
    return math::XYZPoint(0.f, 0.f, 0.f);
  }
}

std::pair<float, float> HGCalLayerClusterHeterogenousProducer::calculateTime(
    std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
    const std::vector<std::pair<DetId, float>>& hitsAndFractions,
    size_t sizeCluster) {
  std::pair<float, float> timeCl(-99., -1.);

  if (sizeCluster >= hitsTime_) {
    std::vector<float> timeClhits;
    std::vector<float> timeErrorClhits;

    for (auto const& hit : hitsAndFractions) {
      //time is computed wrt  0-25ns + offset and set to -1 if no time
      const HGCRecHit* rechit = hitmap[hit.first];

      float rhTimeE = rechit->timeError();
      //check on timeError to exclude scintillator
      if (rhTimeE < 0.)
        continue;
      timeClhits.push_back(rechit->time());
      timeErrorClhits.push_back(1. / (rhTimeE * rhTimeE));
    }
    hgcalsimclustertime::ComputeClusterTime timeEstimator;
    timeCl = timeEstimator.fixSizeHighestDensity(timeClhits, timeErrorClhits, hitsTime_);
  }
  return timeCl;
}
void HGCalLayerClusterHeterogenousProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
//   std::cout << "Hterogenous" << detector_ << " :" << " start produce" << std::endl;
  if(detector_ != "EE" && detector_ != "FH"){ //also HFNose?
    return;
  }
  edm::Handle<HGCRecHitCollection> hits;

  std::unique_ptr<std::vector<reco::BasicCluster>> clusters(new std::vector<reco::BasicCluster>);

  edm::ESHandle<CaloGeometry> geom = es.getHandle(caloGeomToken_);
  rhtools_.setGeometry(*geom);
//   std::cout << "Hterogenous" << detector_ << " :" << " setup geometry, maxLayetr_: " << rhtools_.lastLayer(isNose_) << "isNose: " << isNose_  << std::endl;
  maxlayer_ = rhtools_.lastLayer(isNose_);

  cells_.clear();
  // std::cout << "BEFORE resize" << std::endl;


  //make a map detid-rechit
  // NB for the moment just host EE and FH hits
  // timing in digi for BH not implemented for now
  std::unordered_map<uint32_t, const HGCRecHit*> hitmap;
  evt.getByToken(hits_token_, hits);
  //setup cells_
  cells_.resize((*hits).size()); //is this necessary?
  // std::cout << "AFTER resize" << std::endl;
  // std::cout << "BEFORE populate" << std::endl;
  populate(*hits); 

  cells_.shrink_to_fit();
  // std::cout << "size aftrer fit" << cells_.delta.size() << std::endl;
  const int layers = 96; //todo set to maxLayer_
  // const int layers = 47; //todo set to maxLayer_
  // float dc, float rhoc, float outlierDeltaFactor, bool verbose
  // std::cout << "BEFORE creating clue algo" << std::endl;
  CLUEAlgo<TilesConstants, layers> algoStandelone = CLUEAlgo<TilesConstants, layers>(1.3f,9.f,2.f,false); //todo correct setup
  //  std::cout << "BEFORE set points" << std::endl;
  algoStandelone.setPoints(cells_.dim1, cells_.dim2, cells_.layer, cells_.weight, cells_.sigmaNoise);

  hitmap.reserve((*hits).size());
  for (auto const& it : *hits) {
    hitmap[it.detid().rawId()] = &(it);
  }
  // algoStandelone.setPoints(soa.dim1.size(), &soa.dim1[0], &soa.dim2[0], &soa.layer[0], &soa.weight[0]);
  algoStandelone.makeClusters();
  // std::cout << "BEFORE get points" << std::endl;
  algoStandelone.getPoints(cells_.delta, cells_.nearestHigher, cells_.clusterIndex, cells_.rho, cells_.isSeed, cells_.layer, cells_.weight);
  // std::cout << "BEFORE shrink to fit points" << std::endl;
  // cells_.shrink_to_fit();
  // std::cout << "size aftrer fit" << cells_.delta.size() << std::endl;
  // std::cout << "BEFORE create clusters" << std::endl;
  *clusters = createClusters(algoStandelone.getNumberOfClusters());

//  std::cout << "AFTER create clusters" << std::endl;

  std::vector<std::pair<float, float>> times;
  times.reserve(clusters->size());
// std::cout << "BEFORE calculate position" << std::endl;
  for (unsigned i = 0; i < clusters->size(); ++i) {
    const reco::CaloCluster& sCl = (*clusters)[i];
    (*clusters)[i].setPosition(std::move(calculatePosition(hitmap, sCl.hitsAndFractions()))); 
    if (detector_ != "BH") {
      times.push_back(std::move(calculateTime(hitmap, sCl.hitsAndFractions(), sCl.size())));
    } else {
      times.push_back(std::pair<float, float>(-99., -1.));
    }
  }
// std::cout << "AFTER calculate position" << std::endl;

//   for (auto &i: *clusters){
//       std::cout << i.seed().det() << "," << roundTo5(i.x()) << "," << roundTo5(i.y())<< "," << roundTo5(i.z()) << "," << roundTo5(i.eta()) << "," << roundTo5(i.phi()) << std::endl;
//   }

  auto clusterHandle = evt.put(std::move(clusters));

// std::cout << "447" << std::endl;
  if (detector_ == "HFNose") {
    std::unique_ptr<std::vector<float>> layerClustersMask(new std::vector<float>);
    layerClustersMask->resize(clusterHandle->size(), 1.0);
    evt.put(std::move(layerClustersMask), "InitialLayerClustersMask");
  }
  // std::cout << "480" << std::endl;
  auto timeCl = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
  edm::ValueMap<std::pair<float, float>>::Filler filler(*timeCl);
  filler.insert(clusterHandle, times.begin(), times.end());
  filler.fill();
  evt.put(std::move(timeCl), timeClname_);
  // std::cout << "486" << std::endl;

}

void HGCalLayerClusterHeterogenousProducer::setAlgoId() {
  std::cout << "Heterogenous: " << detector_ << std::endl;
  if (detector_ == "HFNose") {
    algoId_ = reco::CaloCluster::hfnose;
  } else if (detector_ == "EE") {
    algoId_ = reco::CaloCluster::hgcal_em;
  } else {  //for FH or BH
    algoId_ = reco::CaloCluster::hgcal_had;
  }
}
#endif  //__RecoLocalCalo_HGCRecProducers_HGCalLayerClusterHeterogenousProducer_H__
