#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/DumpClustersDetails.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/HGCalReco/interface/HGCalSoACellsHostCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAOutHostCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAOut.h"


class HGCalLayerClustersFromAlpakaProducer : public edm::stream::EDProducer<> {
  public:
    HGCalLayerClustersFromAlpakaProducer(edm::ParameterSet const& config)
      :
        getTokenClustersSoA_(consumes(config.getParameter<edm::InputTag>("hgcalOutSoA"))),
        getTokenCellsSoA_(consumes(config.getParameter<edm::InputTag>("hgcalRecHitsSoA"))),
        hitsTime_(config.getParameter<unsigned int>("nHitsTime")),
        thresholdW0_(config.getParameter<std::vector<double>>("thresholdW0"))
  {
    timeClname_ = config.getParameter<std::string>("timeClname");
    caloGeomToken_ = consumesCollector().esConsumes<CaloGeometry, CaloGeometryRecord>();
    positionDeltaRho2_ = config.getParameter<double>("positionDeltaRho2");
    detector_ = config.getParameter<std::string>("detector");
    hits_token_ = consumes<HGCRecHitCollection>(config.getParameter<edm::InputTag>("recHits"));
    if (detector_ == "HFNose") {
      algoId_ = reco::CaloCluster::hfnose;
    } else if (detector_ == "EE") {
      algoId_ = reco::CaloCluster::hgcal_em;
    } else {  //for FH or BH
      algoId_ = reco::CaloCluster::hgcal_had;
    }

    produces<std::vector<float>>("InitialLayerClustersMask");
    produces<std::vector<reco::BasicCluster>>();
    //time for layer clusters
    produces<edm::ValueMap<std::pair<float, float>>>(timeClname_);
  }

    ~HGCalLayerClustersFromAlpakaProducer() override = default;

    void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override {
      edm::ESHandle<CaloGeometry> geom = iSetup.getHandle(caloGeomToken_);
      rhtools_.setGeometry(*geom);

      auto const& clustersSoA = iEvent.get(getTokenClustersSoA_);
      auto const& cells = iEvent.get(getTokenCellsSoA_);

      hits_ = iEvent.getHandle(hits_token_);

      std::unique_ptr<std::vector<reco::BasicCluster>> clusters(new std::vector<reco::BasicCluster>);

      std::vector<std::pair<float, float>> times;
      unsigned int numberOfClusters = clustersSoA.view().numberOfClustersScalar();
      *clusters = createClusters(numberOfClusters, clustersSoA, cells, times);

      if (detector_ == "HFNose") {
        std::unique_ptr<std::vector<float>> layerClustersMask(new std::vector<float>);
        layerClustersMask->resize(clusters->size(), 1.0);
        iEvent.put( std::move(layerClustersMask), "InitialLayerClustersMask");
      }

    
      // if (detector_ == "EE"){
      //   hgcalUtils::DumpClusters dumper;
      //   dumper.dumpInfos(*clusters, true);
      // }

      auto clusterHandle = iEvent.put(std::move(clusters));
      auto timeCl = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
      edm::ValueMap<std::pair<float, float>>::Filler filler(*timeCl);
      filler.insert(clusterHandle, times.begin(), times.end());
      filler.fill();
      iEvent.put(std::move(timeCl), timeClname_);

    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("hgcalOutSoA", edm::InputTag("hgCalLayerClustersSoAProducer"));
      desc.add<edm::InputTag>("hgcalRecHitsSoA", edm::InputTag("hgCalRecHitsSoAProducer"));
      desc.add<std::string>("detector", "EE");
      desc.add<unsigned int>("nHitsTime", 3);
      desc.add<std::string>("timeClname", "timeLayerCluster");
      desc.add<double>("positionDeltaRho2", 1.69);
      desc.add<std::vector<double>>("thresholdW0", {2.9, 2.9, 2.9});
      desc.add<edm::InputTag>("recHits", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    // use device::EDGetToken<T> to read from device memory space
    edm::EDGetTokenT<HGCalSoAOutHostCollection> const getTokenClustersSoA_;
    edm::EDGetTokenT<HGCalSoACellsHostCollection> const getTokenCellsSoA_;

    reco::CaloCluster::AlgoId algoId_;
    std::string detector_;
    edm::Handle<HGCRecHitCollection> hits_;
    edm::EDGetTokenT<HGCRecHitCollection> hits_token_;
    hgcal::RecHitTools rhtools_;
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
    unsigned int hitsTime_;
    std::string timeClname_;

    // for calculate position
    std::vector<double> thresholdW0_;
    double positionDeltaRho2_;

    math::XYZPoint calculatePosition(std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
        const std::vector<std::pair<DetId, float>>& hitsAndFractions) {
      float total_weight = 0.f;
      float maxEnergyValue = 0.f;
      DetId maxEnergyIndex(0);
      float x = 0.f;
      float y = 0.f;

      for (auto const& hit : hitsAndFractions) {
        //time is computed wrt  0-25ns + offset and set to -1 if no time
        const HGCRecHit* rechit = hitmap[hit.first];
        total_weight += rechit->energy();
        if (rechit->energy() > maxEnergyValue) {
          maxEnergyValue = rechit->energy();
          assert(hit.first == rechit->detid().rawId());
          maxEnergyIndex = hit.first;
        }
      }
      float total_weight_log = 0.f;
      assert(maxEnergyIndex);
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
      }
      return math::XYZPoint(0.f, 0.f, 0.f);
    }

    std::vector<reco::BasicCluster> createClusters(int numberOfClusters,
        const HGCalSoAOutHostCollection &clustersSoA,
        const HGCalSoACellsHostCollection &cells,
        std::vector<std::pair<float, float>> &times) {

      /**
      std::unordered_map<uint32_t, const HGCRecHit*> hitmap;
      for (auto const& it : *hits_) {
        hitmap[it.detid()] = &(it);
      }
      */

      auto clustersSoAView = clustersSoA.view();
      auto cellsView = cells.view();

      //std::vector<int> clustersIdx(clustersSoA->metadata().size());
      // Vector that has 1 entry per cluster. The value is equal to the number
      // of rechits/cells composing each cluster.
      std::vector<int> clustersSize(numberOfClusters, 0);
      //std::iota(clustersIdx.begin(), clustersIdx.end(), 0);
      //std::sort(clustersIdx.begin(), clustersIdx.end(),
              //[&clustersSoAView](int a, int b) { return clustersSoAView[a].clusterIndex() < clustersSoAView[b].clusterIndex(); });

      /**
      int currentCluster = -1;
      size_t i = 0;
      while (i < clustersIdx.size()) {
        int current = clustersSoAView[i].clusterIndex();
        int count = 0;
        size_t j = i + 1;
        // Count consecutive identical entries
        while (j < clustersIdx.size() && clustersSoAView[j].clusterIndex() == current) {
          ++j;
        }
        // Update count and move to the next distinct entry
        if (current >= 0) {
          clustersSize[current] = (j - i);
        }
        i = j;
      }
      */

      std::vector<reco::BasicCluster> clusters;
      clusters.resize(numberOfClusters);
      std::vector<float> total_weight;
      total_weight.resize(numberOfClusters, 0.0f);
      std::vector<float> total_weight_log;
      total_weight_log.resize(numberOfClusters, 0.0f);
      std::vector<float> energy_max;
      energy_max.resize(numberOfClusters, 0.0f);
      std::vector<int> energy_max_idx;
      energy_max_idx.resize(numberOfClusters, -1);
      /**
      std::set<int> usedIdx;
      std::unordered_map<int, int> originalClIdx_CondensedClIx;
      int orderedClIdx = 0;
      */
      for (int i = 0; i < cells->metadata().size(); i++) {
        auto clusterSoAV = clustersSoAView[i];
        auto cellV = cellsView[i];
        if (clusterSoAV.clusterIndex() == -1)
          continue;
        auto globalClusterIdx = clusterSoAV.clusterIndex();
        clustersSize[globalClusterIdx]++;
        clusters[globalClusterIdx].setEnergy(clusters[globalClusterIdx].energy() + cellV.weight());
        clusters[globalClusterIdx].addHitAndFraction(cellV.detid(), 1.f);
        total_weight[globalClusterIdx] += cellV.weight();
        if (cellV.weight() > energy_max[globalClusterIdx]) {
          energy_max[globalClusterIdx] = cellV.weight();
          energy_max_idx[globalClusterIdx] = i;
        }
        /**
        if (usedIdx.find(globalClusterIdx) != usedIdx.end()){
          auto const & currentIdx = originalClIdx_CondensedClIx[globalClusterIdx];
          clusters[currentIdx].setEnergy(clusters[currentIdx].energy() + cellV.weight());
          clusters[currentIdx].addHitAndFraction(cellV.detid(), 1.f);
          assert( cellV.detid() !=0 );
        } else {
          originalClIdx_CondensedClIx[globalClusterIdx] = orderedClIdx;
          std::vector<std::pair<DetId, float>> thisCluster;
          thisCluster.emplace_back(cellV.detid(), 1.f);
          math::XYZPoint position = math::XYZPoint(0.f, 0.f, 0.f);
          clusters[orderedClIdx] = reco::BasicCluster(cellV.weight(), position, reco::CaloID::DET_HGCAL_ENDCAP, std::move(thisCluster), algoId_ );
          usedIdx.emplace(globalClusterIdx);
          orderedClIdx++;
        }
        */
        if (clusterSoAV.isSeed()) {
          clusters[globalClusterIdx].setSeed(cellV.detid());
          clusters[globalClusterIdx].setCaloId(reco::CaloID::DET_HGCAL_ENDCAP);
          clusters[globalClusterIdx].setAlgoId(algoId_);
        }
      }
      std::vector<float> timeCls(std::accumulate(clustersSize.begin(), clustersSize.end(), 0));
      std::vector<float> timeClsError(std::accumulate(clustersSize.begin(), clustersSize.end(), 0));
      std::vector<int> incremental_clustersSize(clustersSize.size());
      std::partial_sum(clustersSize.begin(), clustersSize.end(), incremental_clustersSize.begin());
      std::vector<int> current_clustersSize(clustersSize.size(), 0);
      // Compute the cluster position
      std::vector<float> x_pos;
      x_pos.resize(numberOfClusters, 0.0f);
      std::vector<float> y_pos;
      y_pos.resize(numberOfClusters, 0.0f);
      std::vector<float> z_pos;
      z_pos.resize(numberOfClusters, 0.0f);
      for (int i = 0; i < cells->metadata().size(); i++) {
        auto clusterSoAV = clustersSoAView[i];
        auto cellV = cellsView[i];
        if (clusterSoAV.clusterIndex() == -1)
          continue;
        auto globalClusterIdx = clusterSoAV.clusterIndex();
        int offset = 0;
        if (globalClusterIdx > 0)
          offset = incremental_clustersSize[globalClusterIdx-1];
        int idx = offset + current_clustersSize[globalClusterIdx];
        timeCls[idx] = cellV.time();
        timeClsError[idx] = cellV.time_error();
        current_clustersSize[globalClusterIdx]++;
        auto dim1_i = cellV.dim1() - cellsView[energy_max_idx[globalClusterIdx]].dim1();
        auto dim2_i = cellV.dim2() - cellsView[energy_max_idx[globalClusterIdx]].dim2();
        if ((dim1_i * dim1_i + dim2_i * dim2_i) > positionDeltaRho2_)
          continue;
        //TODO(rovere): check the next assumpion of using fixed 0 index for the thresholds
        float Wi = std::max(thresholdW0_[0] + std::log(cellV.weight() / total_weight[globalClusterIdx]), 0.);
        x_pos[globalClusterIdx] += cellV.dim1() * Wi;
        y_pos[globalClusterIdx] += cellV.dim2() * Wi;
        z_pos[globalClusterIdx] = cellV.z();
        total_weight_log[globalClusterIdx] += Wi;
      }

      times.reserve(clusters.size());
      for (unsigned i = 0; i < clusters.size(); ++i) {
        const reco::CaloCluster& sCl = clusters[i];
        assert(sCl.hitsAndFractions().size() > 0 );
        float inv_tot_weight = 1.f / total_weight_log[i];
        clusters[i].setPosition(math::XYZPoint(x_pos[i] * inv_tot_weight,
                                               y_pos[i] * inv_tot_weight,
                                               z_pos[i]));
        /**
        clusters[i].setPosition(std::move(calculatePosition(hitmap, sCl.hitsAndFractions())));
        */
        if (detector_ != "BH") {
          //times.push_back(std::move(calculateTime(hitmap, sCl.hitsAndFractions(), sCl.size())));
          int start = 0;
          if (i > 0)
            start = incremental_clustersSize[i-1];
          times.push_back(std::move(calculateTime(sCl.size(), &timeCls[start], &timeClsError[start] )));
        } else {
          times.push_back(std::pair<float, float>(-99., -1.));
        }
      }

      return clusters;
    }

    std::pair<float, float> calculateTime(
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

    std::pair<float, float> calculateTime(
        size_t sizeCluster,
        const float * cluster_time,
        const float * cluster_time_error) {
      std::pair<float, float> timeCl(-99., -1.);

      if (sizeCluster >= hitsTime_) {
        std::vector<float> timeClhits(sizeCluster);
        std::vector<float> timeErrorClhits(sizeCluster);

        size_t count = 0;
        while(count < sizeCluster) {
          if (*cluster_time_error >= 0.) {
            timeClhits.push_back(*cluster_time);
            timeErrorClhits.push_back(1. / (*cluster_time_error)*(*cluster_time_error));
          }
          cluster_time++;
          cluster_time_error++;
          count++;
        }
        hgcalsimclustertime::ComputeClusterTime timeEstimator;
        timeCl = timeEstimator.fixSizeHighestDensity(timeClhits, timeErrorClhits, hitsTime_);
      }
      return timeCl;
    }


};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCalLayerClustersFromAlpakaProducer);
