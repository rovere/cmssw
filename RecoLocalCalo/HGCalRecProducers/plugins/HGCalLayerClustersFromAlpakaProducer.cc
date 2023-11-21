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

#define DEBUG_CLUSTERS_ALPAKA 0

class HGCalLayerClustersFromAlpakaProducer : public edm::stream::EDProducer<> {
  public:
    HGCalLayerClustersFromAlpakaProducer(edm::ParameterSet const& config)
      :
        getTokenClustersSoA_(consumes(config.getParameter<edm::InputTag>("hgcalOutSoA"))),
        getTokenCellsSoA_(consumes(config.getParameter<edm::InputTag>("hgcalRecHitsSoA"))),
        // hostToken_{produces()},
        // timeClToken_{produces()},
        // layerClustersMaskToken_{produces()},
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
      //std::cout << fmt::format("Creating {} legacy clusters", numberOfClusters) << std::endl;
      *clusters = createClusters(numberOfClusters, clustersSoA, cells, times);

      //std::cout << "Before create mask" << std::endl;
      if (detector_ == "HFNose") {
        std::unique_ptr<std::vector<float>> layerClustersMask(new std::vector<float>);
        layerClustersMask->resize(clusters->size(), 1.0);
        iEvent.put( std::move(layerClustersMask), "InitialLayerClustersMask");
      }

#if DEBUG_CLUSTERS_ALPAKA
      hgcalUtils::DumpCellsSoA dumperCellsSoA;
      dumperCellsSoA.dumpInfos(cells);

      hgcalUtils::DumpClusters dumper;
      dumper.dumpInfos(*clusters, true);

      hgcalUtils::DumpClustersSoA dumperSoA;
      dumperSoA.dumpInfos(clustersSoA);
#endif

      //std::cout << "Before put clusters" << std::endl;
      auto clusterHandle = iEvent.put(std::move(clusters));
      auto timeCl = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
      edm::ValueMap<std::pair<float, float>>::Filler filler(*timeCl);
      filler.insert(clusterHandle, times.begin(), times.end());
      filler.fill();
      //std::cout << "Before put time" << std::endl;
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

    // edm::EDPutTokenT<std::unique_ptr<std::vector<reco::BasicCluster>>> const hostToken_;
    // edm::EDPutTokenT<std::unique_ptr<edm::ValueMap<std::pair<float, float>>>> const timeClToken_;
    // edm::EDPutTokenT<std::unique_ptr<std::vector<float>>> const layerClustersMaskToken_;

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

      std::unordered_map<uint32_t, const HGCRecHit*> hitmap;
      for (auto const& it : *hits_) {
        hitmap[it.detid()] = &(it);
      }

      auto clustersSoAView = clustersSoA.view();
      auto cellsView = cells.view();

      std::vector<reco::BasicCluster> clusters;
      clusters.resize(numberOfClusters);
      std::set<int> usedIdx;
      std::unordered_map<int, int> originalClIdx_CondensedClIx;
      int orderedClIdx = 0;
      for (int i = 0; i < cells->metadata().size(); i++) {
        auto clusterSoAV = clustersSoAView[i];
        auto cellV = cellsView[i];
        if (clusterSoAV.clusterIndex() == -1)
          continue;
        auto globalClusterIdx = clusterSoAV.clusterIndex();
        if (usedIdx.find(globalClusterIdx) != usedIdx.end()){
          //std::cout << fmt::format("Re-Converting Original Cl {} into legacy {}", globalClusterIdx, originalClIdx_CondensedClIx[globalClusterIdx]) << std::endl;
          auto const & currentIdx = originalClIdx_CondensedClIx[globalClusterIdx];
          clusters[currentIdx].setEnergy(clusters[currentIdx].energy() + cellV.weight());
          clusters[currentIdx].addHitAndFraction(cellV.detid(), 1.f);
          assert( cellV.detid() !=0 );
        } else {
          originalClIdx_CondensedClIx[globalClusterIdx] = orderedClIdx;
          //std::cout << fmt::format("Converting Original Cl {} into legacy {}", globalClusterIdx, originalClIdx_CondensedClIx[globalClusterIdx]) << std::endl;
          std::vector<std::pair<DetId, float>> thisCluster;
          thisCluster.emplace_back(cellV.detid(), 1.f);
          math::XYZPoint position = math::XYZPoint(0.f, 0.f, 0.f);
          clusters[orderedClIdx] = reco::BasicCluster(cellV.weight(), position, reco::CaloID::DET_HGCAL_ENDCAP, std::move(thisCluster), algoId_ );
          usedIdx.emplace(globalClusterIdx);
          orderedClIdx++;
        }
        if (clusterSoAV.isSeed())
          clusters[originalClIdx_CondensedClIx[globalClusterIdx]].setSeed(cellV.detid());
      }

      times.reserve(clusters.size());
      for (unsigned i = 0; i < clusters.size(); ++i) {
        const reco::CaloCluster& sCl = clusters[i];
        //std::cout << fmt::format("Cluster {} has {} cells", i, sCl.hitsAndFractions().size()) << std::endl;
        assert(sCl.hitsAndFractions().size() > 0 );
        //if (sCl.hitsAndFractions().size() == 0)
          //continue;
        clusters[i].setPosition(std::move(calculatePosition(hitmap, sCl.hitsAndFractions())));
        if (detector_ != "BH") {
          times.push_back(std::move(calculateTime(hitmap, sCl.hitsAndFractions(), sCl.size())));
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


};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCalLayerClustersFromAlpakaProducer);
