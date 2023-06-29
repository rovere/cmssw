// #include "FWCore/Framework/interface/Event.h"
// #include "FWCore/Framework/interface/ESHandle.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/HGCalReco/interface/HGCalSoACellsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoACellsDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAOutHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAOutDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAOut.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalLayerClustersAlpakaProducer : public stream::EDProducer<> {
    public:
      HGCalLayerClustersAlpakaProducer(edm::ParameterSet const& config)
        : 
          getTokenOutSoA_(consumes(config.getParameter<edm::InputTag>("hgcalOutSoA"))),
          getTokenInSoA_(consumes(config.getParameter<edm::InputTag>("hgcalRecHitsSoA"))),
          hostToken_{produces()},
          timeClToken_{produces()},
          layerClustersMaskToken_{produces()},
          hitsTime_(config.getParameter<unsigned int>("nHitsTime")),
          thresholdW0_(config.getParameter<std::vector<double>>("thresholdW0"))
          {
            caloGeomToken_ = consumesCollector().esConsumes<CaloGeometry, CaloGeometryRecord>();
            positionDeltaRho2_ = config.getParameter<double>("positionDeltaRho2");
            detector_ = config.getParameter<std::string>("detecotor");
            hits_token_ = consumes<HGCRecHitCollection>(config.getParameter<edm::InputTag>("recHits"));
            if (detector_ == "HFNose") {
                algoId_ = reco::CaloCluster::hfnose;
            } else if (detector_ == "EE") {
                algoId_ = reco::CaloCluster::hgcal_em;
            } else {  //for FH or BH
                algoId_ = reco::CaloCluster::hgcal_had;
            }
        }

      ~HGCalLayerClustersAlpakaProducer() override = default;

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {

      if constexpr (std::is_same_v<ALPAKA_ACCELERATOR_NAMESPACE::Device, alpaka_common::DevHost>) {

        edm::ESHandle<CaloGeometry> geom = iSetup.getHandle(caloGeomToken_);
        rhtools_.setGeometry(*geom);

        auto const& cells = iEvent.get(getTokenOutSoA_);
        auto const& inSoA = iEvent.get(getTokenInSoA_);

        hits_ = iEvent.getHandle(hits_token_);
        
        std::unique_ptr<std::vector<reco::BasicCluster>> clusters(new std::vector<reco::BasicCluster>);
        //todo get number of clusters and change 100 to that number
        std::vector<std::pair<float, float>> times;
        *clusters = createClusters(100, inSoA.view().cellsCout(), cells, inSoA, times);


        if (detector_ == "HFNose") {
          std::unique_ptr<std::vector<float>> layerClustersMask(new std::vector<float>);
          layerClustersMask->resize(clusters->size(), 1.0);
          iEvent.emplace(layerClustersMaskToken_, std::move(layerClustersMask));// todo add lable "InitialLayerClustersMask");
        }
        auto timeCl = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
        edm::ValueMap<std::pair<float, float>>::Filler filler(*timeCl);
        //todo not working with cluster it needs cluster handle (because of id)
        // filler.insert(clusters, times.begin(), times.end());
        // filler.fill();

        iEvent.emplace(timeClToken_, std::move(timeCl)); // todo add lable timeClname

        iEvent.emplace(hostToken_, std::move(clusters));
        // auto clusterHandle = iEvent.put(hostToken_, std::move(clusters));
        // auto clusterHandle = iEvent.put(std::move(clusters)); 

      }
    }

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
        edm::ParameterSetDescription desc;
        desc.add<edm::InputTag>("hgcalOutSoA", edm::InputTag("TO BE DEFINED"));
        desc.add<std::string>("detector", "EE");
        descriptions.addWithDefaultLabel(desc);
        desc.add<unsigned int>("nHitsTime", 3);
        desc.add<double>("positionDeltaRho2", 1.69);
        desc.add<std::vector<double>>("thresholdW0", {2.9, 2.9, 2.9});
      }

    private:
      // use device::EDGetToken<T> to read from device memory space
      device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::PortableCollection<HGCalCellsOutSoA>> const getTokenOutSoA_;
      device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::PortableCollection<HGCalCellsSoA>> const getTokenInSoA_;
      edm::EDPutTokenT<std::unique_ptr<std::vector<reco::BasicCluster>>> const hostToken_;
      edm::EDPutTokenT<std::unique_ptr<edm::ValueMap<std::pair<float, float>>>> const timeClToken_;
      edm::EDPutTokenT<std::unique_ptr<std::vector<float>>> const layerClustersMaskToken_;
      reco::CaloCluster::AlgoId algoId_;
      std::string detector_;
      edm::Handle<HGCRecHitCollection> hits_;
      edm::EDGetTokenT<HGCRecHitCollection> hits_token_;
      hgcal::RecHitTools rhtools_;
      edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
      unsigned int hitsTime_;

       // for calculate position
      std::vector<double> thresholdW0_;
      double positionDeltaRho2_;

      math::XYZPoint calculatePosition(std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
       const std::vector<std::pair<DetId, float>>& hitsAndFractions) {
        float total_weight = 0.f;
        float maxEnergyValue = 0.f;
        DetId maxEnergyIndex;
        float x = 0.f;
        float y = 0.f;

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
        }
        return math::XYZPoint(0.f, 0.f, 0.f);
      }
      std::vector<reco::BasicCluster> createClusters(int numberOfClusters, int numberOfCells, const HGCalSoAOutDeviceCollection &cells, const HGCalSoACellsDeviceCollection &inSoA, std::vector<std::pair<float, float>> &times){
        std::unordered_map<uint32_t, const HGCRecHit*> hitmap;
        for (auto const& it : *hits_) {
            hitmap[it.detid().rawId()] = &(it);
        }

        auto cellsView = cells.view();
        auto inSoAView = inSoA.view();

        std::vector<reco::BasicCluster> clusters;
        clusters.resize(numberOfClusters);
        std::set<int> usedIdx;
        for (int i = 0; i < numberOfCells; i++){
            auto cells = cellsView[i];
            auto inSoA = inSoAView[i];
            if (cells.clusterIndex() == -1)
                continue;
            int globalClusterIdx = cells.clusterIndex();
            if (usedIdx.find(globalClusterIdx) != usedIdx.end()){
                clusters[globalClusterIdx].setEnergy(clusters[globalClusterIdx].energy() + inSoA.weight());
                clusters[globalClusterIdx].addHitAndFraction(inSoA.detid(), 1.f);
            }
            else{
                std::vector<std::pair<DetId, float>> thisCluster;
                thisCluster.emplace_back(inSoA.detid(), 1.f);
                math::XYZPoint position = math::XYZPoint(0.f, 0.f, 0.f);
                clusters[globalClusterIdx] = reco::BasicCluster(inSoA.weight(), position, reco::CaloID::DET_HGCAL_ENDCAP, std::move(thisCluster), algoId_ );
                usedIdx.emplace(globalClusterIdx);
            }
            if (cells.isSeed())
                clusters[globalClusterIdx].setSeed(inSoA.detid());
        }

        times.reserve(clusters.size());
        for (unsigned i = 0; i < clusters.size(); ++i) {
            const reco::CaloCluster& sCl = clusters[i];
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
            // hgcalsimclustertime::ComputeClusterTime timeEstimator;  //todo this can not compile
            // timeCl = timeEstimator.fixSizeHighestDensity(timeClhits, timeErrorClhits, hitsTime_);
        }
        return timeCl;
    }
        
  
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalLayerClustersAlpakaProducer);