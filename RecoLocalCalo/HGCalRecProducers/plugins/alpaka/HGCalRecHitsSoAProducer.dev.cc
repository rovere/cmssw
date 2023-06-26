#include "DataFormats/PortableTestObjects/interface/alpaka/TestDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaTest/interface/AlpakaESTestRecords.h"
#include "HeterogeneousCore/AlpakaTest/interface/alpaka/AlpakaESTestData.h"

#include "DataFormats/HGCalReco/interface/HGCalSoACellsHostCollection.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalRecHitsSoAProducer : public global::EDProducer<> {
  public:
    HGCalRecHitsSoAProducer(edm::ParameterSet const& config)
        : caloGeomToken_(consumesCollector().esConsumes<CaloGeometry, CaloGeometryRecord>()),
          detector_(config.getParameter<std::string>("detector")),
         
          esToken_(esConsumes(config.getParameter<edm::ESInputTag>("eventSetupSource"))),
          deviceToken_{produces()},
          size_{config.getParameter<edm::ParameterSet>("size").getParameter<int32_t>(
          EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE))} {
            hits_token_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("recHits"));
            isNose_ = false;
            if (detector_ == "HFNose") {
              isNose_ = true;
            }
          }

    void produce(edm::StreamID, device::Event& iEvent, device::EventSetup const& iSetup) const override {

      edm::Handle<HGCRecHitCollection> hits;

      edm::ESHandle<CaloGeometry> geom = es.getHandle(caloGeomToken_);
      rhtools_.setGeometry(*geom);
      maxlayer_ = rhtools_.lastLayer(isNose_);  
      cells_.clear();

      evt.getByToken(hits_token_, hits);
      populate(*hits);

      
      // [[maybe_unused]] auto const& esData = iSetup.getData(esToken_);

      // portabletest::TestDeviceCollection deviceProduct{size_, iEvent.queue()};

      // // run the algorithm, potentially asynchronously
      // algo_.fill(iEvent.queue(), deviceProduct);

      // iEvent.emplace(deviceToken_, std::move(deviceProduct));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add("eventSetupSource", edm::ESInputTag{});

      edm::ParameterSetDescription psetSize;
      psetSize.add<int32_t>("alpaka_serial_sync");
      psetSize.add<int32_t>("alpaka_cuda_async");
      desc.add("size", psetSize);

      descriptions.addWithDefaultLabel(desc);
    }

  private:
    hgcal::RecHitTools rhtools_;
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
    bool isNose_;
    std::string detector_;
    HGCalSoACellsHostCollection cells_;
    edm::EDGetTokenT<HGCRecHitCollection> hits_token_;

    bool dependSensor_;
    bool initialized_;

    unsigned maxNumberOfThickIndices_;
    std::vector<std::vector<double>> thresholds_;
    std::vector<std::vector<double>> v_sigmaNoise_;
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
        if  (detector_ == "BH") {
          cells_[index].dim1() = position.eta();
          cells_[index].dim2() = position.phi();
        }  // else, isSilicon == true and eta phi values will not be used
        else {
          cells_[index].dim1() = position.x();
          cells_[index].dim2() = position.y();
        }
        cells_[index].wight() = hgrh.energy();
        cells_[index].sigmaNoise() = sigmaNoise;
        cells_[index].layer() = layer;
        index ++;
      }

      cells_.resize(index);
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
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalRecHitsSoAProducer);