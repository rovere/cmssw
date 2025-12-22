#include <alpaka/alpaka.hpp>


#include <TFormula.h>
#include "CommonTools/Utils/interface/FormulaEvaluator.h"

#include "DataFormats/TrackSoA/interface/TracksHost.h"
#include "DataFormats/TrackSoA/interface/alpaka/TracksSoACollection.h"
#include "DataFormats/TrackSoA/interface/TracksDevice.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/RunningAverage.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoUtilities.h"

#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/PixelSeeding/interface/alpaka/CAGeometrySoACollection.h"
#include "RecoTracker/PixelSeeding/interface/CAGeometryHost.h"
#include "CAHitNtupletGenerator.h"

#include "HeterogeneousCore/AlpakaCore/interface/MoveToDeviceCache.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "RecoTracker/PixelSeeding/interface/CAGeometrySoA.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include <iomanip>

#define GPU_DEBUG

namespace reco {
  struct CAGeometryParams {
    //Constructor from ParameterSet
    CAGeometryParams(edm::ParameterSet const& iConfig)
        : caThetaCuts_(iConfig.getParameter<std::vector<double>>("caThetaCuts")),
          caDCACuts_(iConfig.getParameter<std::vector<double>>("caDCACuts")),
          isStacked_(iConfig.getParameter<std::vector<int32_t>>("isStacked")),
          pairGraph_(iConfig.getParameter<std::vector<unsigned int>>("pairGraph")),
          startingPairs_(iConfig.getParameter<std::vector<unsigned int>>("startingPairs")),
          phiCuts_(iConfig.getParameter<std::vector<int>>("phiCuts")),
          ptCuts_(iConfig.getParameter<std::vector<double>>("ptCuts")),
          minInner_(iConfig.getParameter<std::vector<double>>("minInner")),
          maxInner_(iConfig.getParameter<std::vector<double>>("maxInner")),
          minOuter_(iConfig.getParameter<std::vector<double>>("minOuter")),
          maxOuter_(iConfig.getParameter<std::vector<double>>("maxOuter")),
          maxDZ_(iConfig.getParameter<std::vector<double>>("maxDZ")),
          minDZ_(iConfig.getParameter<std::vector<double>>("minDZ")),
          maxDR_(iConfig.getParameter<std::vector<double>>("maxDR")),
          cellZ0Cuts_(iConfig.getParameter<std::vector<double>>("cellZ0Cuts")) {
      startNoBPix1_ = false;
      for (const unsigned int& i : startingPairs_) {
        if (pairGraph_[2 * i] > 0) {
          startNoBPix1_ = true;
          break;
        }
      }
    }

    // Layers params
    const std::vector<double> caThetaCuts_;
    const std::vector<double> caDCACuts_;
    const std::vector<int32_t> isStacked_;
    const std::vector<int> isBarrel_;

    // Cells params
    const std::vector<unsigned int> pairGraph_;
    const std::vector<unsigned int> startingPairs_;
    const std::vector<int> phiCuts_;
    const std::vector<double> ptCuts_;
    const std::vector<double> minInner_;
    const std::vector<double> maxInner_;
    const std::vector<double> minOuter_;
    const std::vector<double> maxOuter_;
    const std::vector<double> maxDZ_;
    const std::vector<double> minDZ_;
    const std::vector<double> maxDR_;
    const std::vector<double> cellZ0Cuts_;

    bool startNoBPix1_;

    mutable edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tokenGeometry_;
    mutable edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tokenTopology_;
  };

}  // namespace reco

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class CAHitNtupletAlpaka
      : public stream::EDProducer<edm::GlobalCache<::reco::CAGeometryParams>,
                                  edm::RunCache<cms::alpakatools::MoveToDeviceCache<Device, ::reco::CAGeometryHost>>> {
    using HitsConstView = ::reco::TrackingRecHitConstView;
    using HitsOnDevice = reco::TrackingRecHitsSoACollection;
    using HitsOnHost = ::reco::TrackingRecHitHost;

    using TkSoAHost = ::reco::TracksHost;
    using TkSoADevice = reco::TracksSoACollection;

    using Algo = CAHitNtupletGenerator<TrackerTraits>;

    using CAGeometryCache = cms::alpakatools::MoveToDeviceCache<Device, ::reco::CAGeometryHost>;
    using Rotation = SOARotation<float>;
    using Frame = SOAFrame<float>;

  public:
    explicit CAHitNtupletAlpaka(const edm::ParameterSet& iConfig, const ::reco::CAGeometryParams* iCache);
    ~CAHitNtupletAlpaka() override = default;

    void produce(device::Event& iEvent, const device::EventSetup& es) override;

    static void globalEndJob(::reco::CAGeometryParams const*) { /* Do nothing */ };
    static void globalEndRun(edm::Run const& iRun,
                             edm::EventSetup const&,
                             RunContext const* iContext) { /* Do nothing */ };

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    static std::shared_ptr<CAGeometryCache> globalBeginRun(edm::Run const& iRun,
                                                           edm::EventSetup const& iSetup,
                                                           GlobalCache const* iCache) {
      assert(iCache->maxDR_.size() == iCache->minInner_.size());
      assert(iCache->maxDR_.size() == iCache->maxInner_.size());
      assert(iCache->maxDR_.size() == iCache->minOuter_.size());
      assert(iCache->maxDR_.size() == iCache->maxOuter_.size());
      assert(iCache->maxDR_.size() == iCache->maxDZ_.size());
      assert(iCache->maxDR_.size() == iCache->minDZ_.size());
      assert(iCache->maxDR_.size() == iCache->phiCuts_.size());
      assert(iCache->maxDR_.size() == iCache->ptCuts_.size());
      assert(iCache->maxDR_.size() == iCache->cellZ0Cuts_.size());

      assert(iCache->caThetaCuts_.size() == iCache->caDCACuts_.size());
      assert(iCache->caThetaCuts_.size() == iCache->isStacked_.size());

      int n_layers = iCache->caThetaCuts_.size();
      int n_pairs = iCache->pairGraph_.size() / 2;
      int n_modules = 0;
      int n_pixel_modules = 0;

#ifdef GPU_DEBUG
      std::cout << "No. Layers to be used = " << n_layers << std::endl;
      std::cout << "No. Pairs to be used = " << n_pairs << std::endl;
#endif

      assert(int(n_pairs) == int(iCache->maxDR_.size()));
      assert(int(*std::max_element(iCache->startingPairs_.begin(), iCache->startingPairs_.end())) < n_pairs);
      assert(int(*std::max_element(iCache->pairGraph_.begin(), iCache->pairGraph_.end())) < n_layers);

      auto const& trackerGeometry = iSetup.getData(iCache->tokenGeometry_);
      auto const& trackerTopology = iSetup.getData(iCache->tokenTopology_);
      auto const& dets = trackerGeometry.detUnits();

#ifdef GPU_DEBUG
      auto subSystem = 0;
      auto subSystemName = GeomDetEnumerators::tkDetEnum[subSystem];
      std::cout
          << "========================================================================================================="
          << std::endl;
#endif

      auto oldLayer = 0u;
      auto layerCount = 0u;

      std::vector<bool> layerIsBarrel(n_layers);
      std::vector<int> layerStarts(n_layers + 1);
      //^ why n_layers + 1? This is a cumulative sum of the number
      // of modules each layer has. And we need the  extra spot
      // at the end to hold the total number of modules.

      std::vector<int> moduleToindexInDets;

      auto isPinPSinOTBarrel = [&](DetId detId) {
        // Select only P-hits from the OT barrel
        return (trackerGeometry.getDetectorType(detId) == TrackerGeometry::ModuleType::Ph2PSP &&
                detId.subdetId() == StripSubdetector::TOB);
      };
      auto isPixel = [&](DetId detId) {
        auto subId = detId.subdetId();
        return (subId == PixelSubdetector::PixelBarrel || subId == PixelSubdetector::PixelEndcap);
      };
      auto isBarrel = [&](DetId detId) {
        auto subId = detId.subdetId();
        auto subDetector = trackerGeometry.geomDetSubDetector(subId);
        return GeomDetEnumerators::isBarrel(subDetector);
      };

      auto isOTBarrel = [&](DetId detId) { return detId.subdetId() == StripSubdetector::TOB; };

      auto isPSP = [&](const TrackerGeometry& tg, DetId detId) {
        return tg.getDetectorType(detId) == TrackerGeometry::ModuleType::Ph2PSP;
      };

      auto isPSS = [&](const TrackerGeometry& tg, DetId detId) {
        return tg.getDetectorType(detId) == TrackerGeometry::ModuleType::Ph2PSS;
      };

      auto is2S = [&](const TrackerGeometry& tg, DetId detId) {
        // Adjust if your CMSSW uses a different name for 2S modules.
        return tg.getDetectorType(detId) == TrackerGeometry::ModuleType::Ph2SS;
      };

      auto isSelectedOTBarrel = [&](const TrackerGeometry& tg, DetId detId) {
        return isOTBarrel(detId) && (isPSP(tg, detId) || isPSS(tg, detId) || is2S(tg, detId));
      };

      // Robust inner/outer classification for tilted sensors:
      // Use sign of normal x position (position from origin).
      // If negative, the plane normal points toward the beamline => "inner-facing".
      auto innerOuterFromOrientation = [&](const TrackerTopology & tTopo, DetId detId, const GlobalPoint& pos, const GlobalVector& nrm) {
        const double dot = pos.x() * nrm.x() + pos.y() * nrm.y() + pos.z() * nrm.z();
        const bool normalOut = (dot > 0.0);
        // Which member of the stack is this?
        const bool isLower = tTopo.isLower(detId);  // <-- this is the key topo query

        // Map to inner(0)/outer(1)
        // If normal points outward, lower sensor is inner and upper is outer.
        // If normal points inward, mapping flips.
        const int innerOuter = normalOut ? (isLower ? 0 : 1) : (isLower ? 1 : 0);
        return innerOuter;
      };

      // Collect all selected OT barrel sensors with (layer, innerOuter).
      struct ModInfo {
        int detUnitIndex;
        uint32_t rawId;
        int layer;        // real TOB layer
        int innerOuter;   // 0 inner-facing, 1 outer-facing
      };

      std::unordered_map<int, int> moduleIndexToOffset_;  // detUnit index -> offset in SoA 
      std::vector<ModInfo> mods;
      mods.reserve(trackerGeometry.detUnits().size());

      // loop over all detector modules and build the CA layers
      int counter = 0;
      for (auto& det : dets) {
        DetId detid = det->geographicalId();
        auto layer = trackerTopology.layer(detid);
        // Logic:
        // - if we are not inside pixels, we need to ignore anything **but** the OT.
        // - for the time being, this is assuming that the CA extension will
        //   only cover the OT barrel part, and will ignore the OT forward.

#ifdef GPU_DEBUG
        auto subId = detid.subdetId();
        if (subSystemName != trackerGeometry.geomDetSubDetector(subId)) {
          subSystemName = trackerGeometry.geomDetSubDetector(subId);
          std::cout << " ===================== Subsystem: " << subSystemName << " on layer " << layer << std::endl;
        }
#endif

        // Modules of the pixel layers
        if (isPixel(detid)) {
          if (layer != oldLayer) {
#ifdef GPU_DEBUG
            std::cout
              << "PixelLayer "
              << "CA="     << std::setw(2) << layerCount
              << "  subL=" << std::setw(2) << layer
              << "  mod="  << std::setw(5) << n_modules
              << "  cnt="  << std::setw(5) << counter
              << "  type=" << std::left << std::setw(6)
              << (isBarrel(detid) ? "barrel" : "endcap")
              << std::right
              << std::endl;
#endif
            layerIsBarrel[layerCount] = isBarrel(detid);
            layerStarts[layerCount++] = n_modules;
            if (layerCount >= layerStarts.size())
              break;
            oldLayer = layer;
          }
          moduleToindexInDets.push_back(counter);
          n_modules++;
          n_pixel_modules++;
        }

        // if we are using the CA extension for Phase-2,
        // we also have to collect the modules from the considered OT layers
        if constexpr (std::is_same_v<pixelTopology::Phase2OT, TrackerTraits>) {
          // Modules of the considered OT layers
          if (isPinPSinOTBarrel(detid)) {
            if (layer != oldLayer) {
#ifdef GPU_DEBUG
              std::cout
                << "OTLayer   "
                << "CA="     << std::setw(2) << layerCount
                << "  subL=" << std::setw(2) << layer
                << "  mod="  << std::setw(5) << n_modules
                << "  type=" << std::left << std::setw(6)
                << (isBarrel(detid) ? "barrel" : "endcap")
                << std::right
                << std::endl;
#endif
              layerIsBarrel[layerCount] = isBarrel(detid);
              layerStarts[layerCount++] = n_modules;
              if (layerCount >= layerStarts.size())
                break;
              oldLayer = layer;
            }
            moduleToindexInDets.push_back(counter);
            n_modules++;
          }
        }
        // if we are using the CA extension for Phase-2,
        // we also have to collect the modules from the considered OT layers
        if constexpr (std::is_same_v<pixelTopology::Phase2OTFull, TrackerTraits>) {
          uint32_t rawId = detid.rawId();

          if (isSelectedOTBarrel(trackerGeometry, rawId)) {
            moduleToindexInDets.push_back(counter);

            const int layer = trackerTopology.getOTLayerNumber(detid);

            const auto& surf = det->surface();
            const GlobalPoint pos = surf.position();
            const GlobalVector nrm = surf.normalVector();

            const int innerOuter = innerOuterFromOrientation(trackerTopology, rawId, pos, nrm);

            mods.push_back(ModInfo{det->index(), rawId, layer, innerOuter});
          }
        }
        counter++;
      } // end of loop over detUnits from trackerGeometry

      // if we are using the CA extension for Phase-2,
      // we also have to collect the modules from the considered OT layers
      if constexpr (std::is_same_v<pixelTopology::Phase2OTFull, TrackerTraits>) {
        // Sort by (real layer, inner/outer) to ensure grouping:
        //   layer L: all inner-facing first, then all outer-facing
        std::sort(mods.begin(), mods.end(), [](ModInfo const& a, ModInfo const& b) {
          if (a.layer != b.layer)
            return a.layer < b.layer;
          if (a.innerOuter != b.innerOuter)
            return a.innerOuter < b.innerOuter;
          return a.rawId < b.rawId;
        });
        int prevL = -1, prevIO = -1;
        for (size_t i = 0; i < mods.size(); ++i) {
          moduleIndexToOffset_[mods[i].detUnitIndex] = static_cast<int>(i) + n_pixel_modules;
          if (mods[i].layer != prevL || mods[i].innerOuter != prevIO) {
#ifdef GPU_DEBUG
            
          std::cout
            << "Group "
            << "off="   << std::setw(5) << i
            << "  lyr=" << std::setw(2) << mods[i].layer
            << "  io="  << std::setw(1) << mods[i].innerOuter   // 0=inner, 1=outer
            << "  CA="  << std::setw(2) << layerCount
            << "  mod0="<< std::setw(5) << n_modules
            << std::endl;
#endif
            layerIsBarrel[layerCount] = isBarrel(mods[i].rawId);
            layerStarts[layerCount++] = n_modules;
            prevL = mods[i].layer;
            prevIO = mods[i].innerOuter;
          }
          n_modules++;
        }
      }

#ifdef GPU_DEBUG
      std::cout << "Full CA LayerStart: " << n_layers << " layers with " << n_modules << " modules in total."
                << std::endl;
#endif
      layerStarts[n_layers] = n_modules;

      reco::CAGeometryHost product{{{n_layers + 1, n_pairs, n_modules}}, cms::alpakatools::host()};

      auto layerSoA = product.view();
      auto cellSoA = product.view<::reco::CAGraphSoA>();
      auto modulesSoA = product.view<::reco::CAModulesSoA>();

      // TODO the ordering of the insertion of modules should follow the modules numbering of the OT SoA and it does not
      if constexpr (std::is_same_v<pixelTopology::Phase2OTFull, TrackerTraits>) {
        for (int i = 0; i < n_pixel_modules; ++i) {
          auto idx = moduleToindexInDets[i];
          auto det = dets[idx];
          auto vv = det->surface().position();
          auto rr = Rotation(det->surface().rotation());
          modulesSoA[i].detFrame() = Frame(vv.x(), vv.y(), vv.z(), rr);
#ifdef GPU_DEBUG
  //        auto const& detUnits = det->components();
  //        for (auto& detUnit : detUnits) {
  //          DetId unitDetId(detUnit->geographicalId());
  //        }
          std::cout
            << "Frame "
            << "idx="   << std::setw(5) << idx
            << "  soa=" << std::setw(5) << i
            << "  det=" << std::setw(10) << det->geographicalId()
            << "  pos=" << std::setw(15) << vv
            << "\nrot=" << std::setw(25) << det->surface().rotation()
            << "\nz-rot=" << std::setw(8) << std::fixed << std::setprecision(2)
            << atan2(det->surface().normalVector().perp(),
                     det->surface().normalVector().z()) * 180. / M_PI
            << std::endl;
#endif
        }
        for (size_t i = 0; i < mods.size(); ++i) {
          auto idx = mods[i].detUnitIndex;
          auto soaIdx = moduleIndexToOffset_[idx];
          auto det = dets[idx];
          auto vv = det->surface().position();
          auto rr = Rotation(det->surface().rotation());
          modulesSoA[soaIdx].detFrame() = Frame(vv.x(), vv.y(), vv.z(), rr);
#ifdef GPU_DEBUG
  //        auto const& detUnits = det->components();
  //        for (auto& detUnit : detUnits) {
  //          DetId unitDetId(detUnit->geographicalId());
  //        }
          std::cout
            << "Frame "
            << "idx="   << std::setw(5) << idx
            << "  soa=" << std::setw(5) << soaIdx
            << "  det=" << std::setw(10) << det->geographicalId()
            << "  pos=" << std::setw(15) << vv
            << "\nrot=" << std::setw(25) << det->surface().rotation()
            << "\nz-rot=" << std::setw(8) << std::fixed << std::setprecision(2)
            << atan2(det->surface().normalVector().perp(),
                     det->surface().normalVector().z()) * 180. / M_PI
            << std::endl;
#endif
        }
      } else {
        for (int i = 0; i < n_modules; ++i) {
          auto idx = moduleToindexInDets[i];
          auto det = dets[idx];
          auto vv = det->surface().position();
          auto rr = Rotation(det->surface().rotation());
          modulesSoA[i].detFrame() = Frame(vv.x(), vv.y(), vv.z(), rr);
#ifdef GPU_DEBUG
  //        auto const& detUnits = det->components();
  //        for (auto& detUnit : detUnits) {
  //          DetId unitDetId(detUnit->geographicalId());
  //        }
          std::cout
            << "Frame "
            << "idx="   << std::setw(5) << idx
            << "  soa=" << std::setw(5) << i
            << "  det=" << std::setw(10) << det->geographicalId()
            << "  pos=" << std::setw(15) << vv
            << "\nrot=" << std::setw(25) << det->surface().rotation()
            << "\nz-rot=" << std::setw(8) << std::fixed << std::setprecision(2)
            << atan2(det->surface().normalVector().perp(),
                     det->surface().normalVector().z()) * 180. / M_PI
            << std::endl;
#endif
        }
      }

      for (int i = 0; i < n_layers; ++i) {
        layerSoA.layerStarts()[i] = layerStarts[i];
        layerSoA.caThetaCut()[i] = iCache->caThetaCuts_[i];
        layerSoA.caDCACut()[i] = iCache->caDCACuts_[i];
        layerSoA.isStacked()[i] = iCache->isStacked_[i];
        layerSoA.isBarrel()[i] = layerIsBarrel[i];
      }

      layerSoA.layerStarts()[n_layers] = layerStarts[n_layers];

      for (int i = 0; i < n_pairs; ++i) {
        cellSoA.graph()[i] = {{uint32_t(iCache->pairGraph_[2 * i]), uint32_t(iCache->pairGraph_[2 * i + 1])}};
        cellSoA.phiCuts()[i] = iCache->phiCuts_[i];
        // convert ptCut in curvature radius in cm
        // 1 GeV track has 1 GeV/c / (e * 3.8T) ~ 87 cm radius in a 3.8T field
        const float minRadius = iCache->ptCuts_[i] * 87.78f;
        // Use minRadius^2/4 in the CA to avoid sqrt
        const float minRadius2T4 = 4.f * minRadius * minRadius;
        cellSoA.ptCuts()[i] = minRadius2T4;
        cellSoA.minInner()[i] = iCache->minInner_[i];
        cellSoA.maxInner()[i] = iCache->maxInner_[i];
        cellSoA.minOuter()[i] = iCache->minOuter_[i];
        cellSoA.maxOuter()[i] = iCache->maxOuter_[i];
        cellSoA.maxDZ()[i] = iCache->maxDZ_[i];
        cellSoA.minDZ()[i] = iCache->minDZ_[i];
        cellSoA.maxDR()[i] = iCache->maxDR_[i];
        cellSoA.cellZ0Cuts()[i] = iCache->cellZ0Cuts_[i];
        cellSoA.startingPair()[i] = false;
      }

      for (const unsigned int& i : iCache->startingPairs_)
        cellSoA.startingPair()[i] = true;

      return std::make_shared<CAGeometryCache>(std::move(product));
    }

    static std::unique_ptr<::reco::CAGeometryParams> initializeGlobalCache(edm::ParameterSet const& iConfig) {
      return std::make_unique<::reco::CAGeometryParams>(iConfig.getParameterSet("geometry"));
    }

  private:
    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> tokenField_;
    const device::EDGetToken<HitsOnDevice> tokenHit_;
    const device::EDPutToken<TkSoADevice> tokenTrack_;

    const ::reco::FormulaEvaluator maxNumberOfDoublets_;
    const ::reco::FormulaEvaluator maxNumberOfTuples_;

    Algo deviceAlgo_;
  };

  template <typename TrackerTraits>
  CAHitNtupletAlpaka<TrackerTraits>::CAHitNtupletAlpaka(const edm::ParameterSet& iConfig,
                                                        const ::reco::CAGeometryParams* iCache)
      : EDProducer(iConfig),
        tokenField_(esConsumes()),
        tokenHit_(consumes(iConfig.getParameter<edm::InputTag>("pixelRecHitSrc"))),
        tokenTrack_(produces()),
        maxNumberOfDoublets_(iConfig.getParameter<std::string>("maxNumberOfDoublets")),
        maxNumberOfTuples_(iConfig.getParameter<std::string>("maxNumberOfTuples")),
        deviceAlgo_(iConfig) {
    iCache->tokenGeometry_ = esConsumes<edm::Transition::BeginRun>();
    iCache->tokenTopology_ = esConsumes<edm::Transition::BeginRun>();
  }

  template <typename TrackerTraits>
  void CAHitNtupletAlpaka<TrackerTraits>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("pixelRecHitSrc", edm::InputTag("siPixelRecHitsPreSplittingAlpaka"));

    Algo::fillPSetDescription(desc);
    descriptions.addWithDefaultLabel(desc);
  }

  template <typename TrackerTraits>
  void CAHitNtupletAlpaka<TrackerTraits>::produce(device::Event& iEvent, const device::EventSetup& es) {
    auto bf = 1. / es.getData(tokenField_).inverseBzAtOriginInGeV();

    auto const& geometry = runCache()->get(iEvent.queue());
    auto const& hits = iEvent.get(tokenHit_);

    /// Don't bother if no hits on BPix1 and no good graph for that
    /// (so no staring pair without BPix1 as first layer).
    /// TODO: this could be extended to a more general check for
    /// no hits on any of the starting layers.

    if (globalCache()->startNoBPix1_ or hits.offsetBPIX2() > 0) {
      std::array<double, 1> nHitsV = {{double(hits.nHits())}};
      std::array<double, 1> emptyV;

      uint32_t const maxTuples = maxNumberOfTuples_.evaluate(nHitsV, emptyV);
      uint32_t const maxDoublets = maxNumberOfDoublets_.evaluate(nHitsV, emptyV);

      iEvent.emplace(tokenTrack_,
                     deviceAlgo_.makeTuplesAsync(hits, geometry, bf, maxDoublets, maxTuples, iEvent.queue()));

    } else {
      edm::LogWarning("CAHitNtupletAlpaka") << "No hit on BPix1 (" << hits.offsetBPIX2()
                                            << ") and all the starting pairs has BPix1 as inner layer.\nIt's useless "
                                            << "to run the CA. Returning with 0 tracks!";
      auto& queue = iEvent.queue();
      reco::TracksSoACollection tracks({{0, 0}}, queue);
      auto ntracks_d = cms::alpakatools::make_device_view(queue, tracks.view().nTracks());
      alpaka::memset(queue, ntracks_d, 0);
      iEvent.emplace(tokenTrack_, std::move(tracks));
    }
  }

  using CAHitNtupletAlpakaPhase1 = CAHitNtupletAlpaka<pixelTopology::Phase1>;
  using CAHitNtupletAlpakaHIonPhase1 = CAHitNtupletAlpaka<pixelTopology::HIonPhase1>;
  using CAHitNtupletAlpakaPhase2 = CAHitNtupletAlpaka<pixelTopology::Phase2>;
  using CAHitNtupletAlpakaPhase2OT = CAHitNtupletAlpaka<pixelTopology::Phase2OT>;
  using CAHitNtupletAlpakaPhase2OTFull = CAHitNtupletAlpaka<pixelTopology::Phase2OTFull>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"

DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletAlpakaPhase1);
DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletAlpakaHIonPhase1);
DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletAlpakaPhase2);
DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletAlpakaPhase2OT);
DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletAlpakaPhase2OTFull);
