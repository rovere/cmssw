// Phase2OTRecHitsFullSoAConverter.cc
//
// Convert Phase2 OT rechits (PS + 2S) into a TrackingRecHitHost SoA,
// ordering OT barrel sensors as:
//   for each real TOB layer L: [all INNER-facing sensors] then [all OUTER-facing sensors]
//
// INNER/OUTER is decided from the sensor plane orientation:
//   innerOuter = (normal · position < 0) ? INNER : OUTER
// which is robust against tilted sensors.
//
// Notes:
//  - Assumes TrackerGeometry::detUnits() provides sensor-level GeomDetUnits.
//  - You may need to adjust the 2S module enum name (Ph2SS here) for your CMSSW.

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Math/interface/approx_atan2.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

//#define HITS_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class Phase2OTRecHitsFullSoAConverter : public stream::EDProducer<> {
    using Hits = ::reco::TrackingRecHitHost;
    using HMSstorage = std::vector<uint32_t>;

  public:
    explicit Phase2OTRecHitsFullSoAConverter(const edm::ParameterSet& iConfig);
    ~Phase2OTRecHitsFullSoAConverter() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    void beginRun(edm::Run const& run, edm::EventSetup const& setup) override;

  private:
    void produce(device::Event& iEvent, const device::EventSetup& es) override;

    // ES
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomTokenRun_;
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoTokenRun_;

    // ED
    const edm::EDGetTokenT<Phase2TrackerRecHit1DCollectionNew> recHitToken_;
    const edm::EDGetTokenT<::reco::BeamSpot> beamSpotToken_;
    const edm::EDGetTokenT<Hits> pixelHitsSoA_;

    // outputs
    const edm::EDPutTokenT<Hits> stripSoA_;
    const edm::EDPutTokenT<HMSstorage> hitModuleStart_;

    // cached at beginRun
    int modulesInPixel_ = 0;

    // Selected OT barrel sensors (PS + 2S), and their ordering
    std::unordered_map<uint32_t, bool> detIdIsSelectedOTBarrel_;   // rawId -> selected
    std::unordered_map<uint32_t, int> rawIdToDetUnitIndex_;        // rawId -> detUnit->index()
    std::vector<int> orderedModules_;                              // detUnit->index(), sorted by (layer, inner/outer)
    std::unordered_map<int, int> moduleIndexToOffset_;             // detUnit index -> offset in orderedModules_

    // debug / tagging (optional)
    std::unordered_map<uint32_t, uint16_t> detIdToRealLayer_;      // rawId -> TOB layer
    std::unordered_map<uint32_t, uint8_t> detIdToInnerOuter_;      // rawId -> 0 inner, 1 outer

  private:
    static inline bool isPh2Pixel(DetId detId) {
      auto subId = detId.subdetId();
      return (subId == PixelSubdetector::PixelBarrel || subId == PixelSubdetector::PixelEndcap);
    }

    static inline bool isOTBarrel(DetId detId) { return detId.subdetId() == StripSubdetector::TOB; }

    static inline bool isPSP(const TrackerGeometry* tg, DetId detId) {
      return tg->getDetectorType(detId) == TrackerGeometry::ModuleType::Ph2PSP;
    }

    static inline bool isPSS(const TrackerGeometry* tg, DetId detId) {
      return tg->getDetectorType(detId) == TrackerGeometry::ModuleType::Ph2PSS;
    }

    static inline bool is2S(const TrackerGeometry* tg, DetId detId) {
      // Adjust if your CMSSW uses a different name for 2S modules.
      return tg->getDetectorType(detId) == TrackerGeometry::ModuleType::Ph2SS;
    }

    static inline bool isSelectedOTBarrel(const TrackerGeometry* tg, DetId detId) {
      return isOTBarrel(detId) && (isPSP(tg, detId) || isPSS(tg, detId) || is2S(tg, detId));
    }

    // Robust inner/outer classification for tilted sensors:
    // Use sign of normal x position (position from origin).
    // If negative, the plane normal points toward the beamline => "inner-facing".
    static inline int innerOuterFromOrientation(const TrackerTopology & tTopo, DetId detId, const GlobalPoint& pos, const GlobalVector& nrm) {
      const double dot = pos.x() * nrm.x() + pos.y() * nrm.y() + pos.z() * nrm.z();
      const bool normalOut = (dot > 0.0);
      // Which member of the stack is this?
      const bool isLower = tTopo.isLower(detId);  // <-- this is the key topo query

      // Map to inner(0)/outer(1)
      // If normal points outward, lower sensor is inner and upper is outer.
      // If normal points inward, mapping flips.
      const int innerOuter = normalOut ? (isLower ? 0 : 1) : (isLower ? 1 : 0);
      return innerOuter;
    }
  };

  Phase2OTRecHitsFullSoAConverter::Phase2OTRecHitsFullSoAConverter(const edm::ParameterSet& iConfig)
      : stream::EDProducer<>(iConfig),
        geomToken_(esConsumes()),
        geomTokenRun_(esConsumes<edm::Transition::BeginRun>()),
        topoTokenRun_(esConsumes<edm::Transition::BeginRun>()),
        recHitToken_(consumes(iConfig.getParameter<edm::InputTag>("otRecHitSource"))),
        beamSpotToken_(consumes<::reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
        pixelHitsSoA_(consumes(iConfig.getParameter<edm::InputTag>("pixelRecHitSoASource"))),
        stripSoA_(produces()),
        hitModuleStart_(produces()) {}

  void Phase2OTRecHitsFullSoAConverter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("pixelRecHitSoASource", edm::InputTag("hltPhase2SiPixelRecHitsSoA"));
    desc.add<edm::InputTag>("otRecHitSource", edm::InputTag("hltSiPhase2RecHits"));
    desc.add<edm::InputTag>("beamSpot", edm::InputTag("hltOnlineBeamSpot"));
    descriptions.addWithDefaultLabel(desc);
  }

  void Phase2OTRecHitsFullSoAConverter::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {
    modulesInPixel_ = 0;
    detIdIsSelectedOTBarrel_.clear();
    rawIdToDetUnitIndex_.clear();
    orderedModules_.clear();
    moduleIndexToOffset_.clear();
    detIdToRealLayer_.clear();
    detIdToInnerOuter_.clear();

    const auto* trackerGeometry = &iSetup.getData(geomTokenRun_);
    const auto& tTopo = iSetup.getData(topoTokenRun_);

    // Collect all selected OT barrel sensors with (layer, innerOuter).
    struct ModInfo {
      int detUnitIndex;
      uint32_t rawId;
      int layer;        // real TOB layer
      int innerOuter;   // 0 inner-facing, 1 outer-facing
    };

    std::vector<ModInfo> mods;
    mods.reserve(trackerGeometry->detUnits().size());

    for (auto const& detUnit : trackerGeometry->detUnits()) {
      DetId detId(detUnit->geographicalId());
      uint32_t rawId = detId.rawId();

      if (isPh2Pixel(detId)) {
        modulesInPixel_++;
      }

      if (!isSelectedOTBarrel(trackerGeometry, detId)) {
        detIdIsSelectedOTBarrel_[rawId] = false;
        continue;
      }

      detIdIsSelectedOTBarrel_[rawId] = true;
      rawIdToDetUnitIndex_[rawId] = detUnit->index();

      const int layer = tTopo.getOTLayerNumber(detId);

      const auto& surf = detUnit->surface();
      const GlobalPoint pos = surf.position();
      const GlobalVector nrm = surf.normalVector();

      const int innerOuter = innerOuterFromOrientation(tTopo, rawId, pos, nrm);

      detIdToRealLayer_[rawId] = static_cast<uint16_t>(layer);
      detIdToInnerOuter_[rawId] = static_cast<uint8_t>(innerOuter);

      mods.push_back(ModInfo{detUnit->index(), rawId, layer, innerOuter});
    }

    // Sort by (real layer, inner/outer) to ensure grouping:
    //   layer L: all inner-facing first, then all outer-facing
    std::sort(mods.begin(), mods.end(), [](ModInfo const& a, ModInfo const& b) {
      if (a.layer != b.layer)
        return a.layer < b.layer;
      if (a.innerOuter != b.innerOuter)
        return a.innerOuter < b.innerOuter;
      return a.rawId < b.rawId;
    });

    orderedModules_.reserve(mods.size());
    for (size_t i = 0; i < mods.size(); ++i) {
      orderedModules_.push_back(mods[i].detUnitIndex);
      moduleIndexToOffset_[mods[i].detUnitIndex] = static_cast<int>(i);
      LogDebug("Phase2OTRecHitsFullSoAConverter") << "After Sorting " << mods[i].detUnitIndex << " " << orderedModules_.size()
        << " on layer " << mods[i].layer << " innerOuter " << mods[i].innerOuter <<'\n';
    }

    LogDebug("Phase2OTRecHitsFullSoAConverter")
        << "modulesInPixel_=" << modulesInPixel_ << "\n"
        << "selected OT barrel sensor detUnits=" << orderedModules_.size() << "\n";

#ifdef HITS_DEBUG
    int prevL = -1, prevIO = -1;
  for (size_t i = 0; i < mods.size(); ++i) {
    if (mods[i].layer != prevL || mods[i].innerOuter != prevIO) {
      std::cout << "Group starts at offset " << i << " : layer=" << mods[i].layer
        << " innerOuter=" << mods[i].innerOuter << " (0=inner,1=outer)\n";
      prevL = mods[i].layer;
      prevIO = mods[i].innerOuter;
    }
  }
#endif
}

  void Phase2OTRecHitsFullSoAConverter::produce(device::Event& iEvent, device::EventSetup const& iSetup) {
    auto queue = iEvent.queue();

    auto const& bs = iEvent.get(beamSpotToken_);
    const auto* trackerGeometry = &iSetup.getData(geomToken_);
    const auto& otHits = iEvent.get(recHitToken_);
    const auto& pixelHitsSoA = iEvent.get(pixelHitsSoA_);
    int nPixelHits = pixelHitsSoA.view().metadata().size();

    // Count total selected OT hits (PS + 2S in TOB)
    int nSelectedHits = 0;
    for (auto const& detSet : otHits) {
      for (auto const& recHit : detSet) {
        uint32_t rawId = recHit.geographicalId().rawId();
        auto it = detIdIsSelectedOTBarrel_.find(rawId);
        if (it != detIdIsSelectedOTBarrel_.end() && it->second) {
          nSelectedHits++;
        }
      }
    }

    LogDebug("Phase2OTRecHitsFullSoAConverter")
        << "nPixelHits=" << nPixelHits << "\n"
        << "nSelectedOTHits=" << nSelectedHits << "\n"
        << "nSelectedOTModules=" << orderedModules_.size() << "\n";

    // Allocate SoA: (selected OT hits, selected OT modules)
    Hits outHitsSoA(queue, nSelectedHits, orderedModules_.size());
    auto& moduleView = outHitsSoA.view<::reco::HitModuleSoA>();

    // hits per module (in our orderedModules_ order)
    std::vector<int> hitsPerModule(orderedModules_.size(), 0);

    for (auto const& detSet : otHits) {
      if (detSet.empty())
        continue;

      uint32_t rawId = detSet.begin()->geographicalId().rawId();
      auto itSel = detIdIsSelectedOTBarrel_.find(rawId);
      if (itSel == detIdIsSelectedOTBarrel_.end() || !itSel->second)
        continue;

      auto itIdx = rawIdToDetUnitIndex_.find(rawId);
      if (itIdx == rawIdToDetUnitIndex_.end())
        continue;

      int detUnitIndex = itIdx->second;
      auto itOff = moduleIndexToOffset_.find(detUnitIndex);
      if (itOff == moduleIndexToOffset_.end())
        continue;

      hitsPerModule[itOff->second] = static_cast<int>(detSet.size());
    }

    // Build cumulative sums to assign moduleStart()
    std::vector<int> cumulative(hitsPerModule.size(), 0);
    std::partial_sum(hitsPerModule.begin(), hitsPerModule.end(), cumulative.begin());

    // moduleStart is shifted by nPixelHits (OT hits come after pixel hits in the consumer convention)
    if (!orderedModules_.empty()) {
      moduleView[0].moduleStart() = nPixelHits;
      LogDebug("Phase2OTRecHitsFullSoAConverter")
            << "Module start: 0 with hits: " << moduleView[0].moduleStart() << '\n';
      for (size_t i = 1; i < cumulative.size(); ++i) {
        moduleView[i].moduleStart() = cumulative[i - 1] + nPixelHits;
        LogDebug("Phase2OTRecHitsFullSoAConverter")
            << "Module start: " << i << " with hits: " << moduleView[i].moduleStart() << '\n';
      }
      moduleView[orderedModules_.size()].moduleStart() = cumulative.back() + nPixelHits;
    } else {
      moduleView[0].moduleStart() = nPixelHits;
    }

    // Fill hits
    for (auto const& detSet : otHits) {
      if (detSet.empty())
        continue;

      uint32_t rawId = detSet.begin()->geographicalId().rawId();

      auto itSel = detIdIsSelectedOTBarrel_.find(rawId);
      if (itSel == detIdIsSelectedOTBarrel_.end() || !itSel->second)
        continue;

      auto itIdx = rawIdToDetUnitIndex_.find(rawId);
      if (itIdx == rawIdToDetUnitIndex_.end())
        continue;

      int detUnitIndex = itIdx->second;

      auto itOff = moduleIndexToOffset_.find(detUnitIndex);
      if (itOff == moduleIndexToOffset_.end())
        continue;

      int offset = itOff->second;

      // index within OT-only block
      int moduleHitIndex = (offset == 0 ? 0 : cumulative[offset - 1]);

      auto det = trackerGeometry->idToDet(DetId(rawId));
      if (det == nullptr) {
        edm::LogWarning("Phase2OTRecHitsFullSoAConverter") << "Null det for rawId=" << rawId;
        continue;
      }

      for (auto const& recHit : detSet) {
        uint32_t hitRawId = recHit.geographicalId().rawId();
        auto itS = detIdIsSelectedOTBarrel_.find(hitRawId);
        if (itS == detIdIsSelectedOTBarrel_.end() || !itS->second)
          continue;

        int idx = moduleHitIndex++;
        assert(idx < nSelectedHits);

        auto hit = outHitsSoA.view()[idx];

        hit.xLocal() = recHit.localPosition().x();
        hit.yLocal() = recHit.localPosition().y();
        hit.xerrLocal() = recHit.localPositionError().xx();
        hit.yerrLocal() = recHit.localPositionError().yy();

        auto globalPosition = det->toGlobal(recHit.localPosition());
        double gx = globalPosition.x() - bs.x0();
        double gy = globalPosition.y() - bs.y0();
        double gz = globalPosition.z() - bs.z0();

        hit.xGlobal() = gx;
        hit.yGlobal() = gy;
        hit.zGlobal() = gz;
        hit.rGlobal() = std::sqrt(gx * gx + gy * gy);
        hit.iphi() = unsafe_atan2s<7>(gy, gx);

        hit.chargeAndStatus().charge = 0;
        hit.chargeAndStatus().status = {false, false, false, false, 0};
        hit.clusterSizeX() = -1;
        hit.clusterSizeY() = -1;

        // Pixel modules are [0..modulesInPixel_-1], OT modules follow in our (layer,innerOuter) order.
        hit.detectorIndex() = modulesInPixel_ + offset;

        LogDebug("Phase2OTRecHitsFullSoAConverter")
            << "Filled OT hit idx=" << idx << " detIndex=" << hit.detectorIndex() << " rawId=" << hitRawId
            << " layer=" << detIdToRealLayer_[hitRawId] << " innerOuter=" << int(detIdToInnerOuter_[hitRawId])
            << " Local (x, y) with (xx, yy) --> (" << recHit.localPosition().x() << ", "
            << recHit.localPosition().y() << ") with (" << recHit.localPositionError().xx() << ", "
            << recHit.localPositionError().yy() << ")" << '\n'
            << "Global           (x, y, z) --> (" << globalPosition.x() << ", " << globalPosition.y() << ", "
            << globalPosition.z() << ")" << '\n'
            << "Corrected Global (x, y, z) --> (" << gx << ", " << gy << ", " << gz << ")" << '\n'
            << gx
            << "\n";
      }  // end of loop over RecHits in detSet
    }  // end of loop over detSet in OT hit container

#ifdef HITS_DEBUG
    for (int h = 0; h < moduleView.metadata().size(); ++h) {
      std::cout << "module " << h << " start=" << moduleView[h].moduleStart() << "\n";
    }
#endif

    // Produce moduleStart vector (legacy-style consumer convenience)
    HMSstorage moduleStartVec(moduleView.metadata().size());
    std::memcpy(moduleStartVec.data(),
                moduleView.moduleStart().data(),
                sizeof(uint32_t) * moduleView.metadata().size());
    iEvent.emplace(hitModuleStart_, std::move(moduleStartVec));

    // Produce OT hits SoA (host); framework handles host->device if needed downstream
    iEvent.emplace(stripSoA_, std::move(outHitsSoA));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(Phase2OTRecHitsFullSoAConverter);

