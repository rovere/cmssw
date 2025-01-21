#include "SimDataFormats/TrackingAnalysis/interface/SimDoublets.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

namespace simdoublets {

  // function that gets the global position of a RecHit based on its Cluster position
  GlobalPoint getGlobalHitPosition(SiPixelRecHitRef const& recHit,
                                   const TrackerGeometry* trackerGeometry,
                                   GlobalVector const& beamSpotPosition,
                                   bool const useClusterLocalPosition = false) {
    // get DetUnit of the RecHit
    DetId detIdObject(recHit->geographicalId());
    const GeomDetUnit* geomDetUnit = trackerGeometry->idToDetUnit(detIdObject);

    LocalPoint localPosition;

    // if use local position of the cluster is set
    if (useClusterLocalPosition) {
      // get position from cluster local position and DetUnit global position
      auto cluster = recHit->cluster();
      MeasurementPoint measurementPoint(cluster->x(), cluster->y());
      localPosition = geomDetUnit->topology().localPosition(measurementPoint);
    }
    // else, use the local position of the RecHit
    else {
      localPosition = recHit->localPositionFast();
    }

    return (geomDetUnit->surface().toGlobal(localPosition) - beamSpotPosition);
  }

  // function that determines the number of skipped layers for a given pair of layer IDs
  // layerIds cover the ranges:
  // 0 to 3 (barrel), 4 to 15 (forward), 16 to 27 (backward)
  int getNumSkippedLayers(std::pair<uint8_t, uint8_t> const& layerIds) {
    bool innerInBarrel = (layerIds.first < 4);
    bool outerInBarrel = (layerIds.second < 4);
    bool innerInBackward = (layerIds.first > 15);
    bool outerInBackward = (layerIds.second > 15);
    bool innerInForward = (!innerInBarrel) && (!innerInBackward);
    bool outerInForward = (!outerInBarrel) && (!outerInBackward);

    // Possibility 0: invalid case (outer layer is not the outer one), set to -1
    if (layerIds.first >= layerIds.second) {
      return -1;
    }
    // Possibility 1: both RecHits lie in the same detector part (barrel, forward or backward)
    else if ((innerInBarrel && outerInBarrel) || (innerInForward && outerInForward) ||
             (innerInBackward && outerInBackward)) {
      return (layerIds.second - layerIds.first - 1);
    }
    // Possibility 2: the inner RecHit is in the barrel while the outer is in either forward or backward
    else if (innerInBarrel) {
      return (outerInForward) ? (layerIds.second - 4) : (layerIds.second - 16);
    }
    // Possibility 3: invalid case (one is forward and the other in backward), set to -1
    else {
      return -1;
    }
  }

  // function that, for a pair of two layers, gives a unique pair Id (innerLayerId * 100 + outerLayerId)
  int getLayerPairId(std::pair<uint8_t, uint8_t> const& layerIds) {
    // calculate the unique layer pair Id as (innerLayerId * 100 + outerLayerId)
    return (layerIds.first * 100 + layerIds.second);
  }

}  // end namespace simdoublets

// SimDoublets::Doublet class member function
// ------------------------------------------

// constructor with explicit setting of useClusterLocalPosition_
SimDoublets::Doublet::Doublet(SimDoublets const& simDoublets,
                              size_t const innerIndex,
                              size_t const outerIndex,
                              const TrackerGeometry* trackerGeometry,
                              bool useClusterLocalPosition)
    : trackerGeometry_(trackerGeometry),
      useClusterLocalPosition_(useClusterLocalPosition),
      trackingParticleRef_(simDoublets.trackingParticle()),
      beamSpotPosition_(simDoublets.beamSpotPosition()) {
  // fill recHits and layers
  recHitRefs_ = std::make_pair(simDoublets.recHits(innerIndex), simDoublets.recHits(outerIndex));
  layerIds_ = std::make_pair(simDoublets.layerIds(innerIndex), simDoublets.layerIds(outerIndex));

  // determine number of skipped layers
  numSkippedLayers_ = simdoublets::getNumSkippedLayers(layerIds_);

  // determine Id of the layer pair
  layerPairId_ = simdoublets::getLayerPairId(layerIds_);
}

// constrcutor with automatic setting of useClusterLocalPosition_
// checks if RecHit localPosition is meaningful: if yes, prefer RecHit localPosition over cluster
SimDoublets::Doublet::Doublet(SimDoublets const& simDoublets,
                              size_t const innerIndex,
                              size_t const outerIndex,
                              const TrackerGeometry* trackerGeometry)
    : SimDoublets::Doublet::Doublet(simDoublets, innerIndex, outerIndex, trackerGeometry, true) {
  // make sure, that there are RecHits
  if (simDoublets.numRecHits() == 0) {
    return;
  }

  // check if the local position of the RecHit makes sense or is always (0,0,0)
  // background is that the local position of RecHit is transient, and therefore not saved to ROOT files
  // for the check look at the local position of the first RecHit
  SiPixelRecHitRef recHit = simDoublets.recHits(0);
  if ((recHit->localPositionFast().x() == 0) && (recHit->localPositionFast().y() == 0)) {
    // if local position is 0, use clusters instead
    useClusterLocalPosition_ = true;
  } else {
    // else, prefer to use the RecHit position (as this is what's used in reco)
    useClusterLocalPosition_ = false;
  }
}

GlobalPoint SimDoublets::Doublet::innerGlobalPos() const {
  // get the inner RecHit's global position
  return simdoublets::getGlobalHitPosition(
      recHitRefs_.first, trackerGeometry_, beamSpotPosition_, useClusterLocalPosition_);
}

GlobalPoint SimDoublets::Doublet::outerGlobalPos() const {
  // get the outer RecHit's global position
  return simdoublets::getGlobalHitPosition(
      recHitRefs_.second, trackerGeometry_, beamSpotPosition_, useClusterLocalPosition_);
}

// SimDoublets class member function
// ---------------------------------

// method to sort the RecHits according to the position
void SimDoublets::sortRecHits(const TrackerGeometry* trackerGeometry) {
  // get the vector of squared magnitudes of the global RecHit positions
  std::vector<double> recHitMag2;
  recHitMag2.reserve(layerIdVector_.size());
  for (const auto& recHit : recHitRefVector_) {
    Global3DPoint globalPosition = simdoublets::getGlobalHitPosition(recHit, trackerGeometry, beamSpotPosition_, true);
    recHitMag2.push_back(globalPosition.mag2());
  }

  // find the permutation vector that sort the magnitudes
  std::vector<std::size_t> sortedPerm(recHitMag2.size());
  std::iota(sortedPerm.begin(), sortedPerm.end(), 0);
  std::sort(sortedPerm.begin(), sortedPerm.end(), [&](std::size_t i, std::size_t j) {
    return (recHitMag2[i] < recHitMag2[j]);
  });

  // create the sorted recHitRefVector and the sorted layerIdVector accordingly
  SiPixelRecHitRefVector sorted_recHitRefVector;
  sorted_recHitRefVector.reserve(sortedPerm.size());
  for (size_t i : sortedPerm) {
    sorted_recHitRefVector.push_back(recHitRefVector_[i]);
  }
  std::vector<uint8_t> sorted_layerIdVector(sortedPerm.size());
  std::transform(sortedPerm.begin(), sortedPerm.end(), sorted_layerIdVector.begin(), [&](std::size_t i) {
    return layerIdVector_[i];
  });

  // swap them with the class member
  recHitRefVector_.swap(sorted_recHitRefVector);
  layerIdVector_.swap(sorted_layerIdVector);

  // set sorted bool to true
  recHitsAreSorted_ = true;
}

// method to produce the true doublets on the fly
std::vector<SimDoublets::Doublet> SimDoublets::getSimDoublets(const TrackerGeometry* trackerGeometry) const {
  // create output vector for the doublets
  std::vector<SimDoublets::Doublet> doubletVector;

  // FIXME maybe change to assertion or error?
  // confirm that the RecHits are sorted
  if (!recHitsAreSorted_) {
    return doubletVector;
  }

  // loop over the RecHits/layer Ids
  for (size_t i = 0; i < layerIdVector_.size(); i++) {
    uint8_t innerLayerId = layerIdVector_[i];
    uint8_t outerLayerId{};
    size_t outerLayerStart{layerIdVector_.size()};

    // find the next layer Id + at which hit this layer starts
    for (size_t j = i + 1; j < layerIdVector_.size(); j++) {
      if (innerLayerId != layerIdVector_[j]) {
        outerLayerId = layerIdVector_[j];
        outerLayerStart = j;
        break;
      }
    }

    // build the doublets of the inner hit i with all outer hits in the layer outerLayerId
    for (size_t j = outerLayerStart; j < layerIdVector_.size(); j++) {
      // break if the hit doesn't belong to the outer layer anymore
      if (outerLayerId != layerIdVector_[j]) {
        break;
      }

      doubletVector.push_back(SimDoublets::Doublet(*this, i, j, trackerGeometry));
    }
  }  // end loop over the RecHits/layer Ids

  return doubletVector;
}