#ifndef RecoVertex_PixelVertexFinding_plugins_alpaka_vertexFinder_h
#define RecoVertex_PixelVertexFinding_plugins_alpaka_vertexFinder_h

#include <cstddef>
#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackSoA/interface/alpaka/TrackUtilities.h"
#include "DataFormats/VertexSoA/interface/ZVertexDevice.h"
#include "DataFormats/VertexSoA/interface/ZVertexHost.h"
#include "DataFormats/VertexSoA/interface/ZVertexSoA.h"
#include "DataFormats/VertexSoA/interface/alpaka/ZVertexSoACollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoVertex/PixelVertexFinding/interface/PixelVertexWorkSpaceLayout.h"
#include "RecoVertex/PixelVertexFinding/plugins/alpaka/PixelVertexWorkSpaceSoADeviceAlpaka.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexFinder {

  using namespace cms::alpakatools;
  using VtxSoAView = ::reco::ZVertexSoAView;
  using TrkSoAView = ::reco::ZVertexTracksSoAView;
  using WsSoAView = ::vertexFinder::PixelVertexWorkSpaceSoAView;

  class Init {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const &acc, VtxSoAView data, WsSoAView ws) const {
      data.nvFinal() = 0;  // initialization
      ::vertexFinder::init(ws);
    }
  };

  template <typename TrackerTraits>
  class Producer {
    using TkSoAConstView = ::reco::TrackSoAConstView;

  public:
    Producer(bool oneKernel,
             bool useDensity,
             bool useDensityClue,
             bool useDBSCAN,
             bool useIterative,
             bool doSplitting,
             int iminT,      // min number of neighbours to be "core"
             float ieps,     // max absolute distance to cluster
             float ierrmax,  // max error to be "seed"
             float ichi2max,  // max normalized distance to cluster
             float maxChi2ForFirstFit, // Reject outlier tracks that contribute more than this to the chi2 of the initial vertex fit
             float maxChi2ForFinalFit, // Reject outlier tracks that contribute more than this to the chi2 of the final vertex fit
             float maxChi2ForSplit     // Split vertices with a chi2/NDoF greater than this threshold
             )
        : oneKernel_(oneKernel && !(useDBSCAN || useIterative || useDensityClue)),
          useDensity_(useDensity),
          useDensityClue_(useDensityClue),
          useDBSCAN_(useDBSCAN),
          useIterative_(useIterative),
          doSplitting_(doSplitting),
          minT_(iminT),
          eps_(ieps),
          errmax_(ierrmax),
          chi2max_(ichi2max),
          maxChi2ForFirstFit_(maxChi2ForFirstFit),
          maxChi2ForFinalFit_(maxChi2ForFinalFit),
          maxChi2ForSplit_(maxChi2ForSplit) {}

    ~Producer() = default;

    ZVertexSoACollection makeAsync(
        Queue &queue, TkSoAConstView const &tracks_view, int maxVertices, float ptMin, float ptMax) const;

  private:
    const bool oneKernel_;      // run everything (cluster,fit,split,sort) in one kernel. Uses only density clusterizer
    const bool useDensity_;     // use density clusterizer
    const bool useDensityClue_; // use density clusterizer based on CLUE
    const bool useDBSCAN_;      // use DBScan clusterizer
    const bool useIterative_;   // use iterative clusterizer
    const bool doSplitting_;    //run vertex splitting

    const int minT_;       // min number of neighbours to be "core"
    const float eps_;      // max absolute distance to cluster
    const float errmax_;   // max error to be "seed"
    const float chi2max_;  // max normalized distance to cluster

    // Vertex splitting and fitting parameters
    const float maxChi2ForFirstFit_;
    const float maxChi2ForFinalFit_;
    const float maxChi2ForSplit_;

  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexFinder

#endif  // RecoVertex_PixelVertexFinding_plugins_alpaka_vertexFinder_h
