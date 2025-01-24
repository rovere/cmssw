// -*- C++ -*-
//
// Package:    Validation/TrackingMCTruth
// Class:      SimDoubletsAnalyzer
//
/**\class SimDoubletsAnalyzer SimDoubletsAnalyzer.cc Validation/TrackingMCTruth/plugins/SimDoubletsAnalyzer.cc

 Description: DQM analyzer for true RecHit doublets (SimDoublets) of the inner tracker in Phase-2 HLT

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Luca Ferragina, Elena Vernazza, Jan Schulz
//         Created:  Thu, 16 Jan 2025 13:46:21 GMT
//
//

#include <string>
#include <map>
#include <vector>

// user include files
#include "DataFormats/Histograms/interface/MonitorElementCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/approx_atan2.h"
#include "SimDataFormats/TrackingAnalysis/interface/SimDoublets.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"

// -------------------------------------------------------------------------------------------------------------
// class declaration
// -------------------------------------------------------------------------------------------------------------

class SimDoubletsAnalyzer : public DQMEDAnalyzer {
public:
  explicit SimDoubletsAnalyzer(const edm::ParameterSet&);
  ~SimDoubletsAnalyzer() override;

  void dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ------------ member data ------------

  const TrackerGeometry* trackerGeometry_ = nullptr;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geometry_getToken_;
  const edm::EDGetTokenT<SimDoubletsCollection> simDoublets_getToken_;

  // cutting parameters
  std::vector<double> cellMinz_;
  std::vector<double> cellMaxz_;
  std::vector<int> cellPhiCuts_;
  std::vector<double> cellMaxr_;
  int cellMinYSizeB1_;
  int cellMinYSizeB2_;
  int cellMaxDYSize12_;
  int cellMaxDYSize_;
  int cellMaxDYPred_;
  double cellZ0Cut_;
  double cellPtCut_;

  std::string folder_;  // main folder in the DQM file
  int eventCount_ = 0;  // event counter

  // monitor elements (histograms) to be filled
  MonitorElement* h_layerPairs_;
  MonitorElement* h_numSkippedLayers_;
  MonitorElement* h_numTotVsPt_;
  MonitorElement* h_numPassVsPt_;
  MonitorElement* h_numTotVsEta_;
  MonitorElement* h_numPassVsEta_;
  MonitorElement* h_z0_;
  MonitorElement* h_curvatureR_;
  MonitorElement* h_pTFromR_;
  MonitorElement* h_YsizeB1_;
  MonitorElement* h_YsizeB2_;
  MonitorElement* h_DYsize12_;
  MonitorElement* h_DYsize_;
  MonitorElement* h_DYPred_;
  std::vector<MonitorElement*> hVector_dr_;
  std::vector<MonitorElement*> hVector_dphi_;
  std::vector<MonitorElement*> hVector_idphi_;
  std::vector<MonitorElement*> hVector_innerZ_;
  std::vector<MonitorElement*> hVector_clusterSizeY_;
  std::vector<MonitorElement*> hVector_dsizeYonlyBarrel_;
  std::vector<MonitorElement*> hVector_dsizeYinnerBarrel_;
};

namespace simdoublets {
  std::pair<std::string, std::string> getInnerOuterLayerNames(int const layerPairId) {
    // make a string from the Id (int)
    std::string index = std::to_string(layerPairId);
    // determine inner and outer layer name
    std::string innerLayerName;
    std::string outerLayerName;
    if (index.size() < 3) {
      innerLayerName = "0";
      outerLayerName = index;
    } else if (index.size() == 3) {
      innerLayerName = index.substr(0, 1);
      outerLayerName = index.substr(1, 3);
    } else {
      innerLayerName = index.substr(0, 2);
      outerLayerName = index.substr(2, 4);
    }
    if (outerLayerName[0] == '0') {
      outerLayerName = outerLayerName.substr(1, 2);
    }

    return {innerLayerName, outerLayerName};
  }

  // make bins logarithmic
  void BinLogX(TH1* h) {
    TAxis* axis = h->GetXaxis();
    int bins = axis->GetNbins();

    float from = axis->GetXmin();
    float to = axis->GetXmax();
    float width = (to - from) / bins;
    std::vector<float> new_bins(bins + 1, 0);

    for (int i = 0; i <= bins; i++) {
      new_bins[i] = TMath::Power(10, from + i * width);
    }
    axis->Set(bins, new_bins.data());
  }

  // function to produce histogram with log scale on x (taken from MultiTrackValidator)
  template <typename... Args>
  dqm::reco::MonitorElement* make1DLogX(dqm::reco::DQMStore::IBooker& ibook, Args&&... args) {
    auto h = std::make_unique<TH1F>(std::forward<Args>(args)...);
    BinLogX(h.get());
    const auto& name = h->GetName();
    return ibook.book1D(name, h.release());
  }
}  // namespace simdoublets

//
// static data member definitions
//
// map that takes the layerPairId as defined in the SimDoublets
// and gives the position of the histogram in the histogram vector
// NOTE: It is absolutely necessary that the map is sorted here,
// otherwise the histograms will not be labeled corresponding to the correct layer pair but are mixed up
static const std::map<int, int> layerPairId2Index{
    {1, 0},     {4, 1},     {16, 2},    {102, 3},   {104, 4},   {116, 5},   {203, 6},   {204, 7},   {216, 8},
    {405, 9},   {506, 10},  {607, 11},  {708, 12},  {809, 13},  {910, 14},  {1011, 15}, {1617, 16}, {1718, 17},
    {1819, 18}, {1920, 19}, {2021, 20}, {2122, 21}, {2223, 22}, {2, 23},    {5, 24},    {17, 25},   {6, 26},
    {18, 27},   {103, 28},  {105, 29},  {117, 30},  {106, 31},  {118, 32},  {1112, 33}, {1213, 34}, {1314, 35},
    {1415, 36}, {2324, 37}, {2425, 38}, {2526, 39}, {2627, 40}, {406, 41},  {507, 42},  {608, 43},  {709, 44},
    {810, 45},  {911, 46},  {1012, 47}, {1618, 48}, {1719, 49}, {1820, 50}, {1921, 51}, {2022, 52}, {2123, 53},
    {2224, 54}, {1315, 55}, {217, 56},  {205, 57},  {2325, 58}, {1113, 59}, {2426, 60}, {1214, 61}, {2527, 62}};

static const size_t numLayerPairs = layerPairId2Index.size();

// -------------------------------
// constructors and destructor
// -------------------------------
SimDoubletsAnalyzer::SimDoubletsAnalyzer(const edm::ParameterSet& iConfig)
    : geometry_getToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
      simDoublets_getToken_(consumes(iConfig.getParameter<edm::InputTag>("simDoubletsSrc"))),
      cellMinz_(iConfig.getParameter<std::vector<double>>("cellMinz")),
      cellMaxz_(iConfig.getParameter<std::vector<double>>("cellMaxz")),
      cellPhiCuts_(iConfig.getParameter<std::vector<int>>("cellPhiCuts")),
      cellMaxr_(iConfig.getParameter<std::vector<double>>("cellMaxr")),
      cellMinYSizeB1_(iConfig.getParameter<int>("cellMinYSizeB1")),
      cellMinYSizeB2_(iConfig.getParameter<int>("cellMinYSizeB2")),
      cellMaxDYSize12_(iConfig.getParameter<int>("cellMaxDYSize12")),
      cellMaxDYSize_(iConfig.getParameter<int>("cellMaxDYSize")),
      cellMaxDYPred_(iConfig.getParameter<int>("cellMaxDYPred")),
      cellZ0Cut_(iConfig.getParameter<double>("cellZ0Cut")),
      cellPtCut_(iConfig.getParameter<double>("cellPtCut")),
      folder_(iConfig.getParameter<std::string>("folder")) {
  hVector_dr_.resize(numLayerPairs);
  hVector_dphi_.resize(numLayerPairs);
  hVector_idphi_.resize(numLayerPairs);
  hVector_innerZ_.resize(numLayerPairs);
  hVector_clusterSizeY_.resize(numLayerPairs);
  hVector_dsizeYonlyBarrel_.resize(numLayerPairs);
  hVector_dsizeYinnerBarrel_.resize(numLayerPairs);
}

SimDoubletsAnalyzer::~SimDoubletsAnalyzer() {}

// -----------------------
// member functions
// -----------------------

void SimDoubletsAnalyzer::dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  trackerGeometry_ = &iSetup.getData(geometry_getToken_);
}

void SimDoubletsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  eventCount_++;

  // get simDoublets
  SimDoubletsCollection const& simDoubletsCollection = iEvent.get(simDoublets_getToken_);

  // loop over SimDoublets (= loop over TrackingParticles)
  for (auto const& simDoublets : simDoubletsCollection) {
    // get true pT of the TrackingParticle
    auto true_pT = simDoublets.trackingParticle()->pt();
    auto true_eta = simDoublets.trackingParticle()->eta();

    // create the true RecHit doublets of the TrackingParticle
    auto doublets = simDoublets.getSimDoublets(trackerGeometry_);

    // loop over those doublets
    for (auto const& doublet : doublets) {
      // RecHit properties
      auto inner_r = doublet.innerGlobalPos().perp();
      auto inner_z = doublet.innerGlobalPos().z();
      auto inner_phi = doublet.innerGlobalPos().barePhi();  // returns float, whereas .phi() returns phi object
      auto inner_iphi = unsafe_atan2s<7>(doublet.innerGlobalPos().y(), doublet.innerGlobalPos().x());
      auto outer_r = doublet.outerGlobalPos().perp();
      auto outer_z = doublet.outerGlobalPos().z();
      auto outer_phi = doublet.outerGlobalPos().barePhi();
      auto outer_iphi = unsafe_atan2s<7>(doublet.outerGlobalPos().y(), doublet.outerGlobalPos().x());

      auto dz = outer_z - inner_z;
      auto dr = outer_r - inner_r;
      auto dphi = reco::deltaPhi(inner_phi, outer_phi);
      auto idphi = std::min(std::abs(int16_t(outer_iphi - inner_iphi)), std::abs(int16_t(inner_iphi - outer_iphi)));

      // ----------------------------------------------------------
      // layer pair independent plots (main folder)
      // ----------------------------------------------------------

      // outer layer vs inner layer of SimDoublets
      h_layerPairs_->Fill(doublet.innerLayerId(), doublet.outerLayerId());

      // number of skipped layers by SimDoublets
      h_numSkippedLayers_->Fill(doublet.numSkippedLayers());

      // longitudinal impact parameter with respect to the beamspot
      double z0 = std::abs(inner_r * outer_z - inner_z * outer_r) / dr;
      h_z0_->Fill(z0);

      // radius of the circle defined by the two RecHits and the beamspot
      auto curvature = 1.f / 2.f * std::sqrt((dr / dphi) * (dr / dphi) + (inner_r * outer_r));
      h_curvatureR_->Fill(curvature);

      // pT that this curvature radius corresponds to
      auto pT = curvature / 87.78f;
      h_pTFromR_->Fill(pT);

      // ----------------------------------------------------------
      // layer pair dependent plots (sub-folders for layer pairs)
      // ----------------------------------------------------------

      // first, get layer pair Id and exclude layer pairs that are not considered
      int layerPairId = doublet.layerPairId();
      if (layerPairId2Index.find(layerPairId) == layerPairId2Index.end()) {
        continue;
      }

      // get the position of the layer pair in the vectors of histograms
      int layerPairIdIndex = layerPairId2Index.at(layerPairId);

      // dr = (outer_r - inner_r) histogram
      hVector_dr_[layerPairIdIndex]->Fill(dr);

      // dphi histogram
      hVector_dphi_[layerPairIdIndex]->Fill(dphi);
      hVector_idphi_[layerPairIdIndex]->Fill(idphi);

      // z of the inner RecHit histogram
      hVector_innerZ_[layerPairIdIndex]->Fill(inner_z);

      // cluster size in local y histogram
      auto innerClusterSizeY = doublet.innerRecHit()->cluster()->sizeY();
      hVector_clusterSizeY_[layerPairIdIndex]->Fill(innerClusterSizeY);

      // create bool that indicates if the doublet gets cut
      bool doubletGetsCut = false;
      // apply all cuts that do not depend on the cluster size
      // z window cut
      if (inner_z < cellMinz_[layerPairIdIndex] || inner_z > cellMaxz_[layerPairIdIndex]) {
        doubletGetsCut = true;
      }
      // z0cutoff
      if (dr > cellMaxr_[layerPairIdIndex] || dr < 0 || z0 > cellZ0Cut_) {
        doubletGetsCut = true;
      }
      // ptcut
      if (pT < cellPtCut_) {
        doubletGetsCut = true;
      }
      // iphicut
      if (idphi > cellPhiCuts_[layerPairIdIndex]) {
        doubletGetsCut = true;
      }

      // determine the moduleId
      DetId detIdObject(doublet.innerRecHit()->geographicalId());
      const GeomDetUnit* geomDetUnit = trackerGeometry_->idToDetUnit(detIdObject);
      const uint32_t moduleId = geomDetUnit->index();

      // define bools needed to decide on cutting parameters
      const bool innerInB1 = (doublet.innerLayerId() == 0);
      const bool innerInB2 = (doublet.innerLayerId() == 1);
      const bool isOuterLadder = (0 == (moduleId / 8) % 2);  // check if this even makes sense in Phase-2
      const bool innerInBarrel = (doublet.innerLayerId() < 4);
      const bool outerInBarrel = (doublet.outerLayerId() < 4);

      // histograms for clusterCut
      // cluster size in local y
      if (!outerInBarrel) {
        if (innerInB1 && isOuterLadder) {
          h_YsizeB1_->Fill(innerClusterSizeY);
          // apply the cut
          if (innerClusterSizeY < cellMinYSizeB1_) {
            doubletGetsCut = true;
          }
        }
        if (innerInB2) {
          h_YsizeB2_->Fill(innerClusterSizeY);
          // apply the cut
          if (innerClusterSizeY < cellMinYSizeB2_) {
            doubletGetsCut = true;
          }
        }
      }

      // histograms for zSizeCut
      if (innerInBarrel) {
        if (outerInBarrel) {  // onlyBarrel
          auto DYsize = std::abs(innerClusterSizeY - doublet.outerRecHit()->cluster()->sizeY());
          if (innerInB1 && isOuterLadder) {
            hVector_dsizeYonlyBarrel_[layerPairIdIndex]->Fill(DYsize);
            h_DYsize12_->Fill(DYsize);
            // apply the cut
            if (DYsize > cellMaxDYSize12_) {
              doubletGetsCut = true;
            }
          } else if (!innerInB1) {
            hVector_dsizeYonlyBarrel_[layerPairIdIndex]->Fill(DYsize);
            h_DYsize_->Fill(DYsize);
            // apply the cut
            if (DYsize > cellMaxDYSize_) {
              doubletGetsCut = true;
            }
          }
        } else {  // not onlyBarrel
          int DYsizePred =
              std::abs(innerClusterSizeY - int(std::abs(dz / dr) * pixelTopology::Phase2::dzdrFact + 0.5f));
          hVector_dsizeYinnerBarrel_[layerPairIdIndex]->Fill(DYsizePred);
          h_DYPred_->Fill(DYsizePred);
          // apply the cut
          if (DYsizePred > cellMaxDYPred_) {
            doubletGetsCut = true;
          }
        }
      }

      // fill the number histograms
      // histogram of all valid doublets
      h_numTotVsPt_->Fill(true_pT);
      h_numTotVsEta_->Fill(true_eta);
      // fill histogram of doublets that pass all cuts
      if (!doubletGetsCut) {
        h_numPassVsPt_->Fill(true_pT);
        h_numPassVsEta_->Fill(true_eta);
      }
    }  // end loop over those doublets
  }  // end loop over SimDoublets (= loop over TrackingParticles)
}

void SimDoubletsAnalyzer::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  // set some common parameters
  int pTNBins = 50;
  double pTmin = log10(0.01);
  double pTmax = log10(1000);
  int etaNBins = 50;
  double etamin = -4.;
  double etamax = 4.;

  ibook.setCurrentFolder(folder_);

  // ----------------------------------------------------------
  // booking layer pair independent histograms (main folder)
  // ----------------------------------------------------------

  // overview histograms
  h_layerPairs_ = ibook.book2D(
      "layerPairs", "Layer pairs in SimDoublets; Inner layer ID; Outer layer ID", 28, -0.5, 27.5, 28, -0.5, 27.5);
  h_numSkippedLayers_ = ibook.book1D(
      "numSkippedLayers", "Number of skipped layers; Number of skipped layers; Number of SimDoublets", 16, -1.5, 14.5);
  h_numTotVsPt_ = simdoublets::make1DLogX(
      ibook,
      "numTotVsPt",
      "Total number of SimDoublets; True transverse momentum p_{T} [GeV]; Total number of valid SimDoublets",
      pTNBins,
      pTmin,
      pTmax);
  h_numPassVsPt_ = simdoublets::make1DLogX(ibook,
                                           "numPassVsPt",
                                           "Number of passing SimDoublets; True transverse momentum p_{T} [GeV]; "
                                           "Number of valid SimDoublets passing all cuts",
                                           pTNBins,
                                           pTmin,
                                           pTmax);
  h_numTotVsEta_ =
      ibook.book1D("numTotVsEta",
                   "Total number of SimDoublets; True pseudorapidity #eta; Total number of valid SimDoublets",
                   etaNBins,
                   etamin,
                   etamax);
  h_numPassVsEta_ = ibook.book1D(
      "numPassVsEta",
      "Total number of SimDoublets; True pseudorapidity #eta; Number of valid SimDoublets passing all cuts",
      etaNBins,
      etamin,
      etamax);

  // histogram for z0cutoff  (z0Cut)
  h_z0_ = ibook.book1D("z0", "z_{0}; Longitudinal impact parameter z_{0} [cm]; Number of SimDoublets", 51, -1, 50);

  // histograms for ptcut  (ptCut)
  h_curvatureR_ = ibook.book1D(
      "curvatureR", "Curvature from SimDoublet+beamspot; Curvature radius [cm] ; Number of SimDoublets", 100, 0, 1000);
  h_pTFromR_ = simdoublets::make1DLogX(
      ibook,
      "pTFromR",
      "Transverse momentum from curvature; Transverse momentum p_{T} [GeV]; Number of SimDoublets",
      pTNBins,
      pTmin,
      pTmax);

  // histograms for clusterCut  (minYsizeB1 and minYsizeB2)
  h_YsizeB1_ =
      ibook.book1D("YsizeB1", "Cluster size Y (inner from B1); Cluster Size Y; Number of SimDoublets", 51, -1, 50);
  h_YsizeB2_ =
      ibook.book1D("YsizeB2", "Cluster size Y (inner from B2); Cluster Size Y; Number of SimDoublets", 51, -1, 50);

  // histograms for zSizeCut  (maxDYsize12, maxDYsize and maxDYPred)
  h_DYsize12_ = ibook.book1D("DYsize12",
                             "Difference in cluster y-size (inner from B1); Absolute difference in cluster y-size of "
                             "the two RecHits; Number of SimDoublets",
                             31,
                             -1,
                             30);
  h_DYsize_ = ibook.book1D(
      "DYsize",
      "Difference in cluster y-size; Absolute difference in cluster y-size of the two RecHits; Number of SimDoublets",
      31,
      -1,
      30);
  h_DYPred_ = ibook.book1D("DYPred",
                           "Projected difference in cluster z-size; Absolute difference in projected cluster z-size of "
                           "the two RecHits; Number of SimDoublets",
                           201,
                           -1,
                           200);

  // -----------------------------------------------------------------------
  // booking layer pair dependent histograms (sub-folders for layer pairs)
  // -----------------------------------------------------------------------

  // loop through valid layer pairs and add for each one booked hist per vector
  for (auto id = layerPairId2Index.begin(); id != layerPairId2Index.end(); ++id) {
    // get the position of the layer pair in the histogram vectors
    int layerPairIdIndex = id->second;

    // get layer names from the layer pair Id
    auto layerNames = simdoublets::getInnerOuterLayerNames(id->first);
    std::string innerLayerName = layerNames.first;
    std::string outerLayerName = layerNames.second;

    // name the sub-folder for the layer pair "lp_${innerLayerId}_${outerLayerId}"
    std::string subFolderName = "/lp_" + innerLayerName + "_" + outerLayerName;

    // layer mentioning in histogram titles
    std::string layerTitle = "(layers (" + innerLayerName + "," + outerLayerName + "))";

    // set folder to the sub-folder for the layer pair
    ibook.setCurrentFolder(folder_ + subFolderName);

    // histogram for z0cutoff  (maxr)
    hVector_dr_.at(layerPairIdIndex) = ibook.book1D(
        "dr",
        "dr of RecHit pair " + layerTitle + "; dr between outer and inner RecHit [cm]; Number of SimDoublets",
        31,
        -1,
        30);

    // histograms for iphicut  (phiCuts)
    hVector_dphi_.at(layerPairIdIndex) = ibook.book1D(
        "dphi",
        "dphi of RecHit pair " + layerTitle + "; d#phi between outer and inner RecHit [rad]; Number of SimDoublets",
        50,
        -M_PI,
        M_PI);
    hVector_idphi_.at(layerPairIdIndex) =
        ibook.book1D("idphi",
                     "idphi of RecHit pair " + layerTitle +
                         "; Absolute int d#phi between outer and inner RecHit [rad]; Number of SimDoublets",
                     50,
                     0,
                     1000);

    // histogram for z window  (minz and maxz)
    hVector_innerZ_.at(layerPairIdIndex) =
        ibook.book1D("innerZ",
                     "z of the inner RecHit " + layerTitle + "; z of inner RecHit [cm]; Number of SimDoublets",
                     100,
                     -300,
                     300);

    // other histograms
    hVector_dsizeYonlyBarrel_.at(layerPairIdIndex) =
        ibook.book1D("dsizeYonlyBarrel",
                     "Cluster Size Y Difference between outer and inner RecHit " + layerTitle +
                         "; Cluster Size Y Difference ; Number of SimDoublets",
                     51,
                     -1,
                     50);
    hVector_dsizeYinnerBarrel_.at(layerPairIdIndex) =
        ibook.book1D("dsizeYinnerBarrel",
                     "Projected Cluster Size Y Difference between outer and inner RecHit " + layerTitle +
                         "; Cluster Size Y Difference ; Number of SimDoublets",
                     51,
                     -1,
                     50);
    hVector_clusterSizeY_.at(layerPairIdIndex) =
        ibook.book1D("sizeY", "Cluster Size Y " + layerTitle + "; Cluster Size Y ; Number of SimDoublets", 51, -1, 50);
  }
}

void SimDoubletsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("folder", "Tracking/TrackingMCTruth/SimDoublets");
  desc.add<edm::InputTag>("simDoubletsSrc", edm::InputTag("simDoubletsProducer"));

  // cutting parameters
  desc.add<std::vector<double>>("cellMinz", std::vector<double>(55, -20.))->setComment("Minimum z for each layer pair");
  desc.add<std::vector<double>>("cellMaxz", std::vector<double>(55, 20.))->setComment("Maximum z for each layer pair");
  desc.add<std::vector<int>>("cellPhiCuts", std::vector<int>(55, 20))->setComment("Cuts in phi for cells");
  desc.add<std::vector<double>>("cellMaxr", std::vector<double>(55, 20.))->setComment("Cut for dr of cells");
  desc.add<int>("cellMinYSizeB1", 25)->setComment("Minimum cluster size for B1");
  desc.add<int>("cellMinYSizeB2", 15)->setComment("Minimum cluster size for B2");
  desc.add<int>("cellMaxDYSize12", 12)->setComment("Maximum cluster size difference for B1/B2");
  desc.add<int>("cellMaxDYSize", 10)->setComment("Maximum cluster size difference");
  desc.add<int>("cellMaxDYPred", 20)->setComment("Maximum cluster size difference prediction");
  desc.add<double>("cellZ0Cut", 7.5)->setComment("Maximum longitudinal impact parameter");
  desc.add<double>("cellPtCut", 0.85)->setComment("Minimum tranverse momentum");

  descriptions.addWithDefaultLabel(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(SimDoubletsAnalyzer);
