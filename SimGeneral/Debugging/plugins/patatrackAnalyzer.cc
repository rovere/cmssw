// -*- C++ -*-
//
// Package:    patatrackAnalyzer/patatrackAnalyzer
// Class:      patatrackAnalyzer
//
/**\class patatrackAnalyzer patatrackAnalyzer.cc patatrackAnalyzer/patatrackAnalyzer/plugins/patatrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andreas Gruber
//         Created:  Mon, 16 Oct 2023 14:24:35 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "genTau.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class patatrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit patatrackAnalyzer(const edm::ParameterSet&);
  ~patatrackAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void buildSimTau (int, const reco::GenParticle &, const std::vector<CaloParticle> &);
  const edm::EDGetTokenT<std::vector<CaloParticle>> CaloParticle_token;
  const edm::EDGetTokenT<std::vector<SimCluster>> SimClusters_token;
  const edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;
  const edm::EDGetTokenT<std::vector<int>> genBarcodesToken_;

  int recursion_n;
  int daughter_n;
  edm::Ref<reco::CaloParticleCollection> matchedCPs;
  // ----------member data ---------------------------
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
patatrackAnalyzer::patatrackAnalyzer(const edm::ParameterSet& iConfig)
    : CaloParticle_token(consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("CaloParticle"))),
    SimClusters_token(consumes<std::vector<SimCluster>>(iConfig.getParameter<edm::InputTag>("SimClusters"))),
    genParticles_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticles"))),
    genBarcodesToken_(mayConsume<std::vector<int>>(iConfig.getParameter<edm::InputTag>("genParticles")))
    {

    edm::Service<TFileService> fs;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

patatrackAnalyzer::~patatrackAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

void patatrackAnalyzer::buildSimTau (int n, const reco::GenParticle &genPart, const std::vector<CaloParticle> &CaloPartVec) {
    auto &daughters = genPart.daughterRefVector();
    //if(daughters.empty()){
    n++;
    //}

    for (auto daughter = daughters.begin(); daughter != daughters.end(); ++daughter) {
        if (abs((*daughter)->pdgId()) == 15) {
            break;
        }
        std::cout << "recursion #: " << n << std::endl;
        std::cout << (*daughter)->pdgId() << std::endl;
        //recursion_n++;
        buildSimTau(n, *(*daughter), CaloPartVec);
        //std::cout << "Printing daughters of" << ge
    }
}

/*edm::Ref<reco::GenParticleCollection> patatrackAnalyzer::generatorRef_(const SimTrack &simtrk, const GlobalContext &g) const {
    assert(simtrk.genpartIndex() != -1);
    // Note that simtrk.genpartIndex() is the barcode, not the index within GenParticleCollection, so I have to search the particle
    std::vector<int>::const_iterator it;
    if (gg.barcodesAreSorted) {
        it = std::lower_bound(g.genBarcodes->begin(), g.genBarcodes->end(), simtrk.genpartIndex());
    } else {
        it = std::find(g.genBarcodes->begin(), g.genBarcodes->end(), simtrk.genpartIndex());
    }
    // Check that I found something
    // I need to check '*it == simtrk.genpartIndex()' because lower_bound just finds the right spot for an item in a sorted list, not the item
    if ((it != g.genBarcodes->end()) && (*it == simtrk.genpartIndex())) {
        return reco::GenParticleRef(g.gens, it - g.genBarcodes->begin());
    } else {
        return reco::GenParticleRef();
    }
}

int patatrackAnalyzer::findGenPartIdx (const reco::GenParticle &c) {
   
}*/
//
// member functions
//

// ------------ method called for each event  ------------
void patatrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    edm::Handle<std::vector<CaloParticle>> CaloParticle_h;
    iEvent.getByToken(CaloParticle_token, CaloParticle_h);
    const auto& CaloParticle = *CaloParticle_h;   
    edm::Handle<std::vector<SimCluster>> SimClusters_h;
    iEvent.getByToken(SimClusters_token, SimClusters_h);
    const auto& SimClusters = *SimClusters_h;    
    edm::Handle<std::vector<reco::GenParticle>> gen_particles_h;
    iEvent.getByToken(genParticles_, gen_particles_h);
    const auto& genParticles = *gen_particles_h;  
    edm::Handle<std::vector<int>> gen_barcodes_h;
    iEvent.getByToken(genBarcodesToken_, gen_barcodes_h);
    const auto& genBarcodes = *gen_barcodes_h;  

    genTau test;
    std::cout << test.getSize() << std::endl;

    /*for (const auto& genBarcode: genBarcodes) {
        std::cout << "GenParticle Barcode: " << genBarcode << std::endl;
        
                auto &daughters = genParticle.daughterRefVector();
            for (auto daughter = daughters.begin(); daughter != daughters.end(); ++daughter) {
                daughter.genParticleId();
    
    }*/
    int tau_n=0;
    for (const auto& GP: genParticles) {
        if (abs(GP.pdgId()) == 15) {
            recursion_n=0;
            std::cout << "Found Tau #" << tau_n << "!" << std::endl;
            buildSimTau(tau_n, GP, CaloParticle);
            tau_n++;
            }
        }
    
    /*    for (const auto& CP: CaloParticle) {
            std::cout << "SimTrack size: " << CP.g4Tracks().size() << std::endl;
            std::cout << "SimTrack GenPartIndex: " << CP.g4Tracks()[0].genpartIndex() << std::endl;
            assert(CP.g4Tracks()[0].genpartIndex() != -1);
            // Note that simtrk.genpartIndex() is the barcode, not the index within GenParticleCollection, so I have to search the particle
            std::vector<int>::const_iterator it;
            it = std::find(genBarcodes.begin(), genBarcodes.end(), CP.g4Tracks()[0].genpartIndex());
            if ((it != genBarcodes.end())) {
                auto GP = reco::GenParticleRef(gen_particles_h, it - genBarcodes.begin());
                auto genPartReturn = *GP;
                std::cout << genPartReturn.pdgId() << std::endl;

            } else {
                std::cout << "no genParticle found" << std::endl;
            }
        }
    }

    
            int simTrackGenIndex = CP.g4Tracks()[0].genpartIndex();
            std::cout << "simTrack Gen Index: " << simTrackGenIndex << std::endl;
            std::cout << "CP pdgId: " << CP.pdgId() << std::endl;
            if (simTrackGenIndex < static_cast<int>(genParticles.size())) {
                findGenTau(genParticles[simTrackGenIndex]);
                if (CP.g4Tracks()[0].genpartIndex() == 15) {*/

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void patatrackAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void patatrackAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void patatrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(patatrackAnalyzer);
