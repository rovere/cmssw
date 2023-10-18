// -*- C++ -*-
//
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

#define DEBUG 0

struct TauDecay {
  std::vector<std::pair<int, int>> resonances;
  std::vector<std::pair<int, int>> leaves;

  void dump(void) {
    for (auto const & l : leaves) {
      std::cout << "L " << l.first << " " << l.second << std::endl;
    }
    for (auto const & r : resonances) {
      std::cout << "R " << r.first << " " << r.second << std::endl;
    }
  }

  void dumpDecay(const std::pair<int, int> & entry) {
    if (entry.second == -1) { // No intermediate mother.
      std::cout << entry.first << " " << entry.second << std::endl;;
    } else {
      std::cout << entry.first << " " << entry.second << " coming from: ";
      auto const &mother = resonances[entry.second];
      dumpDecay(mother);
    }
  }

  void dumpFullDecay(void) {
    for (auto const & leaf : leaves) {
      dumpDecay(leaf);
    }
  }
};

class SimTauAnalyzer : public edm::one::EDAnalyzer<> {
public:

  explicit SimTauAnalyzer(const edm::ParameterSet&);
  ~SimTauAnalyzer() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void buildSimTau(TauDecay &, uint8_t, int, const reco::GenParticle &, const std::vector<CaloParticle> &);
  // ----------member data ---------------------------
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloParticle_token_;
  const edm::EDGetTokenT<std::vector<SimCluster>> simClusters_token_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token_;
  const edm::EDGetTokenT<std::vector<int>> genBarcodes_token_;
  //edm::Ref<reco::CaloParticleCollection> matchedCPs;
};

//
// static data member definitions
//

//
// constructors and destructor
//
SimTauAnalyzer::SimTauAnalyzer(const edm::ParameterSet& iConfig)
    : caloParticle_token_(consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("CaloParticle"))),
    simClusters_token_(consumes<std::vector<SimCluster>>(iConfig.getParameter<edm::InputTag>("SimClusters"))),
    genParticles_token_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticles"))),
    genBarcodes_token_(consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("genParticles"))) {}


void SimTauAnalyzer::buildSimTau(TauDecay &t,
    uint8_t generation,
    int resonance_idx,
    const reco::GenParticle & gen_particle, // the decaying tau,
    const std::vector<CaloParticle> &CaloPartVec) {

  auto &daughters = gen_particle.daughterRefVector();
  bool is_leaf = (daughters.size() == 0);
  if (is_leaf) {
    if (DEBUG)
      std::cout << " TO BE SAVED " << resonance_idx << " ";
    t.leaves.push_back({gen_particle.pdgId(), resonance_idx});
    return;
  } else if (generation !=0) {
    t.resonances.push_back({gen_particle.pdgId(), resonance_idx});
    resonance_idx = t.resonances.size() - 1;
    if (DEBUG)
      std::cout << " RESONANCE/ITERMEDIATE " << resonance_idx << " ";
  }
  if (DEBUG)
    std::cout << std::endl;

  ++generation;
  std::string separator("");
  for (uint8_t c = 0; c < generation; ++c) {
    separator += "  ";
  }
  for (auto daughter = daughters.begin(); daughter != daughters.end(); ++daughter) {
    auto const & daughter_flags = (*daughter)->statusFlags();
    if (DEBUG) {
      std::cout << separator << " gen " << (int)generation << " " << (*daughter)->pdgId() << " ";
      for (unsigned int bit = 0; bit <= reco::GenStatusFlags::kIsLastCopyBeforeFSR; ++bit) {
        std::cout << daughter_flags.flags_[bit] << " ";
      }
    }
    buildSimTau(t, generation, resonance_idx, *(*daughter), CaloPartVec);
    if (DEBUG)
      std::cout << std::endl;
  }
}

/*edm::Ref<reco::GenParticleCollection> SimTauAnalyzer::generatorRef_(const SimTrack &simtrk, const GlobalContext &g) const {
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

int SimTauAnalyzer::findGenPartIdx (const reco::GenParticle &c) {

}*/
//
// member functions
//

// ------------ method called for each event  ------------
void SimTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<std::vector<CaloParticle>> CaloParticle_h;
  iEvent.getByToken(caloParticle_token_, CaloParticle_h);
  //    Handle<std::vector<SimCluster>> SimClusters_h;
  //    iEvent.getByToken(simClusters_token_, SimClusters_h);
  edm::Handle<std::vector<reco::GenParticle>> gen_particles_h;
  iEvent.getByToken(genParticles_token_, gen_particles_h);
  //    Handle<std::vector<int>> gen_barcodes_h;
  //    iEvent.getByToken(genBarcodes_token_, gen_barcodes_h);

  const auto& caloParticle = *CaloParticle_h;
  //    const auto& SimClusters = *SimClusters_h;
  const auto& genParticles = *gen_particles_h;
  //    const auto& genBarcodes = *gen_barcodes_h;

  int counter=0;
  for (auto const & g : genParticles) {
    auto const & flags = g.statusFlags();
    if (std::abs(g.pdgId()) == 15 and flags.isPrompt() and flags.isDecayedLeptonHadron()) {
      if (DEBUG) {
        std::cout << iEvent.eventAuxiliary().event() << " " << counter++ << " "
          << " " << g.pdgId() << " ";
        for (unsigned int bit = 0; bit <= reco::GenStatusFlags::kIsLastCopyBeforeFSR; ++bit) {
          std::cout << flags.flags_[bit] << " ";
        }
        std::cout << std::endl;
      }
      TauDecay t;
      buildSimTau(t, 0, -1, g, caloParticle);
      t.dumpFullDecay();
      //        if (std::abs(g.pdgId()) == 15) {
      //          auto &daughters = g.daughterRefVector();
      //          for (auto daughter = daughters.begin(); daughter != daughters.end(); ++daughter) {
      //            auto const & daughter_flags = (*daughter)->statusFlags();
      //            std::cout << "\t\t" << (*daughter)->pdgId() << " ";
      //            for (unsigned int bit = 0; bit <= reco::GenStatusFlags::kIsLastCopyBeforeFSR; ++bit) {
      //              std::cout << daughter_flags.flags_[bit] << " ";
      //            }
      //            std::cout << std::endl;
      //          }
      //        }
    }
  }
}
//    /*for (const auto& genBarcode: genBarcodes) {
//        std::cout << "GenParticle Barcode: " << genBarcode << std::endl;
//
//                auto &daughters = genParticle.daughterRefVector();
//            for (auto daughter = daughters.begin(); daughter != daughters.end(); ++daughter) {
//                daughter.genParticleId();
//
//    }*/
//    int generation=0;
//  for (const auto& gen_particle: gen_particles) {
//    if (abs(gen_particle.pdgId()) == 15) {
//
//      // if 15 is in the list of dauthers, it means that the original
//      // Tau went through some brem or similar process. We ignore all the
//      // intermediate products (dauthers_ and just keep on following the tau)
//      auto const &daughters = gen_particles.daughterRefVector();
//      auto intermediate_tau = std::find(daughters.begin(), daughters.end(), [](const auto &d) {return d.pdgId()==15);
//        if (intermediate_tau != daughters.end()) {
//    for (const auto& gen_particle: genParticles) {
//        if (abs(gen_particle.pdgId()) == 15) {
//            recursion_n=0;
//            std::cout << "Found Tau #" << tau_n << "!" << std::endl;
//            buildSimTau(generation, gen_particle, CaloParticle);
//            tau_n++;
//            }
//        }
//
//    /*    for (const auto& CP: CaloParticle) {
//            std::cout << "SimTrack size: " << CP.g4Tracks().size() << std::endl;
//            std::cout << "SimTrack GenPartIndex: " << CP.g4Tracks()[0].genpartIndex() << std::endl;
//            assert(CP.g4Tracks()[0].genpartIndex() != -1);
//            // Note that simtrk.genpartIndex() is the barcode, not the index within GenParticleCollection, so I have to search the particle
//            std::vector<int>::const_iterator it;
//            it = std::find(genBarcodes.begin(), genBarcodes.end(), CP.g4Tracks()[0].genpartIndex());
//            if ((it != genBarcodes.end())) {
//                auto gen_particle = reco::GenParticleRef(gen_particles_h, it - genBarcodes.begin());
//                auto genPartReturn = *gen_particle;
//                std::cout << genPartReturn.pdgId() << std::endl;
//
//            } else {
//                std::cout << "no genParticle found" << std::endl;
//            }
//        }
//    }
//
//
//            int simTrackGenIndex = CP.g4Tracks()[0].genpartIndex();
//            std::cout << "simTrack Gen Index: " << simTrackGenIndex << std::endl;
//            std::cout << "CP pdgId: " << CP.pdgId() << std::endl;
//            if (simTrackGenIndex < static_cast<int>(genParticles.size())) {
//                findGenTau(genParticles[simTrackGenIndex]);
//                if (CP.g4Tracks()[0].genpartIndex() == 15) {*/
//
//}

// ------------ method called once each job just before starting event loop  ------------
void SimTauAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void SimTauAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SimTauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(SimTauAnalyzer);
