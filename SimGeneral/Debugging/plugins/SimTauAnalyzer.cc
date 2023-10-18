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

#define DEBUG 1

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

  void buildSimTau(TauDecay &, uint8_t, int, const reco::GenParticle &, int, const std::vector<CaloParticle> &, const std::vector<int> &);
  // ----------member data ---------------------------
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloParticle_token_;
  const edm::EDGetTokenT<std::vector<SimCluster>> simClusters_token_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token_;
  const edm::EDGetTokenT<std::vector<int>> genBarcodes_token_;
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
    const reco::GenParticle & gen_particle,
    int gen_particle_key,
    const std::vector<CaloParticle> &caloPartVec,
    const std::vector<int> & gen_particle_barcodes) {

  auto &daughters = gen_particle.daughterRefVector();
  bool is_leaf = (daughters.size() == 0);
  if (is_leaf) {
    if (DEBUG)
      std::cout << " TO BE SAVED " << resonance_idx << " ";
    t.leaves.push_back({gen_particle.pdgId(), resonance_idx});
    auto const & gen_particle_barcode = gen_particle_barcodes[gen_particle_key];
    auto const & found_in_caloparticles = std::find_if(
        caloPartVec.begin(),
        caloPartVec.end(),
        [&](const auto &p) {
        return p.g4Tracks()[0].genpartIndex() == gen_particle_barcode;}
        );
    if (found_in_caloparticles != caloPartVec.end()) {
      auto calo_particle_idx = (found_in_caloparticles - caloPartVec.begin());
      if (DEBUG)
        std::cout << " CP " << calo_particle_idx << " " << caloPartVec[calo_particle_idx];
    }
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
    int gen_particle_key = (*daughter).key();
    if (DEBUG) {
      std::cout << separator << " gen " << (int)generation << " " << gen_particle_key << " " << (*daughter)->pdgId() << " ";
      for (unsigned int bit = 0; bit <= reco::GenStatusFlags::kIsLastCopyBeforeFSR; ++bit) {
        std::cout << daughter_flags.flags_[bit] << " ";
      }
    }
    buildSimTau(t, generation, resonance_idx, *(*daughter), gen_particle_key, caloPartVec, gen_particle_barcodes);
    if (DEBUG)
      std::cout << std::endl;
  }
}

// ------------ method called for each event  ------------
void SimTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<std::vector<CaloParticle>> CaloParticle_h;
  iEvent.getByToken(caloParticle_token_, CaloParticle_h);
  edm::Handle<std::vector<reco::GenParticle>> gen_particles_h;
  iEvent.getByToken(genParticles_token_, gen_particles_h);
  Handle<std::vector<int>> gen_barcodes_h;
  iEvent.getByToken(genBarcodes_token_, gen_barcodes_h);

  const auto& caloParticle = *CaloParticle_h;
  const auto& genParticles = *gen_particles_h;
  const auto& genBarcodes = *gen_barcodes_h;

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
      buildSimTau(t, 0, -1, g, -1, caloParticle, genBarcodes);
      t.dumpFullDecay();
    }
  }
}

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
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimTauAnalyzer);
