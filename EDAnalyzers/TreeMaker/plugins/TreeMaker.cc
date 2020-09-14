// -*- C++ -*-
//
// Package:    EDAnalyzers/TreeMaker
// Class:      TreeMaker
//
/**\class TreeMaker TreeMaker.cc EDAnalyzers/TreeMaker/plugins/TreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Sat, 11 May 2019 13:14:55 GMT
//
//


// system include files
# include <memory>

// user include files



# include "CommonTools/UtilAlgos/interface/TFileService.h"
# include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
# include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
# include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
# include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
# include "DataFormats/FWLite/interface/ESHandle.h"
# include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
# include "DataFormats/HGCalReco/interface/Trackster.h"
# include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
# include "DataFormats/HepMCCandidate/interface/GenParticle.h"
# include "DataFormats/JetReco/interface/PFJet.h"
# include "DataFormats/Math/interface/LorentzVector.h"
# include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
# include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
# include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
# include "DataFormats/TrackReco/interface/Track.h"
# include "DataFormats/TrackReco/interface/TrackFwd.h"
# include "DataFormats/VertexReco/interface/Vertex.h"
# include "FWCore/Framework/interface/Event.h"
# include "FWCore/Framework/interface/ESHandle.h"
# include "FWCore/Framework/interface/Frameworkfwd.h"
# include "FWCore/Framework/interface/MakerMacros.h"
# include "FWCore/Framework/interface/one/EDAnalyzer.h"
# include "FWCore/ParameterSet/interface/ParameterSet.h"
# include "FWCore/ServiceRegistry/interface/Service.h"
# include "FWCore/Utilities/interface/InputTag.h"
# include "Geometry/CaloTopology/interface/HGCalTopology.h"
# include "Geometry/Records/interface/IdealGeometryRecord.h"
# include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
# include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
# include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
# include "SimDataFormats/CaloHit/interface/PCaloHit.h"
# include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
# include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

# include "EDAnalyzers/TreeMaker/interface/Common.h"
# include "EDAnalyzers/TreeMaker/interface/Constants.h"
# include "EDAnalyzers/TreeMaker/interface/TreeOutputInfo.h"

# include <CLHEP/Matrix/Matrix.h>
# include <CLHEP/Vector/ThreeVector.h>
# include <CLHEP/Vector/ThreeVector.h>

# include <Compression.h>
# include <TH1F.h>
# include <TH2F.h>
# include <TMatrixD.h>
# include <TTree.h> 
# include <TVector2.h> 
# include <TVectorD.h> 

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



double HGCal_minEta = 1.479;
double HGCal_maxEta = 3.1;

double el_minPt = 10; //15;
double el_maxPt = 99999; //30;

double _largeVal = 999999999;


class TreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
    
    explicit TreeMaker(const edm::ParameterSet&);
    ~TreeMaker();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
    
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    
    double getDeltaPhi(double phi1, double phi2);
    
    std::tuple <TMatrixD, TMatrixD, TVectorD> getMultiClusterPC(
        reco::HGCalMultiCluster *multiCluster,
        std::map <DetId, const HGCRecHit*> m_recHit
    );
    
    
    int minLayer;
    int maxLayer;
    
    
    TreeOutputInfo::TreeOutput *treeOutput;
    
    
    // My stuff //
    bool debug;
    bool isGunSample;
    
    bool storeSimHit;
    bool storeRecHit;
    bool storeHGCALlayerClus;
    bool storeSuperClusTICLclus;
    
    double TICLeleGenMatchDR;
    
    
    // Gen particles //
    edm::EDGetTokenT <std::vector <reco::GenParticle> > tok_genParticle;
    
    
    // Pileup //
    edm::EDGetTokenT <std::vector <PileupSummaryInfo> > tok_pileup;
    
    
    // Rho //
    edm::EDGetTokenT <double> tok_rho;
    
    
    // HGCAL layer clusters //
    edm::EDGetTokenT <std::vector <reco::CaloCluster> > tok_HGCALlayerCluster;
    
    
    // TICL //
    edm::EDGetTokenT <std::vector <ticl::Trackster> > tok_TICLtrackster;
    edm::EDGetTokenT <std::vector <reco::HGCalMultiCluster> > tok_TICLmultiCluster;
    
    
    // Gsf electrons from TICL //
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > tok_gsfEleFromTICL;
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
TreeMaker::TreeMaker(const edm::ParameterSet& iConfig)
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    
    // Compression
    //fs->file().SetCompressionAlgorithm(ROOT::kLZMA);
    //fs->file().SetCompressionLevel(8);
    
    
    treeOutput = new TreeOutputInfo::TreeOutput("tree", fs);
    
    
    minLayer = +9999;
    maxLayer = -9999;
    
    
    // My stuff //
    debug = iConfig.getParameter <bool>("debug");
    isGunSample = iConfig.getParameter <bool>("isGunSample");
    
    storeSimHit = iConfig.getParameter <bool>("storeSimHit");
    storeRecHit = iConfig.getParameter <bool>("storeRecHit");
    storeHGCALlayerClus = iConfig.getParameter <bool>("storeHGCALlayerClus");
    storeSuperClusTICLclus = iConfig.getParameter <bool>("storeSuperClusTICLclus");
    
    TICLeleGenMatchDR = iConfig.getParameter <double>("TICLeleGenMatchDR");
    
    
    // Gen particles //
    tok_genParticle = consumes <std::vector <reco::GenParticle> >(iConfig.getUntrackedParameter <edm::InputTag>("label_genParticle"));
    
    
    // Pileup //
    tok_pileup = consumes <std::vector <PileupSummaryInfo> >(iConfig.getUntrackedParameter <edm::InputTag>("label_pileup"));
    
    
    // Rho //
    tok_rho = consumes <double>(iConfig.getUntrackedParameter <edm::InputTag>("label_rho"));
    
    
    // TICL //
    tok_TICLtrackster = consumes <std::vector <ticl::Trackster> >(iConfig.getUntrackedParameter <edm::InputTag>("label_TICLtrackster"));
    tok_TICLmultiCluster = consumes <std::vector <reco::HGCalMultiCluster> >(iConfig.getUntrackedParameter <edm::InputTag>("label_TICLmultiCluster"));
    
    
    // Gsf electrons from TICL //
    tok_gsfEleFromTICL = consumes <std::vector <reco::GsfElectron> >(iConfig.getUntrackedParameter <edm::InputTag>("label_gsfEleFromTICL"));
}


TreeMaker::~TreeMaker()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
    delete treeOutput;
}


//
// member functions
//


// ------------ method called for each event  ------------
void TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    long long eventNumber = iEvent.id().event();
    //printf("Event %llu \n", eventNumber);
    
    
    treeOutput->clear();
    
    //recHitTools.getEventSetup(iSetup);
    
    
    //////////////////// Run info ////////////////////
    treeOutput->runNumber = iEvent.id().run();
    treeOutput->eventNumber = iEvent.id().event();
    treeOutput->luminosityNumber = iEvent.id().luminosityBlock();
    treeOutput->bunchCrossingNumber = iEvent.bunchCrossing();
    
    
    //////////////////// Gen particle ////////////////////
    edm::Handle <std::vector <reco::GenParticle> > v_genParticle;
    iEvent.getByToken(tok_genParticle, v_genParticle);
    
    std::vector <CLHEP::HepLorentzVector> v_genEl_4mom;
    
    
    for(int iPart = 0; iPart < (int) v_genParticle->size(); iPart++)
    {
        reco::GenParticle part = v_genParticle->at(iPart);
        
        int pdgId = part.pdgId();
        int status = part.status();
        
        // Gen ele
        //if(abs(pdgId) == 11 && status == 1)
        //if(abs(pdgId) == 11 && (part.isHardProcess() || status == 1))
        if(
            abs(pdgId) == 11 && (
                (isGunSample && status == 1) ||
                (!isGunSample && part.isHardProcess())
            )
        )
        {
            //printf("[%llu] Gen electron found: E %0.2f, pT %0.2f, eta %+0.2f \n", eventNumber, part.energy(), part.pt(), part.eta());
            
            printf(
                "[%llu] "
                "Gen-ele found: E %0.2f, pT %0.2f, eta %+0.2f, pz %+0.2f, "
                "\n",
                eventNumber,
                part.energy(), part.pt(), part.eta(), part.pz()
            );
            
            if(fabs(part.eta()) > HGCal_minEta && fabs(part.eta()) < HGCal_maxEta && part.pt() > el_minPt && part.pt() < el_maxPt)
            {
                CLHEP::HepLorentzVector genEl_4mom;
                
                genEl_4mom.setT(part.energy());
                genEl_4mom.setX(part.px());
                genEl_4mom.setY(part.py());
                genEl_4mom.setZ(part.pz());
                
                v_genEl_4mom.push_back(genEl_4mom);
                
                treeOutput->v_genEl_E.push_back(genEl_4mom.e());
                treeOutput->v_genEl_px.push_back(genEl_4mom.px());
                treeOutput->v_genEl_py.push_back(genEl_4mom.py());
                treeOutput->v_genEl_pz.push_back(genEl_4mom.pz());
                treeOutput->v_genEl_pT.push_back(genEl_4mom.perp());
                treeOutput->v_genEl_eta.push_back(genEl_4mom.eta());
                treeOutput->v_genEl_phi.push_back(genEl_4mom.phi());
                
                treeOutput->genEl_n++;
            }
        }
    }
    
    
    // Pileup
    edm::Handle <std::vector <PileupSummaryInfo> > pileUps_reco;
    iEvent.getByToken(tok_pileup, pileUps_reco);
    treeOutput->pileup_n = Common::getPileup(pileUps_reco);
    
    
    // Rho
    edm::Handle <double> handle_rho;
    iEvent.getByToken(tok_rho, handle_rho);
    double rho = *handle_rho;
    
    treeOutput->rho = rho;
    
    
    // Gsf electrons from TICL
    edm::Handle <std::vector <reco::GsfElectron> > v_gsfEleFromTICL;
    iEvent.getByToken(tok_gsfEleFromTICL, v_gsfEleFromTICL);
    
    int nEleFromTICL = v_gsfEleFromTICL->size();
    
    std::vector <CLHEP::HepLorentzVector> v_gsfEleFromTICL_4mom;
    
    
    for(int iEle = 0; iEle < nEleFromTICL; iEle++)
    {
        reco::GsfElectron gsfEle = v_gsfEleFromTICL->at(iEle);
        
        CLHEP::HepLorentzVector gsfEleFromTICL_4mom;
        gsfEleFromTICL_4mom.setT(gsfEle.energy());
        gsfEleFromTICL_4mom.setX(gsfEle.px());
        gsfEleFromTICL_4mom.setY(gsfEle.py());
        gsfEleFromTICL_4mom.setZ(gsfEle.pz());
        
        v_gsfEleFromTICL_4mom.push_back(gsfEleFromTICL_4mom);
    }
    
    
    // TICL-ele gen-matching
    TMatrixD mat_gsfEleFromTICL_genEl_deltaR;
    
    std::vector <int> v_gsfEleFromTICL_matchedGenEl_idx;
    
    std::vector <double> v_gsfEleFromTICL_genEl_minDeltaR = Common::getMinDeltaR(
        v_gsfEleFromTICL_4mom,
        v_genEl_4mom,
        mat_gsfEleFromTICL_genEl_deltaR,
        v_gsfEleFromTICL_matchedGenEl_idx
    );
    
    
    for(int iEle = 0; iEle < nEleFromTICL; iEle++)
    {
        reco::GsfElectron gsfEle = v_gsfEleFromTICL->at(iEle);
        CLHEP::HepLorentzVector gsfEleFromTICL_4mom = v_gsfEleFromTICL_4mom.at(iEle);
        
        
        if(gsfEle.pt() < el_minPt || fabs(gsfEle.eta()) < HGCal_minEta || fabs(gsfEle.eta()) > HGCal_maxEta)
        {
            continue;
        }
        
        
        double matchedGenEl_deltaR = v_gsfEleFromTICL_genEl_minDeltaR.at(iEle);
        
        if(matchedGenEl_deltaR > TICLeleGenMatchDR)
        {
            continue;
        }
        
        printf(
            "[%llu] "
            
            "gsfEleFromTICL %d/%d: "
            "E %0.4f, "
            "pT %0.2f, "
            "eta %+0.2f, "
            //"ambiguous %d, "
            
            //"\n"
            //"\t superClus E %0.2f, "
            //"size %d, "
            //"detId %llu, "
            ////"cell %d (x %+0.2f, y %+0.2f), "
            //"type %d, "
            //"neighbors %d, %d, "
            ////"sector %d, "
            //"dist %0.2e, "
            //
            //"\n"
            //"\t\t superClus seed: E %0.2f, eta %+0.2f, z %+0.2f, size %d"
            //
            //"\n"
            //"\t gsfTrack p %0.2f, "
            //
            //"\n"
            //"\t trackMomentumAtVtx p %0.2f, "// (%0.2f), "
            //"pT %0.2f, "// (%0.2f), "
            
            "\n",
            
            eventNumber,
            
            iEle+1, nEleFromTICL,
            gsfEle.energy(),
            gsfEle.pt(),
            gsfEle.eta()
            //gsfEle.ambiguous(),
            
            //gsfEle.superCluster()->energy(),
            //(int) gsfEle.superCluster()->size(),
            //(long long) centroid_detId.rawId(),
            ////centroid_HGCEEDetId.cell(), centroidCell_pos.x(), centroidCell_pos.y(),
            //topo_HGCalEE.decode(centroid_detId).iType,
            //(int) v_neighbour7_detId.size(), (int) v_neighbour19_detId.size(),
            ////centroid_HGCEEDetId.sector(),
            //dist_min,
            //
            //gsfEle.superCluster()->seed().get()->energy(),
            //gsfEle.superCluster()->seed().get()->eta(),
            //gsfEle.superCluster()->seed().get()->z(),
            //(int) gsfEle.superCluster()->seed().get()->size(),
            //
            //gsfEle.gsfTrack()->p(),
            //gsfEle.trackMomentumAtVtx().r(),// std::sqrt(gsfEle.trackMomentumAtVtx().mag2()),
            //gsfEle.trackMomentumAtVtx().rho()//, std::sqrt(gsfEle.trackMomentumAtVtx().perp2())
        );
        
        int matchedGenEl_idx = v_gsfEleFromTICL_matchedGenEl_idx.at(iEle);
        
        treeOutput->v_gsfEleFromTICL_genEl_minDeltaR.push_back(matchedGenEl_deltaR);
        treeOutput->v_gsfEleFromTICL_nearestGenEl_idx.push_back(matchedGenEl_idx);
        
        double matchedGenEl_energy = -99;
        double matchedGenEl_pT = -99;
        double matchedGenEl_eta = -99;
        double matchedGenEl_phi = -99;
        
        if(matchedGenEl_idx >= 0)
        {
            matchedGenEl_energy = v_genEl_4mom.at(matchedGenEl_idx).e();
            matchedGenEl_pT = v_genEl_4mom.at(matchedGenEl_idx).perp();
            matchedGenEl_eta = v_genEl_4mom.at(matchedGenEl_idx).eta();
            matchedGenEl_phi = v_genEl_4mom.at(matchedGenEl_idx).phi();
        }
        
        treeOutput->v_gsfEleFromTICL_matchedGenEl_E.push_back(matchedGenEl_energy);
        treeOutput->v_gsfEleFromTICL_matchedGenEl_pT.push_back(matchedGenEl_pT);
        treeOutput->v_gsfEleFromTICL_matchedGenEl_eta.push_back(matchedGenEl_eta);
        treeOutput->v_gsfEleFromTICL_matchedGenEl_phi.push_back(matchedGenEl_phi);
        
        
        treeOutput->v_gsfEleFromTICL_E.push_back(gsfEle.energy());
        treeOutput->v_gsfEleFromTICL_px.push_back(gsfEle.px());
        treeOutput->v_gsfEleFromTICL_py.push_back(gsfEle.py());
        treeOutput->v_gsfEleFromTICL_pz.push_back(gsfEle.pz());
        
        treeOutput->v_gsfEleFromTICL_pT.push_back(gsfEle.pt());
        treeOutput->v_gsfEleFromTICL_eta.push_back(gsfEle.eta());
        treeOutput->v_gsfEleFromTICL_phi.push_back(gsfEle.phi());
        
        treeOutput->v_gsfEleFromTICL_ET.push_back(gsfEle.et());
        
        treeOutput->gsfEleFromTICL_n++;
        
        
        //std::vector <DetId> v_SC_seedId = gsfEle.superCluster()->getSeedIds();
        //
        //printf("SC_nCluster %d, SC_nSeed %d \n", (int) gsfEle.superCluster()->clusters().size(), (int) v_SC_seedId.size());
        //printf("SC_seedIds: ");
        //
        //for(int iSeed = 0; iSeed < (int) v_SC_seedId.size(); iSeed++)
        //{
        //    printf("[%u] ", v_SC_seedId.at(iSeed).rawId());
        //}
        //
        //printf("\n");
        
        edm::PtrVector <reco::CaloCluster> v_superClus_clus = gsfEle.superCluster()->clusters();
        
        printf("SC_clusterIds: ");
        
        for(int iCluster = 0; iCluster < (int) v_superClus_clus.size(); iCluster++)
        {
            const reco::CaloCluster *cluster = v_superClus_clus[iCluster].get();
            
            printf("[%u] ", cluster->seed().rawId());
        }
        
        printf("\n");
    }
    
    // Fill tree
    treeOutput->fill();
    
    //#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    //ESHandle<SetupData> pSetup;
    //iSetup.get<SetupRecord>().get(pSetup);
    //#endif
    
    printf("\n\n");
    
    fflush(stdout);
    fflush(stderr);
}


double TreeMaker::getDeltaPhi(double phi1, double phi2)
{
    double deltaPhi = phi1 - phi2;
    
    deltaPhi = (deltaPhi > +M_PI)? (deltaPhi - 2*M_PI): deltaPhi;
    deltaPhi = (deltaPhi < -M_PI)? (deltaPhi + 2*M_PI): deltaPhi;
    
    //deltaPhi = (deltaPhi > +M_PI)? (2*M_PI - deltaPhi): deltaPhi;
    //deltaPhi = (deltaPhi < -M_PI)? (2*M_PI + deltaPhi): deltaPhi;
    
    return deltaPhi;
}


// ------------ method called once each job just before starting event loop  ------------
void
TreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TreeMaker::endJob()
{
    
    printf("minLayer = %d, maxLayer = %d \n", minLayer, maxLayer);
    
    
    fflush(stdout);
    fflush(stderr);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeMaker);
