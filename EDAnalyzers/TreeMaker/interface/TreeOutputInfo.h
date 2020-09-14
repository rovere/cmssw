# ifndef TreeOutputInfo_H
# define TreeOutputInfo_H


# include <iostream>
# include <map>
# include <stdlib.h>
# include <string>
# include <type_traits>
# include <utility>
# include <vector>

# include <TH1F.h>
# include <TH2F.h>
# include <TMatrixD.h>
# include <TROOT.h>
# include <TTree.h> 
# include <TVectorD.h> 

# include "EDAnalyzers/TreeMaker/interface/Constants.h"


namespace TreeOutputInfo
{
    class TreeOutput
    {
        public :
        
        
        TTree *tree;
        
        
        // Run info //
        ULong64_t runNumber;
        ULong64_t eventNumber;
        ULong64_t luminosityNumber;
        ULong64_t bunchCrossingNumber;
        
        
        // Gen electron //
        int genEl_n;
        std::vector <double> v_genEl_E;
        std::vector <double> v_genEl_px;
        std::vector <double> v_genEl_py;
        std::vector <double> v_genEl_pz;
        std::vector <double> v_genEl_pT;
        std::vector <double> v_genEl_eta;
        std::vector <double> v_genEl_phi;
        
        
        // Pileup //
        int pileup_n;
        
        
        // Rho //
        double rho;
        
        
        // Tracksters //
        int trackster_n;
        std::vector <double> v_trackster_E;
        std::vector <double> v_trackster_x;
        std::vector <double> v_trackster_y;
        std::vector <double> v_trackster_z;
        std::vector <double> v_trackster_eta;
        std::vector <double> v_trackster_phi;
        std::vector <double> v_trackster_ET;
        
        
        double gsfEleFromTICL_n;
        std::vector <double> v_gsfEleFromTICL_E;
        std::vector <double> v_gsfEleFromTICL_px;
        std::vector <double> v_gsfEleFromTICL_py;
        std::vector <double> v_gsfEleFromTICL_pz;
        std::vector <double> v_gsfEleFromTICL_pT;
        std::vector <double> v_gsfEleFromTICL_eta;
        std::vector <double> v_gsfEleFromTICL_phi;
        std::vector <double> v_gsfEleFromTICL_ET;
        
        std::vector <double> v_gsfEleFromTICL_genEl_minDeltaR;
        std::vector <double> v_gsfEleFromTICL_nearestGenEl_idx;
        std::vector <double> v_gsfEleFromTICL_matchedGenEl_E;
        std::vector <double> v_gsfEleFromTICL_matchedGenEl_pT;
        std::vector <double> v_gsfEleFromTICL_matchedGenEl_eta;
        std::vector <double> v_gsfEleFromTICL_matchedGenEl_phi;
        
        
        char name[500];
        
        
        TreeOutput(std::string details, edm::Service<TFileService> fs)
        {
            printf("Loading custom ROOT dictionaries. \n");
            gROOT->ProcessLine(".L EDAnalyzers/TreeMaker/interface/CustomRootDict.cc+");
            printf("Loaded custom ROOT dictionaries. \n");
            
            tree = fs->make<TTree>(details.c_str(), details.c_str());
            
            
            // Run info //
            tree->Branch("runNumber", &runNumber);
            tree->Branch("eventNumber", &eventNumber);
            tree->Branch("luminosityNumber", &luminosityNumber);
            tree->Branch("bunchCrossingNumber", &bunchCrossingNumber);
            
            
            // Gen electron //
            sprintf(name, "genEl_n");
            tree->Branch(name, &genEl_n);
            
            sprintf(name, "genEl_E");
            tree->Branch(name, &v_genEl_E);
            
            sprintf(name, "genEl_px");
            tree->Branch(name, &v_genEl_px);
            
            sprintf(name, "genEl_py");
            tree->Branch(name, &v_genEl_py);
            
            sprintf(name, "genEl_pz");
            tree->Branch(name, &v_genEl_pz);
            
            sprintf(name, "genEl_pT");
            tree->Branch(name, &v_genEl_pT);
            
            sprintf(name, "genEl_eta");
            tree->Branch(name, &v_genEl_eta);
            
            sprintf(name, "genEl_phi");
            tree->Branch(name, &v_genEl_phi);

            
            // Pileup //
            sprintf(name, "pileup_n");
            tree->Branch(name, &pileup_n);
            
            
            // Rho //
            sprintf(name, "rho");
            tree->Branch(name, &rho);
            
            
            //
            sprintf(name, "gsfEleFromTICL_n");
            tree->Branch(name, &gsfEleFromTICL_n);
            
            sprintf(name, "gsfEleFromTICL_E");
            tree->Branch(name, &v_gsfEleFromTICL_E);
            
            sprintf(name, "gsfEleFromTICL_px");
            tree->Branch(name, &v_gsfEleFromTICL_px);
            
            sprintf(name, "gsfEleFromTICL_py");
            tree->Branch(name, &v_gsfEleFromTICL_py);
            
            sprintf(name, "gsfEleFromTICL_pz");
            tree->Branch(name, &v_gsfEleFromTICL_pz);
            
            sprintf(name, "gsfEleFromTICL_pT");
            tree->Branch(name, &v_gsfEleFromTICL_pT);
            
            sprintf(name, "gsfEleFromTICL_eta");
            tree->Branch(name, &v_gsfEleFromTICL_eta);
            
            sprintf(name, "gsfEleFromTICL_phi");
            tree->Branch(name, &v_gsfEleFromTICL_phi);
            
            sprintf(name, "gsfEleFromTICL_ET");
            tree->Branch(name, &v_gsfEleFromTICL_ET);
            
            sprintf(name, "gsfEleFromTICL_genEl_minDeltaR");
            tree->Branch(name, &v_gsfEleFromTICL_genEl_minDeltaR);
            
            sprintf(name, "gsfEleFromTICL_nearestGenEl_idx");
            tree->Branch(name, &v_gsfEleFromTICL_nearestGenEl_idx);
            
            sprintf(name, "gsfEleFromTICL_matchedGenEl_E");
            tree->Branch(name, &v_gsfEleFromTICL_matchedGenEl_E);
            
            sprintf(name, "gsfEleFromTICL_matchedGenEl_pT");
            tree->Branch(name, &v_gsfEleFromTICL_matchedGenEl_pT);
            
            sprintf(name, "gsfEleFromTICL_matchedGenEl_eta");
            tree->Branch(name, &v_gsfEleFromTICL_matchedGenEl_eta);
            
            sprintf(name, "gsfEleFromTICL_matchedGenEl_phi");
            tree->Branch(name, &v_gsfEleFromTICL_matchedGenEl_phi);
        }
        
        
        void fill()
        {
            tree->Fill();
        }
        
        
        void clear()
        {
            // Gen electron //
            genEl_n = 0;
            v_genEl_E.clear();
            v_genEl_px.clear();
            v_genEl_py.clear();
            v_genEl_pz.clear();
            v_genEl_pT.clear();
            v_genEl_eta.clear();
            v_genEl_phi.clear();
            
            
            // Pileup //
            pileup_n = 0;
            
            
            // Rho //
            rho = 0;
            
            
            //
            gsfEleFromTICL_n = 0;
            v_gsfEleFromTICL_E.clear();
            v_gsfEleFromTICL_px.clear();
            v_gsfEleFromTICL_py.clear();
            v_gsfEleFromTICL_pz.clear();
            v_gsfEleFromTICL_pT.clear();
            v_gsfEleFromTICL_eta.clear();
            v_gsfEleFromTICL_phi.clear();
            v_gsfEleFromTICL_ET.clear();
            
            v_gsfEleFromTICL_genEl_minDeltaR.clear();
            v_gsfEleFromTICL_nearestGenEl_idx.clear();
            v_gsfEleFromTICL_matchedGenEl_E.clear();
            v_gsfEleFromTICL_matchedGenEl_pT.clear();
            v_gsfEleFromTICL_matchedGenEl_eta.clear();
            v_gsfEleFromTICL_matchedGenEl_phi.clear();
        }
    };
}


# endif
