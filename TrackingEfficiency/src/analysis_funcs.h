//  analysis_funcs.h
//  Created by Isaac Mooney on 5/4/17.

#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"

#include <iostream>
#include <sstream>
#include <random>
#include <chrono>

#ifndef analysis_funcs_h
#define analysis_funcs_h

namespace analysis {
    class containers;
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        const int NumEvents     = 100000;
    //CUTS:
        //Particle level:
        const double eta_cut    = 1.0;
        const double pt_cut     = 0.2;
        //Jet level:
        const double R          = 0.4;
        const double jet_eta    = 1 - R;
    
    //EFFICIENCY:
        const double efficiency = 0.8;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //seeds random number generator
    std::mt19937 rands();
    //initializes Pythia
    void pythia_init(Pythia8::Pythia & pythia);
    //checks that particles pass particle-level cuts
    bool constituent_cuts(const Pythia8::Particle particle);
    //check that pseudojets pass jet-level cuts
    bool jet_cuts(const fastjet::PseudoJet pseudojet);
    //tracking efficiency selection
    bool efficiency_cut(const double effic_num, const Pythia8::Particle particle);
    //fills particle vectors after appropriate cuts are passed
    void add_particles(std::uniform_real_distribution<> dis, std::mt19937 g, Pythia8::Event event, containers * container);
    //possibly implementing later: clusters the particles into jets
    //void cluster(containers *container, const fastjet::JetDefinition jet_def);
    //a test to see that pT spectrum is steeply falling if 1st and 2nd leading jets are removed
    void without_lead2_jets(const std::vector<fastjet::PseudoJet> effic_jets, double & without_leadsublead, TTree* ttests);
    //examines the pT difference between efficiency-corrected, and non-efficiency-corrected 'leading' jets
    void geometric_pt_diff(const std::vector<fastjet::PseudoJet> effic_jets, const std::vector<fastjet::PseudoJet> cut2_jets, double & ptdiff, int & num_diff);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~CONTAINER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    class containers {
    private:
        //MOVE THINGS HERE LATER AND WRITE GETTERS? SEEMS LIKE A PAIN
    public:
        //trees
        TTree *tuncut, *tcuncut, *tcut, *tccut, *tcut2, *tccut2, *teffic, *tceffic, *tlead, *tdiffs, *ttests, *tuncutcons, *tcuncutcons;
        
        //branches
        double uncut_px, uncut_py, uncut_pz, uncut_E, uncut_phi, uncut_eta, uncut_leadPt, uncut_Pt;
        double c_uncut_px, c_uncut_py, c_uncut_pz, c_uncut_E, c_uncut_phi, c_uncut_eta, c_uncut_leadPt, c_uncut_Pt;
        double uncut_cons_px, uncut_cons_py, uncut_cons_pz, uncut_cons_E, uncut_cons_phi, uncut_cons_eta, uncut_cons_Pt;
        double c_uncut_cons_px, c_uncut_cons_py, c_uncut_cons_pz, c_uncut_cons_E, c_uncut_cons_phi, c_uncut_cons_eta, c_uncut_cons_Pt;
        double cut_px, cut_py, cut_pz, cut_E, cut_phi, cut_eta, cut_leadPt, cut_Pt;
        double c_cut_px, c_cut_py, c_cut_pz, c_cut_E, c_cut_phi, c_cut_eta, c_cut_leadPt, c_cut_Pt;
        double cut2_px, cut2_py, cut2_pz, cut2_E, cut2_phi, cut2_eta, cut2_leadPt, cut2_Pt;
        double c_cut2_px, c_cut2_py, c_cut2_pz, c_cut2_E, c_cut2_phi, c_cut2_eta, c_cut2_leadPt, c_cut2_Pt;
        double effic_px, effic_py, effic_pz, effic_E, effic_phi, effic_eta, effic_leadPt, effic_Pt;
        double without_leadsublead;
        double c_effic_px, c_effic_py, c_effic_pz, c_effic_E, c_effic_phi, c_effic_eta, c_effic_leadPt, c_effic_Pt;
        double ptdiff, c_ptdiff;
        int num_diff, c_num_diff;
        
        //particles
        std::vector<fastjet::PseudoJet> uncut_part, c_uncut_part, cut_part, cut2_part, c_cut_part, c_cut2_part, effic_part, c_effic_part;
        //LATER EXTEND TO INCLUDE JET VECTORS IN ADDITION TO THE ABOVE PARTICLE VECTORS? TRICKY THOUGH
        
        //jets
        //std::vector<fastjet::PseudoJet> uncut_jets, c_uncut_jets, cut_jets, cut2_jets, c_cut_jets, c_cut2_jets, effic_jets, c_effic_jets;
        
        //cluster sequences
        //fastjet::ClusterSequence uncut(), cut(), cut2(), effic(), c_uncut(), c_cut(), c_cut2(), c_effic();
        
        
        containers();               //constructor (initializes the trees, reserves space for vectors)
        virtual ~containers() {}    //destructor (shouldn't need to implement?)
        
        
        //getters & setters
        
        void SetBranches();         //setter
        //SetJets();
        
        void Clear();               //clears the particle vectors for each event
        void write();               //writes the trees to file
    };
    
}

#endif /* analysis_funcs_h */
