//Runs Pythia events in some pT hard bin, clusters the final state particles with
//FastJet's anti-kT, applies cuts and mock tracking efficiency, fills trees, adds
//back randomly sampled particles in attempt to reproduce particle-level jets
//from detector-level.
//Isaac Mooney - WSU, 2017

#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "TROOT.h"
#include "TRandom.h"
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

#include "analysis_funcs.h"

using namespace fastjet;
using namespace Pythia8;
using namespace std;
using namespace analysis;

//Passing in any arguments will run the program in pedantic mode, filling all containers
//with all information.

int main(int argc, char ** argv) {
    // Generator. Process selection. STAR initialization. Trees.
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    //Initialize Pythia
    Pythia pythia;
    analysis::pythia_init(pythia);
    
    //-------------------------File(s) and containers-------------------------//
    
    //MIGHT WANT TO JUST USE fout INSTEAD OF fout + fit2d???
    //TFile *fit2d = new TFile("out/fit2D05.root", "RECREATE");
    
    
    TFile *fout = new TFile("out/trackeffic510.root","RECREATE");  // Create output file
    analysis::containers * holder = new containers();             // Create trees and pseudojet vectors
    holder->SetBranches();                                        // Create branches
    
    //------------------------------------------------------------------------//
    
    //Define jet for FastJet clustering
    JetDefinition jet_def(antikt_algorithm, R);
    
    //Histograms:
    TH1 * c_cons_pt = new TH1D("cons_pt","constituent pt", 200,0,20);
    TH1 * ptbefore = new TH1D("ptbefore", "", 100,0,30);
    TH1 * ptprecorrection = new TH1D("ptprecorrection", "", 100,0,30);
    TH1 * numbefore = new TH1D("numbefore","",30,0,30);
    TH1 * numprecorrection = new TH1D("numprecorrection","",30,0,30);
    //creates the 2D histogram which will be fit for a given slice and used for sampling number loss
    TH2 * fit = new TH2D("fit","",20,0,20,10,0,10);
    //creates a better (semi-)1D fit than the naive 1D number loss fit
    TH2 * numlostvdetpt = new TH2D("numlostvdetpt","num lost v. detector pt",100,0,30,20,0,20);
    
    for (unsigned iEvent = 0; iEvent < NumEvents; ++iEvent) {
        holder->Clear();
        
        if (!pythia.next()) continue;
        
        //add_particles fills containers with particles generated from the event, applies cuts,
        //mocks up tracking efficiency.
        analysis::add_particles(pythia.event, holder, argc);
        
        //~~~~~~~~~~~~~~~~~~~~~~~~CLUSTER AND GET JETS~~~~~~~~~~~~~~~~~~~~~~~~//
        
        //pedantic mode, cluster all jets
        ClusterSequence uncut(holder->uncut_part,       jet_def);
        ClusterSequence cut(holder->cut_part,           jet_def);
        ClusterSequence cut2(holder->cut2_part,         jet_def);
        ClusterSequence effic(holder->effic_part,       jet_def);
        vector<PseudoJet> uncut_jets    = sorted_by_pt(uncut.inclusive_jets());
        vector<PseudoJet> cut_jets      = sorted_by_pt(cut.inclusive_jets());
        vector<PseudoJet> cut2_jets     = sorted_by_pt(cut2.inclusive_jets());
        vector<PseudoJet> effic_jets    = sorted_by_pt(effic.inclusive_jets());
        
        ClusterSequence c_uncut(holder->c_uncut_part,   jet_def);
        ClusterSequence c_cut(holder->c_cut_part,       jet_def);
        ClusterSequence c_cut2(holder->c_cut2_part,     jet_def);
        ClusterSequence c_effic(holder->c_effic_part,   jet_def);
        vector<PseudoJet> c_uncut_jets  = sorted_by_pt(c_uncut.inclusive_jets());
        vector<PseudoJet> c_cut_jets    = sorted_by_pt(c_cut.inclusive_jets());
        vector<PseudoJet> c_cut2_jets   = sorted_by_pt(c_cut2.inclusive_jets());
        vector<PseudoJet> c_effic_jets  = sorted_by_pt(c_effic.inclusive_jets());
        
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        
        //cout << "Clustering with " << jet_def.description() << '\n';
        
        if (argc == 2) { //this is pedantic mode, fill all containers
            //~~~~~~~~~~~~~~~~~~~~FILL TREES WITH LEAD PTS~~~~~~~~~~~~~~~~~~~~//
            
            //leading jet pt spectra (total and only charged)
            //(a) without cuts, (b) with particle cuts, (c) with jet cuts, and (d) with cuts + efficiency
            
            //(a):
            if (uncut_jets.size() != 0) {
                holder->uncut_leadPt = uncut_jets[0].pt();
            }
            
            if (c_uncut_jets.size() != 0) {
                holder->c_uncut_leadPt = c_uncut_jets[0].pt();
            }
            
            //(b):
            if (cut_jets.size() != 0) {
                holder->cut_leadPt = cut_jets[0].pt();
            }
            
            if (c_cut_jets.size() != 0) {
                holder->c_cut_leadPt = c_cut_jets[0].pt();
            }
            
            //(c):
            if(cut2_jets.size() != 0) {
                if (jet_cuts(cut2_jets[0])) {
                    //leading jet pt spectrum
                    holder->cut2_leadPt = cut2_jets[0].pt();
                }
            }
            
            if (c_cut2_jets.size() != 0 && jet_cuts(c_cut2_jets[0])) {
                holder->c_cut2_leadPt = c_cut2_jets[0].pt();
            }
            
            //(d):
            if (effic_jets.size() != 0) {
                if (jet_cuts(effic_jets[0])) {
                    holder->effic_leadPt = effic_jets[0].pt();
                }
            }
            
            if (c_effic_jets.size() != 0 && jet_cuts(c_effic_jets[0])) {
                holder->c_effic_leadPt = c_effic_jets[0].pt();
            }
            
            //branches in tree "tlead" have been populated, so fill the tree with these branches
            holder->tlead->Fill();
            
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            
            //~~~~~~~~~~~~~~FILL TREES WITH JETS & CONSTITUENTS~~~~~~~~~~~~~~~//
            
            //pseudojets without cuts
            for (unsigned i = 0; i < uncut_jets.size(); ++ i) {
                holder->uncut_px    = uncut_jets[i].px();
                holder->uncut_py    = uncut_jets[i].py();
                holder->uncut_pz    = uncut_jets[i].pz();
                holder->uncut_E     = uncut_jets[i].e();
                holder->uncut_phi   = uncut_jets[i].phi();
                holder->uncut_eta   = uncut_jets[i].eta();
                holder->uncut_Pt    = uncut_jets[i].pt();
                
                holder->tuncut->Fill();
                vector<PseudoJet> uncut_cons = uncut_jets[i].constituents();
                for (unsigned j = 0; j < uncut_cons.size(); ++ j) {
                    holder->uncut_cons_px   = uncut_cons[j].px();
                    holder->uncut_cons_py   = uncut_cons[j].py();
                    holder->uncut_cons_pz   = uncut_cons[j].pz();
                    holder->uncut_cons_E    = uncut_cons[j].e();
                    holder->uncut_cons_phi  = uncut_cons[j].phi();
                    holder->uncut_cons_eta  = uncut_cons[j].eta();
                    holder->uncut_cons_Pt   = uncut_cons[j].pt();
                    
                    holder->tuncutcons->Fill();
                }
            }
            
            //pseudojets without cuts (charged)
            for (unsigned i = 0; i < c_uncut_jets.size(); ++ i) {
                holder->c_uncut_px  = c_uncut_jets[i].px();
                holder->c_uncut_py  = c_uncut_jets[i].py();
                holder->c_uncut_pz  = c_uncut_jets[i].pz();
                holder->c_uncut_E   = c_uncut_jets[i].e();
                holder->c_uncut_phi = c_uncut_jets[i].phi();
                holder->c_uncut_eta = c_uncut_jets[i].eta();
                holder->c_uncut_Pt  = c_uncut_jets[i].pt();
                
                holder->tcuncut->Fill();
                
                vector<PseudoJet> c_uncut_cons = c_uncut_jets[i].constituents();
                for (unsigned j = 0; j < c_uncut_cons.size(); ++ j) {
                    holder->c_uncut_cons_px     = c_uncut_cons[j].px();
                    holder->c_uncut_cons_py     = c_uncut_cons[j].py();
                    holder->c_uncut_cons_pz     = c_uncut_cons[j].pz();
                    holder->c_uncut_cons_E      = c_uncut_cons[j].e();
                    holder->c_uncut_cons_phi    = c_uncut_cons[j].phi();
                    holder->c_uncut_cons_eta    = c_uncut_cons[j].eta();
                    holder->c_uncut_cons_Pt     = c_uncut_cons[j].pt();
                    
                    holder->tcuncutcons->Fill();
                }
            }
            
            //pseudojets with particle level cuts
            for (unsigned i = 0; i < cut_jets.size(); ++ i) {
                holder->cut_px  = cut_jets[i].px();
                holder->cut_py  = cut_jets[i].py();
                holder->cut_pz  = cut_jets[i].pz();
                holder->cut_E   = cut_jets[i].e();
                holder->cut_phi = cut_jets[i].phi();
                holder->cut_eta = cut_jets[i].eta();
                holder->cut_Pt  = cut_jets[i].pt();
                
                holder->tcut->Fill();
            }
            
            for (unsigned i = 0; i < c_cut_jets.size(); ++ i) {
                holder->c_cut_px    = c_cut_jets[i].px();
                holder->c_cut_py    = c_cut_jets[i].py();
                holder->c_cut_pz    = c_cut_jets[i].pz();
                holder->c_cut_E     = c_cut_jets[i].e();
                holder->c_cut_phi   = c_cut_jets[i].phi();
                holder->c_cut_eta   = c_cut_jets[i].eta();
                holder->c_cut_Pt    = c_cut_jets[i].pt();
                
                holder->tccut->Fill();
            }
            
            //pseudojets with both sets of cuts
            for (unsigned i = 0; i < cut2_jets.size(); ++ i) {
                //pseudojet axis cut
                if (jet_cuts(cut2_jets[i])) {
                    holder->cut2_px     = cut2_jets[i].px();
                    holder->cut2_py     = cut2_jets[i].py();
                    holder->cut2_pz     = cut2_jets[i].pz();
                    holder->cut2_E      = cut2_jets[i].e();
                    holder->cut2_phi    = cut2_jets[i].phi();
                    holder->cut2_eta    = cut2_jets[i].eta();
                    holder->cut2_Pt     = cut2_jets[i].pt();
                    
                    holder->tcut2->Fill();
                    
                    
                    vector<PseudoJet> cut2_cons = cut2_jets[i].constituents();
                    for (unsigned j = 0; j < cut2_cons.size(); ++ j) {
                        holder->cut2_cons_px     = cut2_cons[j].px();
                        holder->cut2_cons_py     = cut2_cons[j].py();
                        holder->cut2_cons_pz     = cut2_cons[j].pz();
                        holder->cut2_cons_E      = cut2_cons[j].e();
                        holder->cut2_cons_phi    = cut2_cons[j].phi();
                        holder->cut2_cons_eta    = cut2_cons[j].eta();
                        holder->cut2_cons_Pt     = cut2_cons[j].pt();
                        
                        holder->tcut2cons->Fill();
                    }
                }
            }
            
            for (unsigned i = 0; i < c_cut2_jets.size(); ++ i) {
                if (jet_cuts(c_cut2_jets[i])) {
                    holder->c_cut2_px   = c_cut2_jets[i].px();
                    holder->c_cut2_py   = c_cut2_jets[i].py();
                    holder->c_cut2_pz   = c_cut2_jets[i].pz();
                    holder->c_cut2_E    = c_cut2_jets[i].e();
                    holder->c_cut2_phi  = c_cut2_jets[i].phi();
                    holder->c_cut2_eta  = c_cut2_jets[i].eta();
                    holder->c_cut2_Pt   = c_cut2_jets[i].pt();
                    
                    holder->tccut2->Fill();
                    
                    vector<PseudoJet> c_cut2_cons = c_cut2_jets[i].constituents();
                    for (unsigned j = 0; j < c_cut2_cons.size(); ++ j) {
                        holder->c_cut2_cons_px     = c_cut2_cons[j].px();
                        holder->c_cut2_cons_py     = c_cut2_cons[j].py();
                        holder->c_cut2_cons_pz     = c_cut2_cons[j].pz();
                        holder->c_cut2_cons_E      = c_cut2_cons[j].e();
                        holder->c_cut2_cons_phi    = c_cut2_cons[j].phi();
                        holder->c_cut2_cons_eta    = c_cut2_cons[j].eta();
                        holder->c_cut2_cons_Pt     = c_cut2_cons[j].pt();
                    }
                    holder->tccut2cons->Fill();
                }
            }
            
            //pseudojets with cuts + charged particle tracking efficiency
            for (unsigned i = 0; i < effic_jets.size(); ++ i) {
                if (jet_cuts(effic_jets[i])) {
                    holder->effic_px    = effic_jets[i].px();
                    holder->effic_py    = effic_jets[i].py();
                    holder->effic_pz    = effic_jets[i].pz();
                    holder->effic_E     = effic_jets[i].e();
                    holder->effic_phi   = effic_jets[i].phi();
                    holder->effic_eta   = effic_jets[i].eta();
                    holder->effic_Pt    = effic_jets[i].pt();
                    
                    holder->teffic->Fill();
                    
                    vector<PseudoJet> effic_cons = effic_jets[i].constituents();
                    for (unsigned j = 0; j < effic_cons.size(); ++ j) {
                        holder->effic_cons_px     = effic_cons[j].px();
                        holder->effic_cons_py     = effic_cons[j].py();
                        holder->effic_cons_pz     = effic_cons[j].pz();
                        holder->effic_cons_E      = effic_cons[j].e();
                        holder->effic_cons_phi    = effic_cons[j].phi();
                        holder->effic_cons_eta    = effic_cons[j].eta();
                        holder->effic_cons_Pt     = effic_cons[j].pt();
                        
                        holder->tefficcons->Fill();
                    }
                }
            }
            
            for (unsigned i = 0; i < c_effic_jets.size(); ++ i) {
                if (jet_cuts(c_effic_jets[i])) {
                    holder->c_effic_px  = c_effic_jets[i].px();
                    holder->c_effic_py  = c_effic_jets[i].py();
                    holder->c_effic_pz  = c_effic_jets[i].pz();
                    holder->c_effic_E   = c_effic_jets[i].e();
                    holder->c_effic_phi = c_effic_jets[i].phi();
                    holder->c_effic_eta = c_effic_jets[i].eta();
                    holder->c_effic_Pt  = c_effic_jets[i].pt();
                    
                    holder->tceffic->Fill();
                
                    vector<PseudoJet> c_effic_cons = c_effic_jets[i].constituents();
                    for (unsigned j = 0; j < c_effic_cons.size(); ++ j) {
                        holder->c_effic_cons_px     = c_effic_cons[j].px();
                        holder->c_effic_cons_py     = c_effic_cons[j].py();
                        holder->c_effic_cons_pz     = c_effic_cons[j].pz();
                        holder->c_effic_cons_E      = c_effic_cons[j].e();
                        holder->c_effic_cons_phi    = c_effic_cons[j].phi();
                        holder->c_effic_cons_eta    = c_effic_cons[j].eta();
                        holder->c_effic_cons_Pt     = c_effic_cons[j].pt();
                        
                        holder->tcefficcons->Fill();
                    }
                }
            }
            
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            
        }//Exiting pedantic conditional. Rest will fill no matter what.
        
        //Filling pre-adding-back histograms with charged leading jet info
        if (c_cut2_jets.size() != 0) {
            if (jet_cuts(c_cut2_jets[0])) {
                ptbefore->Fill(c_cut2_jets[0].pt());
                numbefore->Fill(c_cut2_jets[0].constituents().size());
                //Now fill charged leading jet constituent pt spectrum histogram
                vector<PseudoJet> c_cut2_cons = c_cut2_jets[0].constituents();
                for (unsigned j = 0; j < c_cut2_cons.size(); ++ j) {
                    c_cons_pt->Fill(c_cut2_cons[j].pt());
                }
            }
        }
        
        //Filling 'fit' histograms and other pre-correction histograms
        if (c_effic_jets.size() != 0 && c_cut2_jets.size() != 0) {
            if (jet_cuts(c_effic_jets[0]) && jet_cuts(c_cut2_jets[0])) {
                //leading jet number loss vs. detector-level pT filled here
                numlostvdetpt->Fill(c_effic_jets[0].pt(), c_cut2_jets[0].constituents().size() - c_effic_jets[0].constituents().size());
                //filling 2D fit
                fit->Fill(c_effic_jets[0].constituents().size(), c_cut2_jets[0].constituents().size() - c_effic_jets[0].constituents().size());
            }
            if (jet_cuts(c_effic_jets[0])) {
                //filling pre-adding-back constituent multiplicity and pt for leading jet
                numprecorrection->Fill(c_effic_jets[0].constituents().size());
                ptprecorrection->Fill(c_effic_jets[0].pt());
            }
        }
        
        //Checks pt & num difference between leading jet with efficiency correction and
        //geometrically closest jet without efficiency correction.
        //if pedantic mode on, does this for both charged, and neutral + charged jets
        if (argc == 2) {
            analysis::geometric_diff(effic_jets, cut2_jets, holder->ptdiff, holder->num_diff, holder->num_before, holder->num_after, holder->rel_diff);
        }
        analysis::geometric_diff(c_effic_jets, c_cut2_jets, holder->c_ptdiff, holder->c_num_diff, holder->c_num_before, holder->c_num_after, holder->c_rel_diff);
        
        holder->tdiffs->Fill();
        
    } //End of event loop!
    
    pythia.stat();
    
    cout << "Writing to " << fout->GetName() << '\n';
    
    //tout->Print(); if status of a particular tree is desired.
    //save trees to output root file
    holder->write();
    fout->Write();
    
    //fout->Close();
    
    return 0;
}
