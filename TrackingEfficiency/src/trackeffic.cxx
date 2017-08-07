//Runs events in a narrow pT hard bin, clusters the final state
//particles with anti-kT, applies some non-stringent cuts,
//and fills trees with & without a tracking efficiency selection cut.
//Isaac Mooney - WSU, March 2017, updated May 2017

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

int main() {
    // Generator. Process selection. STAR initialization. Trees.
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    Pythia pythia;
    analysis::pythia_init(pythia);
    
    
    TFile *fit2d = TFile::Open("out/fit2d2232.root");
    TH2* fit = (TH2D*)fit2d->Get("fit");
    /*
    
     TFile *fit2d = new TFile("out/fit2D2232.root", "RECREATE");
     TH2 * fit = new TH2D("fit","",20,0,20,10,0,10);
    */
    
    TFile *fout     = new TFile("out/trackeffic2232.root","RECREATE");  // Create output file
    analysis::containers * test = new containers();                 // Create trees and pseudojet vectors
    test->SetBranches();                                            // Create branches
    
    
    JetDefinition jet_def(antikt_algorithm, R);
    
    // Begin event loop. Generate event. Skip if error. List first one.
    
    int count = 0, c_count = 0;
    
    //TESTS!!!
    TF1 * nb = new TF1("nb","6.44464e-01*ROOT::Math::negative_binomial_pdf(x,7.87891e-01,3.75129e+00)", 0, 30);
    TF1 * nbnewtest = new TF1("nbnewtest","9.99801e-01*ROOT::Math::negative_binomial_pdf(x,6.33335e-01,1.49188e+00)",0,30);
    TH1 * pid = new TH1I("pid", "pid", 10000, 0, 10000);
    TH1 * c_cons_pt = new TH1D("cons_pt","constituent pt", 200,0,20);
    TH1 * ptlost = new TH1D("ptlost","pt loss",300,-30,1);
    TH1 * ptafter = new TH1D("ptafter", "", 100,0,30);
    TH1 * ptbefore = new TH1D("ptbefore", "", 100,0,30);
    TH1 * ptprecorrection = new TH1D("ptprecorrection", "", 100,0,30);
    TH1 * numbefore = new TH1D("numbefore","",30,0,30);
    TH1 * numprecorrection = new TH1D("numprecorrection","",30,0,30);
    TH1 * numafter = new TH1D("numafter","",30,0,30);
    TH1 * numlost = new TH1D("numlost","",10,0,10);
    TH1 * efficpt = new TH1D("efficpt", "", 100,0,30);
    TH1 * ccut2pt = new TH1D("ccut2pt", "", 100,0,30);
    TH1 * temp = new TH1D("temp","",300,-30,0);
    TH1 * ptofrands = new TH1D("ptofrands","",200,0,20);
    TH2 * numlostvdetpt = new TH2D("numlostvdetpt","num lost v. detector pt",100,0,30,20,0,20);
    //TESTS!!!
    TH1D * c_multpartraw = new TH1D("c_part mult", "c_part mult",100,0,300);
    TH1D * c_multpartfinalvisible = new TH1D("c_part mult final visible constituent", "c_part mult final visible constituent",100,0,300);
    TH1D * c_multuncut  = new TH1D("c_uncut","c_uncut",100,0,300);
    //TH1D * c_multcut    = new TH1D("c_cut","c_cut",100,0,300);
    TH1D * c_multcut2   = new TH1D("c_cut2","c_cut2",100,0,300);
    TH1D * c_multeffic  = new TH1D("c_effic","c_effic",100,0,300);
    TH1D * multpartraw = new TH1D("part mult", "part mult",100,0,300);
    TH1D * multpartfinalvisible = new TH1D("part mult final visible constituent", "part mult final visible constituent",100,0,300);
    TH1D * multuncut  = new TH1D("uncut","uncut",100,0,300);
    //TH1D * multcut    = new TH1D("cut","cut",100,0,300);
    TH1D * multcut2   = new TH1D("cut2","cut2",100,0,300);
    TH1D * multeffic  = new TH1D("effic","effic",100,0,300);
    
    ////////////////////////////////////////////////////////////////////////////
    
    //TH1D * proj1  = fit->ProjectionY("proj1", 1, 1);
    TH1D * proj2  = fit->ProjectionY("proj2", 2, 2);
    TH1D * proj3  = fit->ProjectionY("proj3", 3, 3);
    TH1D * proj4  = fit->ProjectionY("proj4", 4, 4);
    TH1D * proj5  = fit->ProjectionY("proj5", 5, 5);
    TH1D * proj6  = fit->ProjectionY("proj6", 6, 6);
    TH1D * proj7  = fit->ProjectionY("proj7", 7, 7);
    TH1D * proj8  = fit->ProjectionY("proj8", 8, 8);
    TH1D * proj9  = fit->ProjectionY("proj9", 9, 9);
    TH1D * proj10 = fit->ProjectionY("proj10", 10, 10);
    TH1D * proj11 = fit->ProjectionY("proj11", 11, 11);
    /*
    TH1D * proj12 = fit->ProjectionY("proj12", 12, 12);
    TH1D * proj13 = fit->ProjectionY("proj13", 13, 13);
    TH1D * proj14 = fit->ProjectionY("proj14", 14, 14);
    TH1D * proj15 = fit->ProjectionY("proj15", 15, 15);
    TH1D * proj16 = fit->ProjectionY("proj16", 16, 16);
    */
    
    ////////////////////////////////////////////////////////////////////////////
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    for (unsigned iEvent = 0; iEvent < NumEvents; ++iEvent) {
        test->Clear();
        
        if (!pythia.next()) continue;
        
        for (unsigned i = 0; i < pythia.event.size(); ++ i) {
            pid->Fill(pythia.event[i].id());
        }
        
        unsigned test_count = 0;
        analysis::add_particles(pythia.event, test, test_count, temp);
        
        ClusterSequence c_cut2(test->c_cut2_part,     jet_def);
        vector<PseudoJet> c_cut2_jets   = sorted_by_pt(c_cut2.inclusive_jets());
        
        if (c_cut2_jets.size() != 0) {
            if (jet_cuts(c_cut2_jets[0])){
                ccut2pt->Fill(c_cut2_jets[0].pt());
                vector<PseudoJet> c_cut2_cons = c_cut2_jets[0].constituents();
                for (unsigned j = 0; j < c_cut2_cons.size(); ++ j) {
                    c_cons_pt->Fill(c_cut2_cons[j].pt());
                }
            }
        }
    }
    
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    for (unsigned iEvent = 0; iEvent < NumEvents; ++iEvent) {
        test->Clear();
        
        if (!pythia.next()) continue;
        
        for (unsigned i = 0; i < pythia.event.size(); ++ i) {
            pid->Fill(pythia.event[i].id());
        }
        
        
        //std::cout << pythia.event.size() << "!!!!\n";
        
        // Find all final particles and fill PseudoJet vectors appropriately
        // SHOULD THINK ABOUT WHETHER PASSING IN ENTIRE 'TEST' IS TOO MEMORY INTENSIVE. IF SO, UNCOMMENT AFTER SEMICOLON AND FIX IN THE CLASS.
        
        //TEST!!!!!
        unsigned count_charged = 0;
        if (pythia.event.size() != 0) {
            multpartraw->Fill(pythia.event.size()/(double) 2);
            for (unsigned i = 0; i < pythia.event.size(); ++ i) {
                if (pythia.event[i].isCharged()) {
                    ++ count_charged;
                }
            }
        }
        c_multpartraw->Fill(count_charged/(double) 2);
        
        unsigned test_count = 0;
        c_multpartfinalvisible->Fill(analysis::add_particles(pythia.event, test, test_count, temp));/*test->uncut_part, test->c_uncut_part, test->cut_part, test->cut2_part,
                     test->c_cut_part, test->c_cut2_part, test->effic_part, test->c_effic_part);*/
        multpartfinalvisible->Fill(test_count);
        
        //WANT TO FIGURE OUT HOW TO WRITE THIS AT SOME POINT (NEED TO GIVE THE CLASS THE CLUSTERSEQUENCE, TOO)
        //analysis::cluster(test, jet_def);
        
        //have to cluster charged jets separately
        
        //~~~~~~~~~~~~~~~~~~~~~~~~CLUSTER AND GET JETS~~~~~~~~~~~~~~~~~~~~~~~~//
        
         ClusterSequence uncut(test->uncut_part,       jet_def);
         ClusterSequence cut(test->cut_part,           jet_def);
         ClusterSequence cut2(test->cut2_part,         jet_def);
         ClusterSequence effic(test->effic_part,       jet_def);
         vector<PseudoJet> uncut_jets    = sorted_by_pt(uncut.inclusive_jets());
         vector<PseudoJet> cut_jets      = sorted_by_pt(cut.inclusive_jets());
         vector<PseudoJet> cut2_jets     = sorted_by_pt(cut2.inclusive_jets());
         vector<PseudoJet> effic_jets    = sorted_by_pt(effic.inclusive_jets());
         
         ClusterSequence c_uncut(test->c_uncut_part,   jet_def);
         ClusterSequence c_cut(test->c_cut_part,       jet_def);
         ClusterSequence c_cut2(test->c_cut2_part,     jet_def);
         ClusterSequence c_effic(test->c_effic_part,   jet_def);
         vector<PseudoJet> c_uncut_jets  = sorted_by_pt(c_uncut.inclusive_jets());
         vector<PseudoJet> c_cut_jets    = sorted_by_pt(c_cut.inclusive_jets());
         vector<PseudoJet> c_cut2_jets   = sorted_by_pt(c_cut2.inclusive_jets());
         vector<PseudoJet> c_effic_jets  = sorted_by_pt(c_effic.inclusive_jets());
        
        //TEST
        /*
        unsigned counter = 0;
        for (unsigned i = 0; i < effic_jets.size(); ++ i) {
            if (jet_cuts(effic_jets[i])) {
                if (cut2_jets[i].user_index()) {
                    ++ counter;
                }
            }
        }
        if (counter == 0) {
            continue;
        }
         */
        //
        
        
        //cout << "Clustering with " << jet_def.description() << '\n';
        
        //~~~~~~~~~~~~~~~~~~~~~~FILL TREES WITH LEAD PTS~~~~~~~~~~~~~~~~~~~~~~//
        
        //leading jet pt spectra (total and only charged)
        //(a) without cuts, (b) with particle cuts, (c) with jet cuts, and (d) with cuts + efficiency
        
        //(a):
        if (uncut_jets.size() != 0) {
            test->uncut_leadPt = uncut_jets[0].pt();
        }
        
        if (c_uncut_jets.size() != 0) {
            test->c_uncut_leadPt = c_uncut_jets[0].pt();
        }
        
        //(b):
        if (cut_jets.size() != 0) {
            test->cut_leadPt = cut_jets[0].pt();
        }
        
        if (c_cut_jets.size() != 0) {
            test->c_cut_leadPt = c_cut_jets[0].pt();
        }
        
        //(c):
        if(cut2_jets.size() != 0) {
            if (jet_cuts(cut2_jets[0])) {
                //leading jet pt spectrum
                test->cut2_leadPt = cut2_jets[0].pt();
            }
        }
        
        if (c_cut2_jets.size() != 0 && jet_cuts(c_cut2_jets[0])) {
            test->c_cut2_leadPt = c_cut2_jets[0].pt();
        }
        
        //(d):
        if (effic_jets.size() != 0) {
            if (jet_cuts(effic_jets[0])) {
                test->effic_leadPt = effic_jets[0].pt();
            }
        }
        
        if (c_effic_jets.size() != 0 && jet_cuts(c_effic_jets[0])) {
            test->c_effic_leadPt = c_effic_jets[0].pt();
        }
        test->tlead->Fill();
        
        //~~~~~~~~~~~~~~~~FILL TREES WITH JETS AND CONSTITUENTS~~~~~~~~~~~~~~~~//
        
        //WRITE ASSIGN FUNCTION:
        //SOMETHING LIKE:
        //  for (i = 0, i < size of jet, ++ i) { assign(variable(s), jet); }
        //  where assign(variable(s), jet) goes through px, py, ... assigning them.
        //  STILL NOT SURE HOW THIS WOULD WORK - NO EASY WAY TO CALL FUNCTION
        
        //pseudojets without cuts
        for (unsigned i = 0; i < uncut_jets.size(); ++ i) {
            test->uncut_px    = uncut_jets[i].px();
            test->uncut_py    = uncut_jets[i].py();
            test->uncut_pz    = uncut_jets[i].pz();
            test->uncut_E     = uncut_jets[i].e();
            test->uncut_phi   = uncut_jets[i].phi();
            test->uncut_eta   = uncut_jets[i].eta();
            test->uncut_Pt    = uncut_jets[i].pt();
            
            test->tuncut->Fill();
            vector<PseudoJet> uncut_cons = uncut_jets[i].constituents();
            for (unsigned j = 0; j < uncut_cons.size(); ++ j) {
                test->uncut_cons_px   = uncut_cons[j].px();
                test->uncut_cons_py   = uncut_cons[j].py();
                test->uncut_cons_pz   = uncut_cons[j].pz();
                test->uncut_cons_E    = uncut_cons[j].e();
                test->uncut_cons_phi  = uncut_cons[j].phi();
                test->uncut_cons_eta  = uncut_cons[j].eta();
                test->uncut_cons_Pt   = uncut_cons[j].pt();
                
                test->tuncutcons->Fill();
            }
        }
        
        //pseudojets without cuts (charged)
        for (unsigned i = 0; i < c_uncut_jets.size(); ++ i) {
            test->c_uncut_px  = c_uncut_jets[i].px();
            test->c_uncut_py  = c_uncut_jets[i].py();
            test->c_uncut_pz  = c_uncut_jets[i].pz();
            test->c_uncut_E   = c_uncut_jets[i].e();
            test->c_uncut_phi = c_uncut_jets[i].phi();
            test->c_uncut_eta = c_uncut_jets[i].eta();
            test->c_uncut_Pt  = c_uncut_jets[i].pt();
            
            test->tcuncut->Fill();
            
            vector<PseudoJet> c_uncut_cons = c_uncut_jets[i].constituents();
            for (unsigned j = 0; j < c_uncut_cons.size(); ++ j) {
                test->c_uncut_cons_px     = c_uncut_cons[j].px();
                test->c_uncut_cons_py     = c_uncut_cons[j].py();
                test->c_uncut_cons_pz     = c_uncut_cons[j].pz();
                test->c_uncut_cons_E      = c_uncut_cons[j].e();
                test->c_uncut_cons_phi    = c_uncut_cons[j].phi();
                test->c_uncut_cons_eta    = c_uncut_cons[j].eta();
                test->c_uncut_cons_Pt     = c_uncut_cons[j].pt();
                
                test->tcuncutcons->Fill();
            }
        }
        
        //pseudojets with particle level cuts
        for (unsigned i = 0; i < cut_jets.size(); ++ i) {
            test->cut_px  = cut_jets[i].px();
            test->cut_py  = cut_jets[i].py();
            test->cut_pz  = cut_jets[i].pz();
            test->cut_E   = cut_jets[i].e();
            test->cut_phi = cut_jets[i].phi();
            test->cut_eta = cut_jets[i].eta();
            test->cut_Pt  = cut_jets[i].pt();
            
            test->tcut->Fill();
        }
        
        for (unsigned i = 0; i < c_cut_jets.size(); ++ i) {
            test->c_cut_px    = c_cut_jets[i].px();
            test->c_cut_py    = c_cut_jets[i].py();
            test->c_cut_pz    = c_cut_jets[i].pz();
            test->c_cut_E     = c_cut_jets[i].e();
            test->c_cut_phi   = c_cut_jets[i].phi();
            test->c_cut_eta   = c_cut_jets[i].eta();
            test->c_cut_Pt    = c_cut_jets[i].pt();
            
            test->tccut->Fill();
        }
        
        //pseudojets with both sets of cuts
        for (unsigned i = 0; i < cut2_jets.size(); ++ i) {
            //pseudojet axis cut
            if (jet_cuts(cut2_jets[i])) {
                test->cut2_px     = cut2_jets[i].px();
                test->cut2_py     = cut2_jets[i].py();
                test->cut2_pz     = cut2_jets[i].pz();
                test->cut2_E      = cut2_jets[i].e();
                test->cut2_phi    = cut2_jets[i].phi();
                test->cut2_eta    = cut2_jets[i].eta();
                test->cut2_Pt     = cut2_jets[i].pt();
                
                test->tcut2->Fill();
            
            
                vector<PseudoJet> cut2_cons = cut2_jets[i].constituents();
                for (unsigned j = 0; j < cut2_cons.size(); ++ j) {
                    test->cut2_cons_px     = cut2_cons[j].px();
                    test->cut2_cons_py     = cut2_cons[j].py();
                    test->cut2_cons_pz     = cut2_cons[j].pz();
                    test->cut2_cons_E      = cut2_cons[j].e();
                    test->cut2_cons_phi    = cut2_cons[j].phi();
                    test->cut2_cons_eta    = cut2_cons[j].eta();
                    test->cut2_cons_Pt     = cut2_cons[j].pt();
                
                    test->tcut2cons->Fill();
                }
                
            }
        }
        
        if (c_cut2_jets.size() != 0) {
            if (jet_cuts(c_cut2_jets[0])) {
                ptbefore->Fill(c_cut2_jets[0].pt());
                numbefore->Fill(c_cut2_jets[0].constituents().size());
            }
        }
        
        for (unsigned i = 0; i < c_cut2_jets.size(); ++ i) {
            if (jet_cuts(c_cut2_jets[i])) {
                test->c_cut2_px   = c_cut2_jets[i].px();
                test->c_cut2_py   = c_cut2_jets[i].py();
                test->c_cut2_pz   = c_cut2_jets[i].pz();
                test->c_cut2_E    = c_cut2_jets[i].e();
                test->c_cut2_phi  = c_cut2_jets[i].phi();
                test->c_cut2_eta  = c_cut2_jets[i].eta();
                test->c_cut2_Pt   = c_cut2_jets[i].pt();
                
                test->tccut2->Fill();
                //ccut2pt->Fill(c_cut2_jets[0].pt());
                
                //unsigned num_loss = nb->GetRandom();
                //Double_t pt_loss = 0;
            
                vector<PseudoJet> c_cut2_cons = c_cut2_jets[i].constituents();
                for (unsigned j = 0; j < c_cut2_cons.size(); ++ j) {
                    test->c_cut2_cons_px     = c_cut2_cons[j].px();
                    test->c_cut2_cons_py     = c_cut2_cons[j].py();
                    test->c_cut2_cons_pz     = c_cut2_cons[j].pz();
                    test->c_cut2_cons_E      = c_cut2_cons[j].e();
                    test->c_cut2_cons_phi    = c_cut2_cons[j].phi();
                    test->c_cut2_cons_eta    = c_cut2_cons[j].eta();
                    test->c_cut2_cons_Pt     = c_cut2_cons[j].pt();
                
                    //c_cons_pt->Fill(c_cut2_cons[j].pt());
                }
                /*
                for (unsigned j = 0; j < num_loss; ++ j) {
                    pt_loss += c_cons_pt->GetRandom();
                }*/
                //ptlost->Fill(pt_loss);
                //ptafter->Fill(c_cut2_jets[i].pt() - pt_loss);
                test->tccut2cons->Fill();
            }
        }
    
        
        //pseudojets with cuts + charged particle tracking efficiency selection
        for (unsigned i = 0; i < effic_jets.size(); ++ i) {
            if (jet_cuts(effic_jets[i])) {
                test->effic_px    = effic_jets[i].px();
                test->effic_py    = effic_jets[i].py();
                test->effic_pz    = effic_jets[i].pz();
                test->effic_E     = effic_jets[i].e();
                test->effic_phi   = effic_jets[i].phi();
                test->effic_eta   = effic_jets[i].eta();
                test->effic_Pt    = effic_jets[i].pt();
                
                test->teffic->Fill();
            
            
                vector<PseudoJet> effic_cons = effic_jets[i].constituents();
                for (unsigned j = 0; j < effic_cons.size(); ++ j) {
                    test->effic_cons_px     = effic_cons[j].px();
                    test->effic_cons_py     = effic_cons[j].py();
                    test->effic_cons_pz     = effic_cons[j].pz();
                    test->effic_cons_E      = effic_cons[j].e();
                    test->effic_cons_phi    = effic_cons[j].phi();
                    test->effic_cons_eta    = effic_cons[j].eta();
                    test->effic_cons_Pt     = effic_cons[j].pt();
                    
                    test->tefficcons->Fill();
                }
            }
        }
        
        for (unsigned i = 0; i < c_effic_jets.size(); ++ i) {
            if (jet_cuts(c_effic_jets[i])) {
                test->c_effic_px  = c_effic_jets[i].px();
                test->c_effic_py  = c_effic_jets[i].py();
                test->c_effic_pz  = c_effic_jets[i].pz();
                test->c_effic_E   = c_effic_jets[i].e();
                test->c_effic_phi = c_effic_jets[i].phi();
                test->c_effic_eta = c_effic_jets[i].eta();
                test->c_effic_Pt  = c_effic_jets[i].pt();
                
                test->tceffic->Fill();

                
                ////
                if (i == 0) {
                    //efficpt->Fill(c_cut2_jets[i].pt());
                    
                    /////
                    /*
                    unsigned num_loss;
                    char name[16];
                    if (c_effic_jets[0].constituents().size() == 1) {
                        sprintf(name, "%s%d","proj",2);
                    }
                    if (c_effic_jets[0].constituents().size() >= 11) {
                        sprintf(name, "%s%d","proj",11);
                    }
                    else {
                        sprintf(name, "%s%d","proj",c_effic_jets[0].constituents().size()); //SHOULD THIS BE +1 OR NOT? BINNING WEIRDNESS UGH // THINK THAT MADE IT WORSE
                    }
                    TH1D *h = (TH1D*)gDirectory->Get(name);
                    
                    if ( h ) { // make sure the Get succeeded
                        if (h->GetEntries() > 0) {
                            num_loss = (unsigned)h->GetRandom();
                        }
                        else {
                            //cout << "whoops\n";
                            //cout << name << '\n';
                            num_loss = 1;
                        }
                    }
                    else {cout << "get failed: " << name << " " << c_effic_jets[0].constituents().size() << '\n'; num_loss = 0;}
                    */
                    
                    /////
                    
                    unsigned num_loss = nbnewtest->GetRandom();
                    /*
                    unsigned num_loss = nb->GetRandom();
                     */
                    Double_t pt_loss = 0;
                    for (unsigned j = 0; j < num_loss; ++ j) {
                        Double_t cconspt = c_cons_pt->GetRandom();
                        ptofrands->Fill(cconspt);
                        pt_loss += cconspt;//c_cons_pt->GetRandom();
                    }
                    ptlost->Fill(-pt_loss);
                    numlost->Fill(num_loss);
                    ptprecorrection->Fill(c_effic_jets[i].pt());
                    numlostvdetpt->Fill(c_effic_jets[i].pt(), num_loss);
                    numprecorrection->Fill(c_effic_jets[i].constituents().size());
                    ptafter->Fill(c_effic_jets[i].pt() + pt_loss);
                    numafter->Fill(c_effic_jets[i].constituents().size() + num_loss);
                }
                ////
            
                
                vector<PseudoJet> c_effic_cons = c_effic_jets[i].constituents();
                for (unsigned j = 0; j < c_effic_cons.size(); ++ j) {
                    test->c_effic_cons_px     = c_effic_cons[j].px();
                    test->c_effic_cons_py     = c_effic_cons[j].py();
                    test->c_effic_cons_pz     = c_effic_cons[j].pz();
                    test->c_effic_cons_E      = c_effic_cons[j].e();
                    test->c_effic_cons_phi    = c_effic_cons[j].phi();
                    test->c_effic_cons_eta    = c_effic_cons[j].eta();
                    test->c_effic_cons_Pt     = c_effic_cons[j].pt();
                
                    test->tcefficcons->Fill();
                }
            }
        }
        /*
        //FILL 2D FIT HISTO:
        if (c_effic_jets.size() != 0 && c_cut2_jets.size() != 0) {
            fit->Fill(c_effic_jets[0].constituents().size(), c_cut2_jets[0].constituents().size() - c_effic_jets[0].constituents().size());
        }
        */
        
        //Checking that pt spectrum without 2 highest-pt jets looks right
        analysis::without_lead2_jets(effic_jets, test->without_leadsublead, test->ttests);
        
        //Checks pt & num difference between leading jet with efficiency correction and
        //geometrically closest jet without efficiency correction.
        analysis::geometric_diff(effic_jets, cut2_jets, test->ptdiff, test->num_diff, test->num_before, test->num_after, test->rel_diff, count);
        analysis::geometric_diff(c_effic_jets, c_cut2_jets, test->c_ptdiff, test->c_num_diff, test->c_num_before, test->c_num_after, test->c_rel_diff, c_count);
        
        //TEST!!!!
        
        if (c_uncut_jets.size() != 0) {
            c_multuncut->Fill(c_uncut_jets[0].constituents().size()/(double) 2);
        }
        /*
        if (c_cut_jets.size() != 0) {
            c_multcut->Fill(c_cut_jets[0].constituents().size()/(double) 2);
        }
        */
        if (c_cut2_jets.size() != 0) {
            if (jet_cuts(c_cut2_jets[0])){
                c_multcut2->Fill(c_cut2_jets[0].constituents().size()/(double) 2);
            }
        }
        
        if (c_effic_jets.size() != 0) {
            if (jet_cuts(c_effic_jets[0])){
                c_multeffic->Fill(c_effic_jets[0].constituents().size()/(double) 2);
            }
        }
        
        //~~neutral + charged
        
        if (uncut_jets.size() != 0) {
            multuncut->Fill(uncut_jets[0].constituents().size()/(double) 2);
        }
        /*
         if (cut_jets.size() != 0) {
         multcut->Fill(cut_jets[0].constituents().size()/(double) 2);
         }
         */
        if (cut2_jets.size() != 0) {
            if (jet_cuts(cut2_jets[0])){
                multcut2->Fill(cut2_jets[0].constituents().size()/(double) 2);
            }
        }
        
        if (effic_jets.size() != 0) {
            if (jet_cuts(effic_jets[0])){
                multeffic->Fill(effic_jets[0].constituents().size()/(double) 2);
            }
        }
        
        ////////
        
        test->tdiffs->Fill();
    // End of event loop
    }
    
    temp->Scale(1/(double) temp->GetEntries());
    ptlost->Scale(1/(double) ptlost->GetEntries());

    
    std::cout << "Failed " << (count/(double) NumEvents)*100 << "% of the time\n";
    std::cout << "Failed " << (c_count/(double) NumEvents)*100 << "% of the time (charged)\n";
    
    pythia.stat();
    
    
    cout << "Writing to " << fout->GetName() << '\n';
    //tout->Print(); if status of a particular tree is desired.
    //save trees to output root file
    test->write();
    /*fit2d->Write();*/fout->Write();
    
    /*fit2d->Close();*///fout->Close();
    
    return 0;
}
