//Program modified from main01.cc in pythia8219 examples folder
//Runs events in a narrow pT hard bin, clusters the final state
//particles with anti-kT, applies some non-stringent cuts,
//and plots histograms/fills trees with & without a tracking
//efficiency selection cut.
//Isaac Mooney - WSU, March 2017 

#include "Pythia8/Pythia.h"
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include <iostream>
#include <sstream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <random>
#include <chrono>

using namespace fastjet;
using namespace Pythia8;
using namespace std;

int main() {
  // Generator. Process selection. STAR initialization. Histograms.
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  //seeding pseudorandom number generator
  auto begin = std::chrono::high_resolution_clock::now();

  //generating the distribution
  std::random_device rd;
  std::mt19937 g(rd());
  auto end = std::chrono::high_resolution_clock::now();
  g.seed(std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  Pythia pythia;
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("PhaseSpace:pTHatMax = 30.");
  pythia.init();

  TFile *fout = new TFile("main01out.root","RECREATE");  // Create output file
  TTree *tout = new TTree("pseudojets","pseudojets"); // Create an output tree

  //create branches on the tree
  double px, py, pz, E, phi, eta, pt, con_pt;
  tout->Branch("px", &px);
  tout->Branch("py", &py);
  tout->Branch("pz", &pz);
  tout->Branch("E", &E);
  tout->Branch("phi", &phi);
  tout->Branch("eta", &eta);
  tout->Branch("pt", &pt);
  tout->Branch("con_pt", &con_pt);

  //create histograms

    //without cuts (jet level)
  TH1D* uncut_px = new TH1D("u_px","px",101,-3,3);
  TH1D* uncut_py = new TH1D("u_py","py",101,-3,3);
  TH1D* uncut_pz = new TH1D("u_pz","pz",101,-3,3);
  TH1D* uncut_E = new TH1D("u_E","E",101,-1,10);
  TH1D* uncut_phi = new TH1D("u_phi", "phi", 100, 0, 6.3);
  TH1D* uncut_eta = new TH1D("u_eta", "eta", 120, -3, 3);
  TH1D* uncut_leadPt = new TH1D("u_leadPt", "pt", 100, 0, 50);
  TH1D* uncut_Pt = new TH1D("u_Pt", "pt", 100, 0, 40);

  //charged, without cuts (jet level)                                                                                                                                                                
  TH1D* c_uncut_px = new TH1D("c_u_px","px",101,-3,3);
  TH1D* c_uncut_py = new TH1D("c_u_py","py",101,-3,3);
  TH1D* c_uncut_pz = new TH1D("c_u_pz","pz",101,-3,3);
  TH1D* c_uncut_E = new TH1D("c_u_E","E",101,-1,10);
  TH1D* c_uncut_phi = new TH1D("c_u_phi", "phi", 100, 0, 6.3);
  TH1D* c_uncut_eta = new TH1D("c_u_eta", "eta", 120, -3, 3);
  TH1D* c_uncut_leadPt = new TH1D("c_u_leadPt", "pt", 100, 0, 50);
  TH1D* c_uncut_Pt = new TH1D("c_u_Pt", "pt", 100, 0, 40);

  //without cuts (constituent level)                                                                                                                                      
  TH1D* uncut_cons_px = new TH1D("u_cons_px","px",101,-3,3);
  TH1D* uncut_cons_py = new TH1D("u_cons_py","py",101,-3,3);
  TH1D* uncut_cons_pz = new TH1D("u_cons_pz","pz",101,-3,3);
  TH1D* uncut_cons_E = new TH1D("u_cons_E","E",101,-1,10);
  TH1D* uncut_cons_phi = new TH1D("u_cons_phi", "phi", 100, 0, 6.3);
  TH1D* uncut_cons_eta = new TH1D("u_cons_eta", "eta", 120, -3, 3);
  TH1D* uncut_cons_Pt = new TH1D("u_cons_Pt", "pt", 100, 0, 40);

  //charged, without cuts (constituent level)                                                                                                                                                          
  TH1D* c_uncut_cons_px = new TH1D("c_u_cons_px","px",101,-3,3);
  TH1D* c_uncut_cons_py = new TH1D("c_u_cons_py","py",101,-3,3);
  TH1D* c_uncut_cons_pz = new TH1D("c_u_cons_pz","pz",101,-3,3);
  TH1D* c_uncut_cons_E = new TH1D("c_u_cons_E","E",101,-1,10);
  TH1D* c_uncut_cons_phi = new TH1D("c_u_cons_phi", "phi", 100, 0, 6.3);
  TH1D* c_uncut_cons_eta = new TH1D("c_u_cons_eta", "eta", 120, -3, 3);
  TH1D* c_uncut_cons_Pt = new TH1D("c_u_cons_Pt", "pt", 100, 0, 40);

    //with particle cuts
    TH1D* cut_px = new TH1D("c_px","px",101,-3,3);
    TH1D* cut_py = new TH1D("c_py","py",101,-3,3);
    TH1D* cut_pz = new TH1D("c_pz","pz",101,-3,3);
    TH1D* cut_E = new TH1D("c_E","E",101,-1,10);
    TH1D* cut_phi = new TH1D("c_phi", "phi", 100, 0, 6.3);
    TH1D* cut_eta = new TH1D("c_eta", "eta", 120, -3, 3);
    TH1D* cut_leadPt = new TH1D("c_leadPt", "pt", 100, 0, 50);
    TH1D* cut_Pt = new TH1D("c_Pt", "pt", 100, 0, 40);

    //with jet cuts
    TH1D* cut2_px = new TH1D("c2_px","px",101,-3,3);
    TH1D* cut2_py = new TH1D("c2_py","py",101,-3,3);
    TH1D* cut2_pz = new TH1D("c2_pz","pz",101,-3,3);
    TH1D* cut2_E = new TH1D("c2_E","E",101,-1,10);
    TH1D* cut2_phi = new TH1D("c2_phi", "phi", 100, 0, 6.3);
    TH1D* cut2_eta = new TH1D("c2_eta", "eta", 120, -3, 3);
    TH1D* cut2_leadPt = new TH1D("c2_leadPt", "pt", 100, 0, 50);
    TH1D* cut2_Pt = new TH1D("c2_Pt", "pt", 100, 0, 40);

    //with tracking efficiency selection applied on top of cuts - jet level
    TH1D* effic_px = new TH1D("e_px","px",101,-3,3);
    TH1D* effic_py = new TH1D("e_py","py",101,-3,3);
    TH1D* effic_pz = new TH1D("e_pz","pz",101,-3,3);
    TH1D* effic_E = new TH1D("e_E","E",101,-1,10);
    TH1D* effic_phi = new TH1D("e_phi", "phi", 100, 0, 6.3);
    TH1D* effic_eta = new TH1D("e_eta", "eta", 120, -3, 3);
    TH1D* effic_leadPt = new TH1D("e_leadPt", "pt", 100, 0, 50);
    TH1D* effic_Pt = new TH1D("e_Pt", "pt", 100, 0, 40);
    TH1D* test_ptdiff = new TH1D("e_ptdiff", "pt", 100, 0, 10);
    TH1D* test_without_leadsublead = new TH1D("e_leadsublead", "pt", 100, 0, 10);
    
    //charged (denoted c_) fraction only, with tracking efficiency selection applied on top of cuts - jet level                                                                                                                   
    TH1D* c_effic_px = new TH1D("c_e_px","px",101,-3,3);
    TH1D* c_effic_py = new TH1D("c_e_py","py",101,-3,3);
    TH1D* c_effic_pz = new TH1D("c_e_pz","pz",101,-3,3);
    TH1D* c_effic_E = new TH1D("c_e_E","E",101,-1,10);
    TH1D* c_effic_phi = new TH1D("c_e_phi", "phi", 100, 0, 6.3);
    TH1D* c_effic_eta = new TH1D("c_e_eta", "eta", 120, -3, 3);
    TH1D* c_effic_leadPt = new TH1D("c_e_leadPt", "pt", 100, 0, 50);
    TH1D* c_effic_Pt = new TH1D("c_e_Pt", "pt", 100, 0, 40);

  //TH1D* orig = new TH1D("orig", "no efficiency correction", 100, 0, 5);
      
    //an event worth of pseudojets
    vector<PseudoJet> uncut_part;
    vector<PseudoJet> cut_part;
    vector<PseudoJet> cut2_part;
    vector<PseudoJet> effic_part;
  
    //Hist mult("charged multiplicity", 100, -0.5, 799.5);
  
    // choose a jet definition
    double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);

    //efficiency (constant for now)
    //double efficiency = 0.8;
  
    // Begin event loop. Generate event. Skip if error. List first one.
    
    for (int iEvent = 0; iEvent < 10000; ++iEvent) {
      uncut_part.clear();
      cut_part.clear();
      cut2_part.clear();
      effic_part.clear();
      if (!pythia.next()) continue;
      //int nCharged = 0;
      // Find all final particles and fill histograms
      for (int i = 0; i < pythia.event.size(); ++i) {
        //only select final state particles within particle-level eta cut
        if (pythia.event[i].isFinal() && pythia.event[i].isVisible()) {
	  //save particle charge info for later in fastjet
	  //set_user_index(pythia.event[i].charge())
	  PseudoJet current = PseudoJet(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
          current.set_user_index(pythia.event[i].charge());  
	  uncut_part.push_back(current);
            if (abs(pythia.event[i].eta()) < 1) {
                //++nCharged;
	      double effic_num = dis(g);
	      if (effic_num > 0.8 && pythia.event[i].isCharged()) {
                    //fill the original jet, skip filling the tracking efficiency applied jet
                    cut_part.push_back(current);
                    cut2_part.push_back(current);
                    continue;
	      }
              else {
		//fill both
		cut_part.push_back(current);
		cut2_part.push_back(current);
		effic_part.push_back(current);
	      }
            }
        }
      }

      // run the clustering, extract the jets
      ClusterSequence uncut(uncut_part, jet_def);
      ClusterSequence cut2(cut2_part, jet_def);
      ClusterSequence effic(effic_part, jet_def);
      vector<PseudoJet> uncut_jets = sorted_by_pt(uncut.inclusive_jets());
      vector<PseudoJet> cut2_jets = sorted_by_pt(cut2.inclusive_jets());
      vector<PseudoJet> effic_jets = sorted_by_pt(effic.inclusive_jets());

      //cout << "Clustering with " << jet_def.description() << '\n';

      //leading jet pt spectra (a) without cuts, (b) with particle cuts, (c) with jet cuts, and (d) with cuts + efficiency
      //(a):
      uncut_leadPt->Fill(uncut_jets[0].pt());
      //(b):
      cut_leadPt->Fill(cut_part[0].pt());
      //(c):
      if (abs(cut2_jets[0].eta()) < (1 - R)) {
        //leading jet pt spectrum
        pt = cut2_jets[0].pt();
        cut2_leadPt->Fill(cut2_jets[0].pt());
        tout->Fill();
      }
      //(d):
      if (abs(effic_jets[0].eta()) < (1 - R)) {
        effic_leadPt->Fill(effic_jets[0].pt());
      }    

      //pseudojets without cuts
      for (int i = 0; i < uncut_jets.size(); ++ i) {
        uncut_px->Fill(uncut_jets[i].px());
        uncut_py->Fill(uncut_jets[i].py());
        uncut_pz->Fill(uncut_jets[i].pz());
        uncut_E->Fill(uncut_jets[i].e());
        uncut_phi->Fill(uncut_jets[i].phi());
        uncut_eta->Fill(uncut_jets[i].eta());
        uncut_Pt->Fill(uncut_jets[i].pt());
	if (effic_jets[i].user_index() != 0) {
	  c_uncut_px->Fill(uncut_jets[i].px());
	  c_uncut_py->Fill(uncut_jets[i].py());
	  c_uncut_pz->Fill(uncut_jets[i].pz());
	  c_uncut_E->Fill(uncut_jets[i].e());
	  c_uncut_phi->Fill(uncut_jets[i].phi());
	  c_uncut_eta->Fill(uncut_jets[i].eta());
	  c_uncut_Pt->Fill(uncut_jets[i].pt());
	}
      
	vector<PseudoJet> uncut_cons = uncut_jets[i].constituents(); 
	for (int j = 0; j < uncut_cons.size(); ++ j) {
	  uncut_cons_px->Fill(uncut_cons[j].px());
	  uncut_cons_py->Fill(uncut_cons[j].py());
	  uncut_cons_pz->Fill(uncut_cons[j].pz());
	  uncut_cons_E->Fill(uncut_cons[j].e());
	  uncut_cons_phi->Fill(uncut_cons[j].phi());
	  uncut_cons_eta->Fill(uncut_cons[j].eta());
	  uncut_cons_Pt->Fill(uncut_cons[j].pt());
	  if (effic_jets[i].user_index() != 0) {
	    c_uncut_cons_px->Fill(uncut_cons[j].px());
	    c_uncut_cons_py->Fill(uncut_cons[j].py());
	    c_uncut_cons_pz->Fill(uncut_cons[j].pz());
	    c_uncut_cons_E->Fill(uncut_cons[j].e());
	    c_uncut_cons_phi->Fill(uncut_cons[j].phi());
	    c_uncut_cons_eta->Fill(uncut_cons[j].eta());
	    c_uncut_cons_Pt->Fill(uncut_cons[j].pt());
	  }
	}
      }
      
      //pseudojets with particle level cuts
      for (int i = 0; i < cut_part.size(); ++ i) {
        cut_px->Fill(cut_part[i].px());
        cut_py->Fill(cut_part[i].py());
        cut_pz->Fill(cut_part[i].pz());
        cut_E->Fill(cut_part[i].e());
        cut_phi->Fill(cut_part[i].phi());
        cut_eta->Fill(cut_part[i].eta());
        cut_Pt->Fill(cut_part[i].pt());
      }
      
      //pseudojets with both sets of cuts
      for (int i = 0; i < cut2_jets.size(); ++ i) {
	//pseudojet axis cut
	if (abs(cut2_jets[i].eta()) < (1 - R)) {
          px = cut2_jets[i].px();
          py = cut2_jets[i].py();
          pz = cut2_jets[i].pz();
          E = cut2_jets[i].e();
          phi = cut2_jets[i].phi();
          eta = cut2_jets[i].eta();
          con_pt = cut2_jets[i].pt();
          tout->Fill();
          
          cut2_px->Fill(cut2_jets[i].px());
          cut2_py->Fill(cut2_jets[i].py());
          cut2_pz->Fill(cut2_jets[i].pz());
          cut2_E->Fill(cut2_jets[i].e());
          cut2_phi->Fill(cut2_jets[i].phi());
          cut2_eta->Fill(cut2_jets[i].eta());
          cut2_Pt->Fill(cut2_jets[i].pt());
	}
      }
    
      //pseudojets with cuts + charged particle tracking efficiency selection
      for (int i = 0; i < effic_jets.size(); ++ i) {
        if (abs(effic_jets[i].eta()) < (1 - R)) {
	  effic_px->Fill(effic_jets[i].px());
	  effic_py->Fill(effic_jets[i].py());
	  effic_pz->Fill(effic_jets[i].pz());
	  effic_E->Fill(effic_jets[i].e());
	  effic_phi->Fill(effic_jets[i].phi());
	  effic_eta->Fill(effic_jets[i].eta());
	  effic_Pt->Fill(effic_jets[i].pt());
	  if (effic_jets[i].user_index() != 0) {
	    c_effic_px->Fill(effic_jets[i].px());
	    c_effic_py->Fill(effic_jets[i].py());
	    c_effic_pz->Fill(effic_jets[i].pz());
	    c_effic_E->Fill(effic_jets[i].e());
	    c_effic_phi->Fill(effic_jets[i].phi());
	    c_effic_eta->Fill(effic_jets[i].eta());
	    c_effic_Pt->Fill(effic_jets[i].pt());
	  }
        }
	if (i != 0 && i != 1) {
	  test_without_leadsublead->Fill(effic_jets[i].pt());
	}
      }
      
      //testing!
      if (abs(uncut_jets[0].eta() - effic_jets[0].eta()) < 0.0001) {
	double ptdiff = uncut_jets[0].pt() - effic_jets[0].pt();
	test_ptdiff->Fill(ptdiff);
      }
      /*for (int i = 0; i < uncut_jets.size(); ++ i) {
	if (uncut_jets[i].eta() == effic_jets[0].eta()) {
	  test_ptdiff->Fill(uncut_jets[i].pt() - effic_jets[0].pt());
	  break;
	}
	}*/

      //mult.fill( nCharged );
    
      // End of event loop. Statistics. Histogram. Done.
    }
    pythia.stat();
    //cout << mult << '\n';
  
    cout << "Writing to " << fout->GetName() << '\n';
    //save histograms to output root file
    //tout->Print();
    tout->Write();
    fout->Write();
    fout->Close();

    return 0;
}
