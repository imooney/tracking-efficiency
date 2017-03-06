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

using namespace fastjet;
using namespace Pythia8;
using namespace std;

int main() {
  // Generator. Process selection. STAR initialization. Histograms.
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
    //without cuts
  TH1D* uncut_px = new TH1D("u_px","px",101,-3,3);
  TH1D* uncut_py = new TH1D("u_py","py",101,-3,3);
  TH1D* uncut_pz = new TH1D("u_pz","pz",101,-3,3);
  TH1D* uncut_E = new TH1D("u_E","E",101,-1,10);
  TH1D* uncut_phi = new TH1D("u_phi", "phi", 100, 0, 6.3);
  TH1D* uncut_eta = new TH1D("u_eta", "eta", 100, -3.15, 3.15);
  TH1D* uncut_leadPt = new TH1D("u_leadPt", "pt", 100, 0, 50);
  TH1D* uncut_consPt = new TH1D("u_consPt", "constituent pt", 100, 0, 5);
    //with cuts
    TH1D* cut_px = new TH1D("c_px","px",101,-3,3);
    TH1D* cut_py = new TH1D("c_py","py",101,-3,3);
    TH1D* cut_pz = new TH1D("c_pz","pz",101,-3,3);
    TH1D* cut_E = new TH1D("c_E","E",101,-1,10);
    TH1D* cut_phi = new TH1D("c_phi", "phi", 100, 0, 6.3);
    TH1D* cut_eta = new TH1D("c_eta", "eta", 100, -3.15, 3.15);
    TH1D* cut_leadPt = new TH1D("c_leadPt", "pt", 100, 0, 50);
    TH1D* cut_consPt = new TH1D("c_consPt", "constituent pt", 100, 0, 5);
    //with tracking efficiency selection applied on top of cuts
    TH1D* effic_px = new TH1D("e_px","px",101,-3,3);
    TH1D* effic_py = new TH1D("e_py","py",101,-3,3);
    TH1D* effic_pz = new TH1D("e_pz","pz",101,-3,3);
    TH1D* effic_E = new TH1D("e_E","E",101,-1,10);
    TH1D* effic_phi = new TH1D("e_phi", "phi", 100, 0, 6.3);
    TH1D* effic_eta = new TH1D("e_eta", "eta", 100, -3.15, 3.15);
    TH1D* effic_leadPt = new TH1D("e_leadPt", "pt", 100, 0, 50);
    TH1D* effic_consPt = new TH1D("e_consPt", "constituent pt", 100, 0, 5);
    
  //TH1D* orig = new TH1D("orig", "no efficiency correction", 100, 0, 5);
      
    //an event worth of pseudojets
    vector<PseudoJet> uncut_part;
    vector<PseudoJet> cut_part;
    vector<PseudoJet> effic_part;
  
  //Hist mult("charged multiplicity", 100, -0.5, 799.5);
  
  // choose a jet definition
  double R = 0.4;
  JetDefinition jet_def(antikt_algorithm, R);
  
  // Begin event loop. Generate event. Skip if error. List first one.
    
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
      uncut_part.clear();
      cut_part.clear();
      effic_part.clear();
      if (!pythia.next()) continue;
      //int nCharged = 0;
      int counter = 0;
      // Find all final particles and fill histograms
      for (int i = 0; i < pythia.event.size(); ++i) {
        //only select final state particles within particle-level eta cut
        if (pythia.event[i].isFinal()) {
            uncut_part.push_back(PseudoJet(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
            if (abs(pythia.event[i].eta()) < 1) {
                //++nCharged;
                //if it's the 5th charged particle, exclude it (so we scale charged particles in histogram by 0.8)
                if (counter == 4 && pythia.event[i].isCharged()) {
                    //reset the counter for later
                    counter = 0;
                    //fill the original jet, skip filling the tracking efficiency applied jet
                    cut_part.push_back(PseudoJet(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
                    continue;
                }
                else {
                    //fill both
                    cut_part.push_back(PseudoJet(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
                    effic_part.push_back(PseudoJet(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
                }
                //every charged particle increments the counter
                if (pythia.event[i].isCharged()) {
                    counter ++;
                }
            }
        }
      }
    
    // run the clustering, extract the jets
    ClusterSequence uncut(uncut_part, jet_def);
    ClusterSequence cut(cut_part, jet_def);
    ClusterSequence effic(effic_part, jet_def);
    vector<PseudoJet> uncut_jets = sorted_by_pt(uncut.inclusive_jets());
    vector<PseudoJet> cut_jets = sorted_by_pt(cut.inclusive_jets());
    vector<PseudoJet> effic_jets = sorted_by_pt(effic.inclusive_jets());

    //cout << "Clustering with " << jet_def.description() << '\n';

    //leading jet pt spectra (a) without cuts, (b) with cuts, and (c) with cuts + efficiency
    //(a):
      uncut_leadPt->Fill(uncut_jets[0].pt());
    //(b):
    if (abs(cut_jets[0].eta()) < (1 - R)) {
        //leading jet pt spectrum
        pt = cut_jets[0].pt();
        cut_leadPt->Fill(cut_jets[0].pt());
        tout->Fill();
    }
    //(c):
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
        uncut_consPt->Fill(uncut_jets[i].pt());
      }
      
    //pseudojets with cuts
    for (int i = 0; i < cut_jets.size(); ++ i) {
      //pseudojet axis cut
      if (abs(cut_jets[i].eta()) < (1 - R)) {
          px = cut_jets[i].px();
          py = cut_jets[i].py();
          pz = cut_jets[i].pz();
          E = cut_jets[i].e();
          phi = cut_jets[i].phi();
          eta = cut_jets[i].eta();
          con_pt = cut_jets[i].pt();
          tout->Fill();
          
          cut_px->Fill(cut_jets[i].px());
          cut_py->Fill(cut_jets[i].py());
          cut_pz->Fill(cut_jets[i].pz());
          cut_E->Fill(cut_jets[i].e());
          cut_phi->Fill(cut_jets[i].phi());
          cut_eta->Fill(cut_jets[i].eta());
          cut_consPt->Fill(cut_jets[i].pt());
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
            effic_consPt->Fill(effic_jets[i].pt());
        }
    }
      
      /*
    //pseudojets without charged particle tracking efficiency selection
    for (int i = 0; i < orig_jets.size(); ++ i) {
      //pseudojet axis cut
      if (abs(orig_jets[i].eta()) < (1 - R)) {
          orig->Fill(orig_jets[i].pt());
      }
    }
       */
    //mult.fill( nCharged );
    
  // End of event loop. Statistics. Histogram. Done.
  }
    pythia.stat();
    //cout << mult << '\n';
  
    cout << "Writing to " << fout->GetName() << '\n';
    //save histograms to output root file
    //tout->Print();
    tout->Write();
    uncut_px->Write();
    uncut_py->Write();
    uncut_pz->Write();
    uncut_E->Write();
    uncut_phi->Write();
    uncut_eta->Write();
    uncut_leadPt->Write();
    uncut_consPt->Write();
    
    cut_px->Write();
    cut_py->Write();
    cut_pz->Write();
    cut_E->Write();
    cut_phi->Write();
    cut_eta->Write();
    cut_leadPt->Write();
    cut_consPt->Write();
    
    effic_px->Write();
    effic_py->Write();
    effic_pz->Write();
    effic_E->Write();
    effic_phi->Write();
    effic_eta->Write();
    effic_leadPt->Write();
    effic_consPt->Write();
    //orig->Write();
  
  return 0;
}
