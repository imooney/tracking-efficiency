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
  // Generator. Process selection. LHC initialization. Histogram.
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
    
  TH1D* orig = new TH1D("orig", "no efficiency correction", 100, 0, 5);
      
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
        if (pythia.event[i].isFinal() && abs(pythia.event[i].eta()) < 1) {
          //++nCharged;
          //if it's the 5th charged particle, exclude it (so we scale charged particles in histogram by 0.8)
          if (counter == 4 && pythia.event[i].isCharged()) {
            //reset the counter for later
            counter = 0;
            //fill the original jet, skip filling the tracking efficiency applied jet
            orig_part.push_back(PseudoJet(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
            continue;
          }
          else {
              //fill both
              orig_part.push_back(PseudoJet(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
              part.push_back(PseudoJet(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
          }
          //every charged particle increments the counter
          if (pythia.event[i].isCharged()) {
              counter ++;
          }
        }
      }
    
    // run the clustering, extract the jets
    ClusterSequence cs(part, jet_def);
    ClusterSequence c_s(orig_part, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    vector<PseudoJet> orig_jets = sorted_by_pt(c_s.inclusive_jets());

    //cout << "Clustering with " << jet_def.description() << '\n';

    //pseudojet axis cut
    if (abs(jets[0].eta()) < (1 - R)) {
      //leading jet pt spectrum
      pt = jets[0].pt();
      leading_pt->Fill(pt);
      tout->Fill();
    }
    //pseudojets with charged particle tracking efficiency selection
    for (int i = 0; i < jets.size(); ++ i) {
      //pseudojet axis cut
      if (abs(jets[i].eta()) < (1 - R)) {
	px = jets[i].px();                                                                                                                  
	py = jets[i].py();                                                                                                                  
	pz = jets[i].pz();                                                                                                                  
	E = jets[i].e();                                                                                                                    
	phi = jets[i].phi();                                                                                                                
	eta = jets[i].eta();
	con_pt = jets[i].pt();
	tout->Fill();                                                                                                                           
	pseudojet_px->Fill(jets[i].px());                                                                                                   
	pseudojet_py->Fill(jets[i].py());      
	pseudojet_pz->Fill(jets[i].pz());                                                                                                   
	pseudojet_E->Fill(jets[i].e());                                                                                                     
	pseudojet_phi->Fill(jets[i].phi());
	pseudojet_eta->Fill(jets[i].eta());
	constituent_pt->Fill(jets[i].pt());
      }
    }
    //pseudojets without charged particle tracking efficiency selection
    for (int i = 0; i < orig_jets.size(); ++ i) {
      //pseudojet axis cut
      if (abs(orig_jets[i].eta()) < (1 - R)) {
	orig->Fill(orig_jets[i].pt());
      }
    }

    //mult.fill( nCharged );
    
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  //cout << mult << '\n';
  
  cout << "Writing to " << fout->GetName() << '\n';
  //save histograms to output root file
  //tout->Print();
  tout->Write();
  pseudojet_px->Write();
  pseudojet_py->Write();
  pseudojet_pz->Write();
  pseudojet_E->Write();
  pseudojet_phi->Write();
  pseudojet_eta->Write();
  leading_pt->Write();
  orig->Write();
  constituent_pt->Write();
  
  return 0;
}
