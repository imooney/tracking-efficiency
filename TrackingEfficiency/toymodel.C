//  toymodel.C
//  Created by Isaac Mooney on 7/26/17.

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include <iostream>
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"

//#include <string>
#include <vector>

using namespace std;

const unsigned num_events = 1e5;
//const unsigned num_tracks = 5;
const double efficiency = 0.8;

//we pretend only charged tracks in the following

void toymodel() {
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    TCanvas * c1 = new TCanvas ("c1","c1",800,600);
    
    TFile *ptsamp = TFile::Open("out/trackeffic.root", "READ"); // Open input files
    TFile *dists = TFile::Open("dists.root", "READ");
    TFile *fout = new TFile("out/toy.root","RECREATE");              // Create output file
    TH1 *conspt = (TH1D*)ptsamp->Get("cons_pt");                // Read in histograms
    TH1 *numdist = (TH1I*)dists->Get("num");
    TH1 *hjetptparticle = new TH1D("ptpart","",500,0,50);
    TH1 *hjetptloss = new TH1D("ptloss","",500,0,50);
    TH1 *hjetptdetector = new TH1D("ptdet","",500,0,50);
    TH1 *hjetptcorrected = new TH1D("ptcorr","",500,0,50);
    TH1 *hnumparticle = new TH1D("numpart","",20,0,20);
    TH1 *hnumloss = new TH1D("numloss","",20,0,20);
    TH1 *hnumdetector = new TH1D("numdet","",20,0,20);
    TH1 *hnumcorrected = new TH1D("numcorr","",20,0,20);
    TH2 *compare2D = new TH2D("compare2D","",500,0,50,500,0,50);
    //TFile *fit2d = new TFile("out/fit2D.root", "RECREATE");
    //TH2 *multvaddedback = new TH2D("multvaddedback","",20,0,20,7,0,7);
    TFile *fit2d = TFile::Open("out/fit2D.root");
    TH2 *multvaddedback = (TH2D*)fit2d->Get("multvaddedback");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2D function for number loss
    
    TH1D * proj1  = multvaddedback->ProjectionY("proj1", 1, 1);
    TH1D * proj2  = multvaddedback->ProjectionY("proj2", 2, 2);
    TH1D * proj3  = multvaddedback->ProjectionY("proj3", 3, 3);
    TH1D * proj4  = multvaddedback->ProjectionY("proj4", 4, 4);
    TH1D * proj5  = multvaddedback->ProjectionY("proj5", 5, 5);
    TH1D * proj6  = multvaddedback->ProjectionY("proj6", 6, 6);
    TH1D * proj7  = multvaddedback->ProjectionY("proj7", 7, 7);
    TH1D * proj8  = multvaddedback->ProjectionY("proj8", 8, 8);
    TH1D * proj9  = multvaddedback->ProjectionY("proj9", 9, 9);
    TH1D * proj10 = multvaddedback->ProjectionY("proj10", 10, 10);
    TH1D * proj11 = multvaddedback->ProjectionY("proj11", 11, 11);
    TH1D * proj12 = multvaddedback->ProjectionY("proj12", 12, 12);
    TH1D * proj13 = multvaddedback->ProjectionY("proj13", 13, 13);
    TH1D * proj14 = multvaddedback->ProjectionY("proj14", 14, 14);
    TH1D * proj15 = multvaddedback->ProjectionY("proj15", 15, 15);
    TH1D * proj16 = multvaddedback->ProjectionY("proj16", 16, 16);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    for(unsigned j = 0; j < num_events; ++ j) {
    //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\ <- one event
    //THIS CAN BE CHANGED BACK TO A GLOBAL NUM_TRACKS ABOVE
        unsigned num_tracks = numdist->GetRandom();
    //initialize a fixed number of tracks (particle level) and 0 tracks (detector level)
    vector<Double_t> tracksparticle(num_tracks, -1);
    vector<Double_t> tracksdetector[0] = {};
    //initialize a number loss and ptloss vector
    unsigned num_loss = 0;
    vector<Double_t> ptloss[0] = {};
        
    tracksdetector->reserve(15);
    ptloss->reserve(15);
    
    //initialize a jet pt which we'll add to with each track pt in the following loop
    Double_t jetptparticle = 0;
    //loop over the number of tracks giving each a random pt from the constituent spectrum
    for (unsigned i = 0; i < tracksparticle.size(); ++ i) {
        tracksparticle[i] = conspt->GetRandom();
        //conspt->GetRandom();
        jetptparticle += tracksparticle[i];
    }
    
    //apply efficiency
    for (unsigned i = 0; i < tracksparticle.size(); ++ i) {
        double effic_num = gRandom->Uniform(0.0, 1.0);
        //if true, get rid of track
        if (effic_num > efficiency) {
            ++ num_loss;
            ptloss->push_back(tracksparticle[i]);
        }
        //else, don't get rid of it
        else {
            tracksdetector->push_back(tracksparticle[i]);
        }
    }
    
    //CHECK:
    if (num_loss != tracksparticle.size() - tracksdetector->size()) {
        cerr << "YA DONE MESSED UP A-ARON!\n";
        exit(1);
    }
        
        unsigned thrown;
        char name[16];
        sprintf(name, "%s%d","proj",tracksdetector->size());
        TH1D *h = (TH1D*)gDirectory->Get(name); // THIS IS THE MISSING LINE
        
        if ( h ) { // make sure the Get succeeded
            //if (h->GetIntegral() == 0) cout << "yup duh" << '\n';
            if (h->GetEntries() > 0) {
              thrown = (unsigned)h->GetRandom();
            }
            else {
                cout << "whoops\n";
                cout << name << '\n';
                thrown = 1;
            }
         }
        else {/*cout << "get failed: " << name << " " << tracksdetector->size() << '\n';*/thrown = 0;}
        
        
    //TF1 * b = new TF1("binomial","ROOT::Math::binomial_pdf(x,0.2,5)",0,10);
    //unsigned thrown = b->GetRandom();
    Double_t thrownpt = 0;
    //add up a total jet pt loss from n (given by the binomial) track pts
    //sampled from constituent spectrum
    for (unsigned i = 0; i < thrown; ++ i) {
        thrownpt += conspt->GetRandom();
    }
    
    //first calculate a jet pt for detector level jets by looping over tracks:
    Double_t jetptdetector = 0;
    for (unsigned i = 0; i < tracksdetector->size(); ++ i) {
        jetptdetector += tracksdetector->at(i);
    }
    Double_t jetptcorrected = jetptdetector + thrownpt;
    
    //fill histograms:
    hjetptparticle->Fill(jetptparticle); hjetptloss->Fill(thrownpt);
    hjetptdetector->Fill(jetptdetector); hjetptcorrected->Fill(jetptcorrected);
    hnumparticle->Fill(tracksparticle.size()); hnumloss->Fill(thrown);
    hnumdetector->Fill(tracksdetector->size()); hnumcorrected->Fill(tracksdetector->size() + thrown);
    compare2D->Fill(jetptparticle, jetptcorrected);
    //multvaddedback->Fill(tracksdetector->size(), tracksparticle.size() - tracksdetector->size());//thrown);
    
    //remember to clear all necessary things
        tracksparticle.clear(); tracksdetector->clear(); ptloss->clear();
    //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
        
        if (j % 1000 == 0) { cout << j << " events\n";}
    }
    //at the end we want to write the distributions to a file
    cout << hjetptparticle->GetMeanError() << '\n'; cout << hjetptcorrected->GetMeanError() << '\n';
    
    fout->Write();//fit2d->Write();
    
    return;
}
