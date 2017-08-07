//  check.C

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include <iostream>
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;

//WE'RE PRETENDING THESE ARE ALL CHARGED

void check() {
    
    TFile *dists = new TFile("dists.root", "RECREATE");
    //TFile *ptsamp = TFile::Open("out/trackeffic.root", "READ");
    
    //TH1 *ptdist     = new TH1D("pt", "pt", 200, 0, 20);
    TH1 *numdist    = new TH1D("num", "num", 30, 0, 30);
    //TH1 *conspt = (TH1D*)ptsamp->Get("cons_pt");
    
    //Double_t pt;
    Int_t num;
    //TRandom3* pt_rand   = new TRandom3();
    TRandom3* num_rand  = new TRandom3();
    //pt_rand->SetSeed(10);
    num_rand->SetSeed(7);
    
    for (unsigned i = 0; i < 1e6; ++ i) {
        //pt = pt_rand->Exp(1);
        //pt = conspt->GetRandom();
        //num = 15;
        num = num_rand->Poisson(5);
        //ptdist->Fill(pt);
        numdist->Fill(num);
    }
    
    //delete pt_rand; pt_rand = NULL;
    delete num_rand; num_rand = NULL;
    
    dists->Write();
    dists->Close();
    
    return;
}
