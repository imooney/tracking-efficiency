//  test.C
//  Created by Isaac Mooney on 7/3/17.

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include <iostream>
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"

#include <vector>

using namespace std;

const unsigned num_events = 1e4;
const unsigned num_tracks = 30;

void test() {
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    TFile *samp = new TFile("dists.root", "READ");
    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    //gStyle->SetOptStat(0);
    TH1 *poisloss = new TH1D("tracks lost","tracks lost", 30, 0, 30);
    TH1 *flatloss = new TH1D("a","b", 30, 0, 30);
    TH1 *ptdist = (TH1D*)samp->Get("pt");
    TH1 *numdist = (TH1D*)samp->Get("num");
    TH1 *pt_loss = new TH1D("pt_loss","pt loss",30,0,30);
    //numdist->GetXaxis()->SetLimits(30,30.5);
    TF1 * p = new TF1("p", "ROOT::Math::poisson_pdf(x,6)", 0, 30);
    
    for (unsigned i = 0; i < num_events; ++ i) {
        const int size = numdist->GetRandom();
        unsigned sampled_num = 10000;
        Double_t pt_lost = 0;
        //cout << size << endl;
        vector<double> po(size, -1 ); vector<double> fl(30, -1 );
        unsigned trackslost = 0;
        
        for (unsigned j = 0; j < po.size(); ++ j) {
            double effic_num = gRandom->Uniform(0.0, 1.0);
            if (effic_num < 0.8) {
                po[j] = ptdist->GetRandom();
            }
        }
        
        for (unsigned j = 0; j < fl.size(); ++ j) {
            double effic_num = gRandom->Uniform(0.0, 1.0);
            if (effic_num < 0.8) {
                fl[j] = ptdist->GetRandom();
            }
        }
        
        
        for (unsigned j = 0; j < po.size(); ++ j) {
            if (po[j] == -1) {
                ++ trackslost;
            }
        }
        poisloss->Fill(trackslost);
        trackslost = 0;
        
        for (unsigned j = 0; j < fl.size(); ++ j) {
            if (fl[j] == -1) {
                ++ trackslost;
            }
        }
        flatloss->Fill(trackslost);
        
        sampled_num = p->GetRandom();
        for(unsigned j = 0; j < sampled_num; ++ j) {
            pt_lost += ptdist->GetRandom();
        }
        pt_loss->Fill(pt_lost);
    }
    poisloss->SetLineColor(4);
    flatloss->SetLineColor(12);
    
    poisloss->Scale(1/poisloss->Integral());
    flatloss->Scale(1/flatloss->Integral());
    
    /*
    TF1 * nb = new TF1("nb","[0]*ROOT::Math::negative_binomial_pdf(x,[1],[2])", 0, 30);
    nb->SetParameter(0, 1);
    nb->FixParameter(1, 0.8);
    nb->SetParameter(2, 24);
    nb->SetLineColor(1);
    
    TF1 * b = new TF1("b", "[0]*ROOT::Math::binomial_pdf(x,[1],[2])", 0, 30);
    b->SetParameter(0, 1);
    b->FixParameter(1, 0.2);
    b->SetParameter(2, 30);
    b->SetLineColor(2);
    
    TF1 * p = new TF1("p","[0]*ROOT::Math::poisson_pdf(x,[1])", 0, 30);
    p->SetParameter(0, 1);
    p->SetParameter(1, 6);
    p->SetLineColor(3);
    */
    //poisloss->Fit("nb", "MEI"); poisloss->Fit("b", "ME"); poisloss->Fit("p", "MEI");
    //flatloss->Fit("nb", "MEI"); flatloss->Fit("b", "ME"); flatloss->Fit("p", "MEI");
    
    c1->SetLogy();
    pt_loss->DrawCopy();
    //poisloss->DrawCopy(); /*nb->DrawCopy("same"); b->DrawCopy("same"); p->DrawCopy("same");*/ //flatloss->DrawCopy("same");
    //try increasing number of tracks to 30 or 40
    
    //TLegend * leg = new TLegend(0.5,0.75,0.75,0.875);
    /*TLegendEntry *pois =*/ //leg->AddEntry(poisloss, "Poisson sampling", "lep");
    /*TLegendEntry *flat =*/ //leg->AddEntry(flatloss, "Known track number", "lep");
    /*TLegendEntry *pleg = *///leg->AddEntry("p", "Poisson", "l");
    /*TLegendEntry *bleg = *///leg->AddEntry("b", "Binomial", "l");
    /*TLegendEntry *nbleg = *///leg->AddEntry("nb", "Negative binomial", "l");
    //leg->SetBorderSize(0);
    //leg->Draw("same");
    
    //samp->Close();
    //c1->SaveAs("testfullrange.pdf");
    
    return;
}
