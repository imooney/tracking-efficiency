//  fittoynumloss.C
//  Created by Isaac Mooney on 7/26/17.

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;

void fittoynumloss() {
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);

    c1->cd();
    
    TFile *myfile = TFile::Open("out/toy.root","READ");
    TH1 *hnumloss = (TH1D*)myfile->Get("numloss");
    TH2 *compare2D = (TH2D*)myfile->Get("compare2D");
    TH2 *multvaddedback = (TH2D*)myfile->Get("multvaddedback");
    gPad->SetLogy();
    
    hnumloss->Scale(1/hnumloss->Integral());
    
    TF1 *bin = new TF1("bin","[0]*ROOT::Math::binomial_pdf(x,[1],[2])",0,20);
    bin->SetParameter(0,1);
    bin->SetParameter(1,0.2);
    bin->SetParameter(2,5);
    hnumloss->Fit("bin","ME");
    
    hnumloss->DrawCopy(); bin->DrawCopy("same");
    
    c2->cd();
    gPad->SetLogy(0);
    gPad->SetLogz();
    
    compare2D->GetXaxis()->SetTitle("Particle level p_{T} [GeV/c]");
    compare2D->GetYaxis()->SetTitle("Corrected detector level p_{T} [GeV/c]");
    
    compare2D->DrawCopy("COLZ");
    TF1 *line45 = new TF1("line45","x",0,50);
    line45->DrawCopy("same");
    
    c3->cd();
    gPad->SetLogz();
    
    multvaddedback->GetXaxis()->SetTitle("Detector level track multiplicity");
    multvaddedback->GetYaxis()->SetTitle("Number of tracks added back");
    
    multvaddedback->DrawCopy("COLZ");
    
    TH1D * projX01 = multvaddedback->ProjectionY("projX01", 1, 1);
    projX01->Scale(1/projX01->Integral());
    /*
    TH1D * projX12 = multvaddedback->ProjectionY("projX12", 2, 2);
    projX12->Scale(1/projX12->Integral());
    TH1D * projX23 = multvaddedback->ProjectionY("projX23", 3, 3);
    projX23->Scale(1/projX23->Integral());
    TH1D * projX34 = multvaddedback->ProjectionY("projX34", 4, 4);
    projX34->Scale(1/projX34->Integral());
    TH1D * projX45 = multvaddedback->ProjectionY("projX45", 5, 5);
    projX45->Scale(1/projX45->Integral());
    TH1D * projX56 = multvaddedback->ProjectionY("projX56", 6, 6);
    projX56->Scale(1/projX56->Integral());
    TH1D * projX67 = multvaddedback->ProjectionY("projX67", 7, 7);
    projX67->Scale(1/projX67->Integral());
     */
    TH1D * projX78 = multvaddedback->ProjectionY("projX78", 8, 8);
    projX78->Scale(1/projX78->Integral());
    /*
    TH1D * projX89 = multvaddedback->ProjectionY("projX89", 9, 9);
    projX89->Scale(1/projX89->Integral());
    TH1D * projX910 = multvaddedback->ProjectionY("projX910", 10, 10);
    projX910->Scale(1/projX910->Integral());
    TH1D * projX1011 = multvaddedback->ProjectionY("projX1011", 11, 11);
    projX1011->Scale(1/projX1011->Integral());
    TH1D * projX1112 = multvaddedback->ProjectionY("projX1112", 12, 12);
    projX1112->Scale(1/projX1112->Integral());
    TH1D * projX1213 = multvaddedback->ProjectionY("projX1213", 13, 13);
    projX1213->Scale(1/projX1213->Integral());
     */
    TH1D * projX1314 = multvaddedback->ProjectionY("projX1314", 14, 14);
    projX1314->Scale(1/projX1314->Integral());
    /*
    TH1D * projX1415 = multvaddedback->ProjectionY("projX1415", 15, 15);
    projX1415->Scale(1/projX1415->Integral());
    TH1D * projX1516 = multvaddedback->ProjectionY("projX1516", 16, 16);
    projX1516->Scale(1/projX1516->Integral());
    */
    
    projX01->SetLineColor(1);
    /*
    projX12->SetLineColor(2);
    projX23->SetLineColor(3);
    projX34->SetLineColor(4);
    projX45->SetLineColor(5);
    projX56->SetLineColor(6);
    projX67->SetLineColor(7);
     */
    projX78->SetLineColor(2);
    /*
    projX89->SetLineColor(9);
    projX910->SetLineColor(10);
    projX1011->SetLineColor(11);
    projX1112->SetLineColor(12);
    projX1213->SetLineColor(13);
     */
    projX1314->SetLineColor(3);
     /*
    projX1415->SetLineColor(15);
    projX1516->SetLineColor(16);
    */
    c4->cd();
    projX01->DrawCopy();
    /*
    projX12->DrawCopy("same");
    projX23->DrawCopy("same");
    projX34->DrawCopy("same");
    projX45->DrawCopy("same");
    projX56->DrawCopy("same");
    projX67->DrawCopy("same");
     */
    projX78->DrawCopy("same");
    /*
    projX89->DrawCopy("same");
    projX910->DrawCopy("same");
    projX1011->DrawCopy("same");
    projX1112->DrawCopy("same");
    projX1213->DrawCopy("same");
     */
    projX1314->DrawCopy("same");
     /*
    projX1415->DrawCopy("same");
    projX1516->DrawCopy("same");
    */
    
    TLegend * test = new TLegend(0.6,0.75,0.8,0.875);
    test->AddEntry("projX01", "multiplicity = 1", "lep");
    test->AddEntry("projX78", "multiplicity = 8", "lep");
    test->AddEntry("projX1314", "multiplicity = 14", "lep");
    
    test->SetBorderSize(0);
    test->Draw("same");
    
    return;
}
