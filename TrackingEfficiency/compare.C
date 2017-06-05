//  compare.C
//  Created by Isaac Mooney on 5/23/17.

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

const string path8 = "out/";
const string path6 = "p6out/";
const string histos = "compare/";

const string type = "leadingJetpTs";
const string branchtype = "c_effic_leadPt";//type;

void compare () {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    TCanvas * c = new TCanvas("c", "canvas", 800, 800);
    TH1D * pTraw8 = new TH1D("h1", "p_{T} Pythia 8 v. Pythia 6", 100, 0, 40);
    TH1D * pTraw6 = new TH1D("h2", "", 100, 0, 40);
    
    TPad * pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->SetGridx();
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    pTraw8->SetStats(0);
    
    TFile *p8 = TFile::Open((path8 + "trackeffic.root").c_str());
    TFile *p6 = TFile::Open((path6 + "pythia6trackeffic.root").c_str());
    
    TTree* p8tree = (TTree*) p8->Get(type.c_str());
    TTree* p6tree = (TTree*) p6->Get(type.c_str());
    
    Double_t p8_Pt, p6_Pt;
    p8tree->SetBranchAddress((branchtype /*+ "Pt"*/).c_str(), &p8_Pt);
    p6tree->SetBranchAddress((branchtype /*+ "Pt"*/).c_str(), &p6_Pt);
    
    for(unsigned i = 0; i < p8tree->GetEntries(); ++ i) {
        p8tree->GetEntry(i);
        pTraw8->Fill(p8_Pt);
    }
    for(unsigned i = 0; i < p6tree->GetEntries(); ++ i) {
        p6tree->GetEntry(i);
        pTraw6->Fill(p6_Pt);
    }
    
    //normalization
    Int_t nEntries8 = pTraw8->GetEntries();
    Int_t nEntries6 = pTraw6->GetEntries();
    pTraw8->Scale(1/(double) nEntries8);
    pTraw6->Scale(1/(double) nEntries6);
    
    pTraw8->Draw();
    pTraw6->Draw("same");
    
    //MAKE TLEGEND!!!
    TLegend * test = new TLegend(0.5,0.75,0.85,0.875);
    TLegendEntry *p8leg = test->AddEntry("h1", "Pythia 8", "lep");
    TLegendEntry *p6leg = test->AddEntry("h2", "Pythia 6", "lep");
    test->SetBorderSize(0);
    test->Draw("same");
    
    //pTraw8->GetYaxis()->SetLabelSize(0.);
    TGaxis *axis = new TGaxis(-5, 20, -5, 220, 20, 220, 510, "");
    axis->SetLabelFont(43);
    axis->SetLabelSize(15);
    axis->Draw();
    
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridx();
    pad2->Draw();
    pad2->cd();
    
    TH1D *ratio = (TH1D*)pTraw8->Clone("ratio");
    ratio->SetLineColor(kBlack);
    //ratio->SetMinimum(0.8);
    ratio->SetMaximum(5);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(pTraw6);
    ratio->SetMarkerStyle(21);
    ratio->Draw("ep");
    
    pTraw8->SetLineColor(kBlue+1);
    pTraw8->SetLineWidth(2);
    
    pTraw8->GetYaxis()->SetTitleSize(20);
    pTraw8->GetYaxis()->SetTitleFont(43);
    pTraw8->GetYaxis()->SetTitleOffset(1.55);
    
    pTraw6->SetLineColor(kRed);
    pTraw6->SetLineWidth(2);
    
    ratio->SetTitle("");
    
    ratio->GetYaxis()->SetTitle("Pythia 8 / Pythia 6");
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.55);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(15);
    
    ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ratio->GetXaxis()->SetTitleSize(18);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(4.);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(15);
    
    //c->SaveAs((histos + type + "_pT.pdf").c_str());
    
    return;
}
