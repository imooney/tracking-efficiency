//  compareConstituents.C
//  Created by Isaac Mooney on 6/12/17.

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
//const string path6 = "p6out/";
const string histos = "compareCons/";

//const string corntree   = "charged";
const string tree1      = "chargedEfficConstituents";
const string tree2      = "chargedCut2Constituents";
const string tree3      = "chargedCut2";

//const string corn       = "c_";
const string branch1    = "c_effic_cons_Pt";
const string branch2    = "c_cut2_cons_Pt";

void compareConstituents () {
    //gStyle->SetOptStat(1);
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    TCanvas * c = new TCanvas("c", "canvas", 800, 800);
    TH1D * effic_pt = new TH1D("h1", "constituent p_{T} w/ v. w/o efficiency (charged)", 100, 0, 31);
    TH1D * cut2_pt = new TH1D("h2", "", 100, 0, 31);
    
    TPad * pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->SetGridx();
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    effic_pt->SetStats(0);
    
    TFile *myfile = TFile::Open((path8 + "trackeffic.root").c_str());
    
    TTree* efficTree    = (TTree*) myfile->Get(tree1.c_str());
    TTree* cut2Tree     = (TTree*) myfile->Get(tree2.c_str());
    TTree* normTree     = (TTree*) myfile->Get(tree3.c_str());
    
    //GET ENTRIES IN NORMTREE LATER
    
    Double_t effic_Pt, cut2_Pt;
    efficTree->SetBranchAddress((branch1).c_str(), &effic_Pt);
    cut2Tree->SetBranchAddress((branch2).c_str(), &cut2_Pt);
    
    
    for(unsigned i = 0; i < efficTree->GetEntries(); ++ i) {
        efficTree->GetEntry(i);
        effic_pt->Fill(effic_Pt);
    }
    for(unsigned i = 0; i < cut2Tree->GetEntries(); ++ i) {
        cut2Tree->GetEntry(i);
        cut2_pt->Fill(cut2_Pt);
    }
    
    //normalization
    Int_t nEntriesCut2 = cut2_pt->GetEntries();
    Int_t nEntriesEffic = nEntriesCut2;//effic_pt->GetEntries();
    effic_pt->Scale(1/(double) nEntriesEffic);
    cut2_pt->Scale(1/(double) nEntriesCut2);
    
    effic_pt->Draw();
    cut2_pt->Draw("same");
    
    //MAKE TLEGEND!!!
    TLegend * test = new TLegend(0.5,0.75,0.85,0.875);
    TLegendEntry *efficleg = test->AddEntry("h1", "After efficiency", "lep");
    TLegendEntry *cut2leg = test->AddEntry("h2", "Before efficiency", "lep");
    test->SetBorderSize(0);
    test->Draw("same");
    
    //effic_pt->GetYaxis()->SetLabelSize(0.);
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
    
    TH1D *ratio = (TH1D*)effic_pt->Clone("ratio");
    ratio->SetLineColor(kBlack);
    //ratio->SetMinimum(0.8);
    ratio->SetMaximum(3);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(cut2_pt);
    ratio->SetMarkerStyle(21);
    ratio->Draw("ep");
    
    effic_pt->SetLineColor(kBlue+1);
    effic_pt->SetLineWidth(2);
    
    effic_pt->GetYaxis()->SetTitleSize(20);
    effic_pt->GetYaxis()->SetTitleFont(43);
    effic_pt->GetYaxis()->SetTitleOffset(1.55);
    
    cut2_pt->SetLineColor(kRed);
    cut2_pt->SetLineWidth(2);
    
    ratio->SetTitle("");
    
    ratio->GetYaxis()->SetTitle("with efficiency / without");
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
    
    c->SaveAs((histos + "c_compareCons_pT.pdf").c_str());
    
    return;
}
