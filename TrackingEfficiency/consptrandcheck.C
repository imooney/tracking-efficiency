//  consptrandcheck.C
//  Created by Isaac Mooney on 7/25/17.

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;

void consptrandcheck() {
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    TFile *samp = TFile::Open("out/trackeffic.root", "READ");
    
    TCanvas *c = new TCanvas("c", "canvas", 800, 800);
    //gStyle->SetOptStat(0);
    TH1 *conspt = (TH1D*)samp->Get("cons_pt");
    TH1 *ptrand = (TH1D*)samp->Get("ptofrands");
    
    TPad * pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->SetGridx();
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    conspt->SetStats(0);
    
    conspt->Scale(1/conspt->Integral());
    ptrand->Scale(1/ptrand->Integral());
    
    //conspt->Fit("expo");
    
    
    conspt->Draw();
    ptrand->Draw("same");
    
    TLegend * leg = new TLegend(0.6,0.75,0.8,0.875);
    leg->AddEntry("cons_pt", "Constituent p_{T}", "lep");
    leg->AddEntry("ptofrands", "Randomly sampled p_{T}", "lep");
    leg->SetBorderSize(0);
    leg->Draw("same");
    
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
    
    TH1D *ratio = (TH1D*)conspt->Clone("ratio");
    ratio->SetLineColor(kBlack);
    //ratio->SetMinimum(-5);
    //ratio->SetMaximum(5);
    //ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(ptrand);
    ratio->SetMarkerStyle(21);
    ratio->Draw("ep");
    
    conspt->SetLineColor(kBlue+1);
    conspt->SetLineWidth(2);
    
    conspt->GetYaxis()->SetTitleSize(20);
    conspt->GetYaxis()->SetTitleFont(43);
    conspt->GetYaxis()->SetTitleOffset(1.55);
    
    ptrand->SetLineColor(kRed);
    ptrand->SetLineWidth(2);
    
    ratio->SetTitle("");
    
    ratio->GetYaxis()->SetTitle("Constituent / Sampled");
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
    
    //samp->Close();
    
    return;
}
