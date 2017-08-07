//  comparemethods.C
//  Created by Isaac Mooney on 7/12/17.

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

const string path = "out/";
const string histos = "compare/";

//const string type = "leadingJetpTs";
//const string branchtype = "c_effic_leadPt";//type;

void comparemethods () {
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    TFile *range2022 = TFile::Open((path + "trackeffic2022.root"/*"toy.root"*/).c_str());
    TFile *range2232 = TFile::Open((path + "trackeffic2232.root").c_str());
    //TFile *myfile = TFile::Open((path + "trackeffic.root").c_str());
    
    TCanvas * c = new TCanvas("c", "canvas", 800, 800);
    //TH1D *corrected = (TH1D*)myfile->Get("ptcorr");
    //TH1D *partlevel = (TH1D*)myfile->Get("ptpart");
    TH1D *corrected2022 = (TH1D*)range2022->Get("ptafter");
    TH1D *partlevel2022 = (TH1D*)range2022->Get("ccut2pt");
    TH1D *corrected2232 = (TH1D*)range2232->Get("ptafter");
    TH1D *partlevel2232 = (TH1D*)range2232->Get("ccut2pt");
    
    corrected2022->Add(corrected2232);
    partlevel2022->Add(partlevel2232);
    
    TPad * pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->SetGridx();
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    corrected2022->SetStats(0);
    
    //normalization
    Int_t nEntriesS = corrected2022->GetEntries();
    Int_t nEntriesE = partlevel2022->GetEntries();
    corrected2022->Scale(1/(double) nEntriesS);
    partlevel2022->Scale(1/(double) nEntriesE);
    
    corrected2022->Draw();
    partlevel2022->Draw("same");
    
    TLegend * test = new TLegend(0.6,0.75,0.8,0.875);
    TLegendEntry *sleg = test->AddEntry("ptafter", "Corrected", "lep");
    TLegendEntry *eleg = test->AddEntry("ccut2pt", "Particle Level", "lep");
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
    
    TH1D *ratio = (TH1D*)corrected2022->Clone("ratio");
    ratio->SetLineColor(kBlack);
    //ratio->SetMinimum(-5);
    //ratio->SetMaximum(5);
    //ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(partlevel2022);
    ratio->SetMarkerStyle(21);
    ratio->Draw("ep");
    
    corrected2022->SetLineColor(kBlue+1);
    corrected2022->SetLineWidth(2);
    
    corrected2022->GetYaxis()->SetTitleSize(20);
    corrected2022->GetYaxis()->SetTitleFont(43);
    corrected2022->GetYaxis()->SetTitleOffset(1.55);
    
    partlevel2022->SetLineColor(kRed);
    partlevel2022->SetLineWidth(2);
    
    ratio->SetTitle("");
    
    ratio->GetYaxis()->SetTitle("Corrected / Particle Level");
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
    
    TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);
    c1->cd();
    
    TH2D *numlostvdetpt2232 = (TH2D*)range2232->Get("numlostvdetpt");
    TH2D *numlostvdetpt2022 = (TH2D*)range2022->Get("numlostvdetpt");
                                                    
    numlostvdetpt2022->Add(numlostvdetpt2232);
    
    gPad->SetLogy(0);
    gPad->SetLogz();
    
    numlostvdetpt2022->Draw("COLZ");
    
    TCanvas * c2 = new TCanvas("c2", "c2", 800, 800);
    c2->cd();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    
    TH1D * proj1  = numlostvdetpt2022->ProjectionY("proj1", 0, 33);
    TH1D * proj2  = numlostvdetpt2022->ProjectionY("proj2", 34, 66);
    TH1D * proj3  = numlostvdetpt2022->ProjectionY("proj3", 67, 100);
    
    proj1->GetXaxis()->SetTitle("Number loss");
    proj1->SetTitle("");
    
    proj1->Scale(1/proj1->GetEntries());
    proj2->Scale(1/proj2->GetEntries());
    proj3->Scale(1/proj3->GetEntries());
    
    TF1 *bin = new TF1("bin","[0]*ROOT::Math::binomial_pdf(x,[1],[2])",0,20);
    bin->SetParameter(0,1);
    bin->SetParameter(1,0.2);
    bin->SetParameter(2,5);
    proj1->Fit("bin","ME");
    
    TF1 *nb = new TF1("nb","[0]*ROOT::Math::negative_binomial_pdf(x,[1],[2])",0,20);
    nb->SetParameter(0,1);
    nb->SetParameter(1,0.8);
    nb->SetParameter(2,4);
    proj1->Fit("nb","MEI");
    
    TF1 *pois = new TF1("pois","[0]*ROOT::Math::poisson_pdf(x,[1])",0,20);
    pois->SetParameter(0,1);
    pois->SetParameter(1,2);
    proj1->Fit("pois","MEI");
    
    
    proj1->SetLineColor(1);
    proj2->SetLineColor(2);
    proj3->SetLineColor(3);
    
    bin->SetLineColor(1);
    nb->SetLineColor(2);
    pois->SetLineColor(3);
    
    proj1->DrawCopy();
    proj2->DrawCopy("same");
    proj3->DrawCopy("same");
    
    bin->DrawCopy("same");
    nb->DrawCopy("same");
    pois->DrawCopy("same");

    TLegend * leg = new TLegend(0.6,0.75,0.8,0.875);
    leg->AddEntry("proj1", "  0 < p_{T}^{det} < 10 [GeV]", "lep");
    leg->AddEntry("proj2", "10 < p_{T}^{det} < 20", "lep");
    leg->AddEntry("proj3", "20 < p_{T}^{det} < 30", "lep");
    leg->AddEntry("bin","Binomial","l");
    leg->AddEntry("nb","Negative Binomial","l");
    leg->AddEntry("pois","Poisson","l");
    
    leg->SetBorderSize(0);
    leg->Draw("same");
    
    return;
}
