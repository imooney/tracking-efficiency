#define DiffClass_cxx
#include "DiffClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void DiffClass::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L DiffClass.C
//      Root > DiffClass t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0) return;
    
    string path = "../histos/diffs/";
    string type = "num_diff";
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    
    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    
    TH1 * numdiffraw = new TH1D((type +"_raw").c_str(), "Number difference after efficiency and raw", 100, 0, 60);
    TH1 * numdiff = new TH1D((type).c_str(), "Number difference after efficiency and before", 100, 0, 16);
    TH1 * cnumdiffraw = new TH1D(("c_" + type +"_raw").c_str(), "Number difference after efficiency and raw (charged)", 100, 0, 60);
    TH1 * cnumdiff = new TH1D(("c_" + type).c_str(), "Number difference after efficiency and before (charged)", 100, 0, 16);
    TH1 * pt_diff = new TH1D("ptdiff", "p_{T} difference between lead jet before/after effiency", 100, -21, 13);
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        numdiffraw->Fill(num_diff_raw);
        cnumdiffraw->Fill(c_num_diff_raw);
        numdiff->Fill(num_diff);
        cnumdiff->Fill(c_num_diff);
        pt_diff->Fill(ptdiff);
    }
    Int_t nEntries = pt_diff->GetEntries();
    
    pt_diff->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    pt_diff->Scale(1/(double) nEntries); c1->SetLogy(); pt_diff->DrawCopy(); c1->SaveAs((path + "ptdiff.pdf").c_str());
    gPad->Modified(); gPad->Update(); c1->SetLogy(0);
    
    numdiffraw->GetXaxis()->SetTitle("Number of jets");
    numdiffraw->Scale(1/(double) nEntries); c1->SetLogy(); numdiffraw->DrawCopy(); c1->SaveAs((path + type + "_raw.pdf").c_str());
    gPad->Modified(); gPad->Update(); c1->SetLogy(0);
    
    numdiff->GetXaxis()->SetTitle("Number of jets");
    numdiff->Scale(1/(double) nEntries); c1->SetLogy(); numdiff->DrawCopy(); c1->SaveAs((path + type + ".pdf").c_str());
    gPad->Modified(); gPad->Update(); c1->SetLogy(0);
    
    cnumdiffraw->GetXaxis()->SetTitle("Number of jets");
    cnumdiffraw->Scale(1/(double) nEntries); c1->SetLogy(); cnumdiffraw->DrawCopy(); c1->SaveAs((path + "c_" + type + "_raw.pdf").c_str());
    gPad->Modified(); gPad->Update(); c1->SetLogy(0);
    
    cnumdiff->GetXaxis()->SetTitle("Number of jets");
    cnumdiff->Scale(1/(double) nEntries); c1->SetLogy(); cnumdiff->DrawCopy(); c1->SaveAs((path + "c_" + type + ".pdf").c_str());
    gPad->Modified(); gPad->Update(); c1->SetLogy(0);
    
}
