#include <iostream>
#include <sstream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"

using namespace std;

void plot () {
  auto c0 = new TCanvas("a","b",0,0,1200,700);
  c0->cd();
  /*
  gROOT->SetStyle("Plain");                       // plain histogram style
  gStyle->SetOptStat("nemruoi");                  // expanded stats box
  gStyle->SetPadTickX(1);                         // tic marks on all axes
  gStyle->SetPadTickY(1);                         //
  gStyle->SetOptFit(1111);                        // the results of the fits
  gStyle->SetPadGridX(kTRUE);                     // draw horizontal and vertical grids
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetPalette(1,0);                        // pretty and useful palette
  gStyle->SetHistLineWidth(2);                    // a thicker histogram line
  gStyle->SetFrameFillColor(10);                  // a different frame colour   
  gStyle->SetTitleFillColor(33);                 // title colour to highlight it
  gStyle->SetTitleW(.77);                         // Title Width
  gStyle->SetTitleH(.07);                        // Title height
  gStyle->SetHistMinimumZero(true);               // Suppresses Histogram Zero Supression
  */
  gPad->SetLogy();
  //open file
  TFile myfile("trackeffic.root");
  myfile->ls();
  //get histograms

  //without cuts (jet level)                                                                                                                                                                                     
  TH1D* uncut_px; TH1D* uncut_py; TH1D* uncut_pz; TH1D* uncut_E;
  TH1D* uncut_phi;  TH1D* uncut_eta; TH1D* uncut_leadPt; TH1D* uncut_Pt;

  myfile.GetObject("u_px", uncut_px);
  myfile.GetObject("u_py", uncut_py);
  myfile.GetObject("u_pz", uncut_pz);
  myfile.GetObject("u_E", uncut_E);
  myfile.GetObject("u_phi", uncut_phi);
  myfile.GetObject("u_eta", uncut_eta);
  myfile.GetObject("u_leadPt", uncut_leadPt);
  myfile.GetObject("u_Pt", uncut_Pt);

  uncut_px->SetDirectory(0); uncut_py->SetDirectory(0); uncut_pz->SetDirectory(0); uncut_E->SetDirectory(0);
  uncut_phi->SetDirectory(0); uncut_eta->SetDirectory(0); uncut_leadPt->SetDirectory(0); uncut_Pt->SetDirectory(0);

  //charged, without cuts (jet level)

  TH1D* c_uncut_px; TH1D* c_uncut_py; TH1D* c_uncut_pz; TH1D* c_uncut_E;
  TH1D* c_uncut_phi;  TH1D* c_uncut_eta; TH1D* c_uncut_leadPt; TH1D* c_uncut_Pt;

  myfile.GetObject("c_u_px", c_uncut_px);
  myfile.GetObject("c_u_py", c_uncut_py);
  myfile.GetObject("c_u_pz", c_uncut_pz);
  myfile.GetObject("c_u_E", c_uncut_E);
  myfile.GetObject("c_u_phi", c_uncut_phi);
  myfile.GetObject("c_u_eta", c_uncut_eta);
  myfile.GetObject("c_u_leadPt", c_uncut_leadPt);
  myfile.GetObject("c_u_Pt", c_uncut_Pt);

  c_uncut_px->SetDirectory(0); c_uncut_py->SetDirectory(0); c_uncut_pz->SetDirectory(0); c_uncut_E->SetDirectory(0);
  c_uncut_phi->SetDirectory(0); c_uncut_eta->SetDirectory(0); c_uncut_leadPt->SetDirectory(0); c_uncut_Pt->SetDirectory(0);

  //without cuts (constituent level)                                                                                                                                                                               
  TH1D* uncut_cons_px; TH1D* uncut_cons_py; TH1D* uncut_cons_pz; TH1D* uncut_cons_E;
  TH1D* uncut_cons_phi;  TH1D* uncut_cons_eta; TH1D* uncut_cons_leadPt; TH1D* uncut_cons_Pt;

  myfile.GetObject("u_cons_px", uncut_cons_px);
  myfile.GetObject("u_cons_py", uncut_cons_py);
  myfile.GetObject("u_cons_pz", uncut_cons_pz);
  myfile.GetObject("u_cons_E", uncut_cons_E);
  myfile.GetObject("u_cons_phi", uncut_cons_phi);
  myfile.GetObject("u_cons_eta", uncut_cons_eta);
  myfile.GetObject("u_cons_Pt", uncut_cons_Pt);

  uncut_cons_px->SetDirectory(0); uncut_cons_py->SetDirectory(0); uncut_cons_pz->SetDirectory(0); uncut_cons_E->SetDirectory(0);
  uncut_cons_phi->SetDirectory(0); uncut_cons_eta->SetDirectory(0); uncut_cons_Pt->SetDirectory(0);

  //charged, without cuts (constituent level)                                                                                                                                                                      
  TH1D* c_uncut_cons_px; TH1D* c_uncut_cons_py; TH1D* c_uncut_cons_pz; TH1D* c_uncut_cons_E;
  TH1D* c_uncut_cons_phi;  TH1D* c_uncut_cons_eta; TH1D* c_uncut_cons_leadPt; TH1D* c_uncut_cons_Pt;

  myfile.GetObject("c_u_cons_px", c_uncut_cons_px);
  myfile.GetObject("c_u_cons_py", c_uncut_cons_py);
  myfile.GetObject("c_u_cons_pz", c_uncut_cons_pz);
  myfile.GetObject("c_u_cons_E", c_uncut_cons_E);
  myfile.GetObject("c_u_cons_phi", c_uncut_cons_phi);
  myfile.GetObject("c_u_cons_eta", c_uncut_cons_eta);
  myfile.GetObject("c_u_cons_leadPt", c_uncut_cons_leadPt);
  myfile.GetObject("c_u_cons_Pt", c_uncut_cons_Pt);

  c_uncut_cons_px->SetDirectory(0); c_uncut_cons_py->SetDirectory(0); c_uncut_cons_pz->SetDirectory(0); c_uncut_cons_E->SetDirectory(0);
  c_uncut_cons_phi->SetDirectory(0); c_uncut_cons_eta->SetDirectory(0); c_uncut_cons_Pt->SetDirectory(0);

  //with particle cuts                                                                     

  TH1D* cut_px; TH1D* cut_py; TH1D* cut_pz; TH1D* cut_E;
  TH1D* cut_phi;  TH1D* cut_eta; TH1D* cut_leadPt; TH1D* cut_Pt;

  myfile.GetObject("cut_px", cut_px);
  myfile.GetObject("cut_py", cut_py);
  myfile.GetObject("cut_pz", cut_pz);
  myfile.GetObject("cut_E", cut_E);
  myfile.GetObject("cut_phi", cut_phi);
  myfile.GetObject("cut_eta", cut_eta);
  myfile.GetObject("cut_leadPt", cut_leadPt);
  myfile.GetObject("cut_Pt", cut_Pt);

  cut_px->SetDirectory(0); cut_py->SetDirectory(0); cut_pz->SetDirectory(0); cut_E->SetDirectory(0);
  cut_phi->SetDirectory(0); cut_eta->SetDirectory(0); cut_leadPt->SetDirectory(0); cut_Pt->SetDirectory(0);

  //with particle cuts (charged)                                                                                                          

  TH1D* c_cut_px; TH1D* c_cut_py; TH1D* c_cut_pz; TH1D* c_cut_E;
  TH1D* c_cut_phi;  TH1D* c_cut_eta; TH1D* c_cut_leadPt; TH1D* c_cut_Pt;

  myfile.GetObject("c_cut_px", c_cut_px);
  myfile.GetObject("c_cut_py", c_cut_py);
  myfile.GetObject("c_cut_pz", c_cut_pz);
  myfile.GetObject("c_cut_E", c_cut_E);
  myfile.GetObject("c_cut_phi", c_cut_phi);
  myfile.GetObject("c_cut_eta", c_cut_eta);
  myfile.GetObject("c_cut_leadPt", c_cut_leadPt);
  myfile.GetObject("c_cut_Pt", c_cut_Pt);

  c_cut_px->SetDirectory(0); c_cut_py->SetDirectory(0); c_cut_pz->SetDirectory(0); c_cut_E->SetDirectory(0);
  c_cut_phi->SetDirectory(0); c_cut_eta->SetDirectory(0); c_cut_leadPt->SetDirectory(0); c_cut_Pt->SetDirectory(0);



  //with jet cuts                                                                   

  TH1D* cut2_px; TH1D* cut2_py; TH1D* cut2_pz; TH1D* cut2_E;
  TH1D* cut2_phi;  TH1D* cut2_eta; TH1D* cut2_leadPt; TH1D* cut2_Pt;

  myfile.GetObject("c2_px", cut2_px);
  myfile.GetObject("c2_py", cut2_py);
  myfile.GetObject("c2_pz", cut2_pz);
  myfile.GetObject("c2_E", cut2_E);
  myfile.GetObject("c2_phi", cut2_phi);
  myfile.GetObject("c2_eta", cut2_eta);
  myfile.GetObject("c2_leadPt", cut2_leadPt);
  myfile.GetObject("c2_Pt", cut2_Pt);

  cut2_px->SetDirectory(0); cut2_py->SetDirectory(0); cut2_pz->SetDirectory(0); cut2_E->SetDirectory(0);
  cut2_phi->SetDirectory(0); cut2_eta->SetDirectory(0); cut2_leadPt->SetDirectory(0); cut2_Pt->SetDirectory(0);

  //with jet cuts (charged)                                         

  TH1D* c_cut2_px; TH1D* c_cut2_py; TH1D* c_cut2_pz; TH1D* c_cut2_E;
  TH1D* c_cut2_phi;  TH1D* c_cut2_eta; TH1D* c_cut2_leadPt; TH1D* c_cut2_Pt;

  myfile.GetObject("c_c2_px", c_cut2_px);
  myfile.GetObject("c_c2_py", c_cut2_py);
  myfile.GetObject("c_c2_pz", c_cut2_pz);
  myfile.GetObject("c_c2_E", c_cut2_E);
  myfile.GetObject("c_c2_phi", c_cut2_phi);
  myfile.GetObject("c_c2_eta", c_cut2_eta);
  myfile.GetObject("c_c2_leadPt", c_cut2_leadPt);
  myfile.GetObject("c_c2_Pt", c_cut2_Pt);

  c_cut2_px->SetDirectory(0); c_cut2_py->SetDirectory(0); c_cut2_pz->SetDirectory(0); c_cut2_E->SetDirectory(0);
  c_cut2_phi->SetDirectory(0); c_cut2_eta->SetDirectory(0); c_cut2_leadPt->SetDirectory(0); c_cut2_Pt->SetDirectory(0);

  //with tracking efficiency selection applied on top of cuts - jet level   

  TH1D* effic_px; TH1D* effic_py; TH1D* effic_pz; TH1D* effic_E;
  TH1D* effic_phi;  TH1D* effic_eta; TH1D* effic_leadPt; TH1D* effic_Pt;
  TH1D* test_ptdiff; TH1D* test_without_leadsublead;

  myfile.GetObject("e_px", effic_px);
  myfile.GetObject("e_py", effic_py);
  myfile.GetObject("e_pz", effic_pz);
  myfile.GetObject("e_E", effic_E);
  myfile.GetObject("e_phi", effic_phi);
  myfile.GetObject("e_eta", effic_eta);
  myfile.GetObject("e_leadPt", effic_leadPt);
  myfile.GetObject("e_Pt", effic_Pt);
  myfile.GetObject("e_ptdiff", test_ptdiff);
  myfile.GetObject("e_leadsublead", test_without_leadsublead);

  effic_px->SetDirectory(0); effic_py->SetDirectory(0); effic_pz->SetDirectory(0); effic_E->SetDirectory(0);
  effic_phi->SetDirectory(0); effic_eta->SetDirectory(0); effic_leadPt->SetDirectory(0); effic_Pt->SetDirectory(0);
  test_ptdiff->SetDirectory(0); test_without_leadsublead->SetDirectory(0);

  //charged (denoted c_) fraction only, with tracking efficiency selection applied on top of cuts - jet level
  TH1D* c_effic_px; TH1D* c_effic_py; TH1D* c_effic_pz; TH1D* c_effic_E;
  TH1D* c_effic_phi;  TH1D* c_effic_eta; TH1D* c_effic_leadPt; TH1D* c_effic_Pt;

  myfile.GetObject("c_e_px", c_effic_px);
  myfile.GetObject("c_e_py", c_effic_py);
  myfile.GetObject("c_e_pz", c_effic_pz);
  myfile.GetObject("c_e_E", c_effic_E);
  myfile.GetObject("c_e_phi", c_effic_phi);
  myfile.GetObject("c_e_eta", c_effic_eta);
  myfile.GetObject("c_e_leadPt", c_effic_leadPt);
  myfile.GetObject("c_e_Pt", c_effic_Pt);

  c_effic_px->SetDirectory(0); c_effic_py->SetDirectory(0); c_effic_pz->SetDirectory(0); c_effic_E->SetDirectory(0);
  c_effic_phi->SetDirectory(0); c_effic_eta->SetDirectory(0); c_effic_leadPt->SetDirectory(0); c_effic_Pt->SetDirectory(0);

  TH1D* num_dif_raw;  TH1D* c_num_dif_raw; TH1D* num_dif; TH1D* c_num_dif;


  myfile.GetObject("num_diff_raw", num_dif_raw);
  myfile.GetObject("c_num_diff_raw", c_num_dif_raw);
  myfile.GetObject("num_diff", num_dif);
  myfile.GetObject("c_num_diff", c_num_dif);

  num_dif_raw->SetDirectory(0); c_num_dif_raw->SetDirectory(0); num_dif->SetDirectory(0); c_num_dif->SetDirectory(0);

  Int_t nJets = num_dif->GetEntries();
  num_dif->Scale(1/(double) nJets);

  TF1 *fnb = new TF1("fnb","[0]*ROOT::Math::negative_binomial_pdf(x,[1],[2])", 0, 18);
  fnb->SetLineColor(2);//fnb->DrawCopy("same");                                                                                                
  //  fnb->SetParameter(0,10);//10);
  //fnb->FixParameter(1,0.79);//0.77);
  //fnb->FixParameter(2,18.8);//9.6);

  num_dif->Fit(fnb,"QR+");

  /*
  uncut_px->Draw(); c0->SaveAs("histos/uncut_px.pdf"); c0->Clear();
  uncut_py->Draw(); c0->SaveAs("histos/uncut_py.pdf"); c0->Clear();
  uncut_pz->Draw(); c0->SaveAs("histos/uncut_pz.pdf"); c0->Clear();
  uncut_E->Draw(); c0->SaveAs("histos/uncut_E.pdf"); c0->Clear();
  uncut_phi->Draw(); c0->SaveAs("histos/uncut_phi.pdf"); c0->Clear();
  uncut_eta->Draw(); c0->SaveAs("histos/uncut_eta.pdf"); c0->Clear();
  uncut_leadPt->Draw(); c0->SaveAs("histos/uncut_leadPt.pdf"); c0->Clear();
  uncut_Pt->Draw(); c0->SaveAs("histos/uncut_Pt.pdf"); c0->Clear();

  c_uncut_px->Draw(); c0->SaveAs("histos/c_uncut_px.pdf"); c0->Clear();
  c_uncut_py->Draw(); c0->SaveAs("histos/c_uncut_py.pdf"); c0->Clear();
  c_uncut_pz->Draw(); c0->SaveAs("histos/c_uncut_pz.pdf"); c0->Clear();
  c_uncut_E->Draw(); c0->SaveAs("histos/c_uncut_E.pdf"); c0->Clear();
  c_uncut_phi->Draw(); c0->SaveAs("histos/c_uncut_phi.pdf"); c0->Clear();
  c_uncut_eta->Draw(); c0->SaveAs("histos/c_uncut_eta.pdf"); c0->Clear();
  c_uncut_leadPt->Draw(); c0->SaveAs("histos/c_uncut_leadPt.pdf"); c0->Clear();
  c_uncut_Pt->Draw(); c0->SaveAs("histos/c_uncut_Pt.pdf"); c0->Clear();

  uncut_cons_px->Draw(); c0->SaveAs("histos/uncut_cons_px.pdf"); c0->Clear();
  uncut_cons_py->Draw(); c0->SaveAs("histos/uncut_cons_py.pdf"); c0->Clear();
  uncut_cons_pz->Draw(); c0->SaveAs("histos/uncut_cons_pz.pdf"); c0->Clear();
  uncut_cons_E->Draw(); c0->SaveAs("histos/uncut_cons_E.pdf"); c0->Clear();
  uncut_cons_phi->Draw(); c0->SaveAs("histos/uncut_cons_phi.pdf"); c0->Clear();
  uncut_cons_eta->Draw(); c0->SaveAs("histos/uncut_cons_eta.pdf"); c0->Clear();
  uncut_cons_Pt->Draw(); c0->SaveAs("histos/uncut_cons_Pt.pdf"); c0->Clear();

  c_uncut_cons_px->Draw(); c0->SaveAs("histos/c_uncut_cons_px.pdf"); c0->Clear();
  c_uncut_cons_py->Draw(); c0->SaveAs("histos/c_uncut_cons_py.pdf"); c0->Clear();
  c_uncut_cons_pz->Draw(); c0->SaveAs("histos/c_uncut_cons_pz.pdf"); c0->Clear();
  c_uncut_cons_E->Draw(); c0->SaveAs("histos/c_uncut_cons_E.pdf"); c0->Clear();
  c_uncut_cons_phi->Draw(); c0->SaveAs("histos/c_uncut_cons_phi.pdf"); c0->Clear();
  c_uncut_cons_eta->Draw(); c0->SaveAs("histos/c_uncut_cons_eta.pdf"); c0->Clear();
  c_uncut_cons_Pt->Draw(); c0->SaveAs("histos/c_uncut_cons_Pt.pdf"); c0->Clear();

  cut_px->Draw(); c0->SaveAs("histos/cut_px.pdf"); c0->Clear();
  cut_py->Draw(); c0->SaveAs("histos/cut_py.pdf"); c0->Clear();
  cut_pz->Draw(); c0->SaveAs("histos/cut_pz.pdf"); c0->Clear();
  cut_E->Draw(); c0->SaveAs("histos/cut_E.pdf"); c0->Clear();
  cut_phi->Draw(); c0->SaveAs("histos/cut_phi.pdf"); c0->Clear();
  cut_eta->Draw(); c0->SaveAs("histos/cut_eta.pdf"); c0->Clear();
  cut_leadPt->Draw(); c0->SaveAs("histos/cut_leadPt.pdf"); c0->Clear();
  cut_Pt->Draw(); c0->SaveAs("histos/cut_Pt.pdf"); c0->Clear();

  c_cut_px->Draw(); c0->SaveAs("histos/c_cut_px.pdf"); c0->Clear();
  c_cut_py->Draw(); c0->SaveAs("histos/c_cut_py.pdf"); c0->Clear();
  c_cut_pz->Draw(); c0->SaveAs("histos/c_cut_pz.pdf"); c0->Clear();
  c_cut_E->Draw(); c0->SaveAs("histos/c_cut_E.pdf"); c0->Clear();
  c_cut_phi->Draw(); c0->SaveAs("histos/c_cut_phi.pdf"); c0->Clear();
  c_cut_eta->Draw(); c0->SaveAs("histos/c_cut_eta.pdf"); c0->Clear();
  c_cut_leadPt->Draw(); c0->SaveAs("histos/c_cut_leadPt.pdf"); c0->Clear();
  c_cut_Pt->Draw(); c0->SaveAs("histos/c_cut_Pt.pdf"); c0->Clear();

  cut2_px->Draw(); c0->SaveAs("histos/cut2_px.pdf"); c0->Clear();
  cut2_py->Draw(); c0->SaveAs("histos/cut2_py.pdf"); c0->Clear();
  cut2_pz->Draw(); c0->SaveAs("histos/cut2_pz.pdf"); c0->Clear();
  cut2_E->Draw(); c0->SaveAs("histos/cut2_E.pdf"); c0->Clear();
  cut2_phi->Draw(); c0->SaveAs("histos/cut2_phi.pdf"); c0->Clear();
  cut2_eta->Draw(); c0->SaveAs("histos/cut2_eta.pdf"); c0->Clear();
  cut2_leadPt->Draw(); c0->SaveAs("histos/cut2_leadPt.pdf"); c0->Clear();
  cut2_Pt->Draw(); c0->SaveAs("histos/cut2_Pt.pdf"); c0->Clear();

  c_cut2_px->Draw(); c0->SaveAs("histos/c_cut2_px.pdf"); c0->Clear();
  c_cut2_py->Draw(); c0->SaveAs("histos/c_cut2_py.pdf"); c0->Clear();
  c_cut2_pz->Draw(); c0->SaveAs("histos/c_cut2_pz.pdf"); c0->Clear();
  c_cut2_E->Draw(); c0->SaveAs("histos/c_cut2_E.pdf"); c0->Clear();
  c_cut2_phi->Draw(); c0->SaveAs("histos/c_cut2_phi.pdf"); c0->Clear();
  c_cut2_eta->Draw(); c0->SaveAs("histos/c_cut2_eta.pdf"); c0->Clear();
  c_cut2_leadPt->Draw(); c0->SaveAs("histos/c_cut2_leadPt.pdf"); c0->Clear();
  c_cut2_Pt->Draw(); c0->SaveAs("histos/c_cut2_Pt.pdf"); c0->Clear();

  effic_px->Draw(); c0->SaveAs("histos/effic_px.pdf"); c0->Clear();
  effic_py->Draw(); c0->SaveAs("histos/effic_py.pdf"); c0->Clear();
  effic_pz->Draw(); c0->SaveAs("histos/effic_pz.pdf"); c0->Clear();
  effic_E->Draw(); c0->SaveAs("histos/effic_E.pdf"); c0->Clear();
  effic_phi->Draw(); c0->SaveAs("histos/effcc_phi.pdf"); c0->Clear();
  effic_eta->Draw(); c0->SaveAs("histos/effic_eta.pdf"); c0->Clear();
  effic_leadPt->Draw(); c0->SaveAs("histos/effic_leadPt.pdf"); c0->Clear();
  effic_Pt->Draw(); c0->SaveAs("histos/effic_Pt.pdf"); c0->Clear();
  test_ptdiff->Draw(); c0->SaveAs("histos/test_ptdiff.pdf"); c0->Clear();
  test_without_leadsublead->Draw(); c0->SaveAs("histos/test_without_leadsublead.pdf"); c0->Clear();

  c_effic_px->Draw(); c0->SaveAs("histos/c_effic_px.pdf"); c0->Clear();
  c_effic_py->Draw(); c0->SaveAs("histos/c_effic_py.pdf"); c0->Clear();
  c_effic_pz->Draw(); c0->SaveAs("histos/c_effic_pz.pdf"); c0->Clear();
  c_effic_E->Draw(); c0->SaveAs("histos/c_effic_E.pdf"); c0->Clear();
  c_effic_phi->Draw(); c0->SaveAs("histos/c_effic_phi.pdf"); c0->Clear();
  c_effic_eta->Draw(); c0->SaveAs("histos/c_effic_eta.pdf"); c0->Clear();
  c_effic_leadPt->Draw(); c0->SaveAs("histos/c_effic_leadPt.pdf"); c0->Clear();
  c_effic_Pt->Draw(); c0->SaveAs("histos/c_effic_Pt.pdf"); c0->Clear();
*/
  num_dif_raw->Draw(); c0->SaveAs("histos/num_diff_raw.pdf"); c0->Clear();
  c_num_dif_raw->Draw(); c0->SaveAs("histos/c_num_diff_raw.pdf"); c0->Clear();
  num_dif->Draw(); c0->SaveAs("histos/num_diff.pdf"); c0->Clear();
  c_num_dif->Draw(); c0->SaveAs("histos/c_num_diff.pdf"); c0->Clear();
  

  return;
}
