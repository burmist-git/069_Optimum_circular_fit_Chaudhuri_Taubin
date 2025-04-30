//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

using namespace std;

TString getObjectName(TString namein, Int_t i);

Int_t plots_diff(){

  TString fileN01;
  TString fileN02;
  TString fileN03;
  //
  fileN01 = "./histOut_usual.root";
  fileN02 = "./histOut_taubin_circle_fit.root";
  //fileN03 = "./histOut_taubin_circle_fit_mod.root";

  TFile *f1 = new TFile(fileN01.Data());
  TFile *f2 = new TFile(fileN02.Data());
  //TFile *f3 = new TFile(fileN03.Data());

  TH1D *h1_01 = (TH1D*)f1->Get("h1_r_pull");
  TH1D *h1_02 = (TH1D*)f2->Get("h1_r_pull");
  //TH1D *h1_03 = (TH1D*)f3->Get("h1_r_pull");
  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  // 
  h1_01->SetLineColor(kBlack);
  h1_01->SetLineWidth(3.0);
  h1_01->SetMarkerColor(kBlack);
  //
  h1_02->SetLineColor(kRed+2);
  h1_02->SetLineWidth(3.0);
  h1_02->SetMarkerColor(kRed+2);
  //
  //h1_03->SetLineColor(kBlue+2);
  //h1_03->SetLineWidth(3.0);
  //h1_03->SetMarkerColor(kBlue+2);
  //
  h1_01->SetMaximum(500);
  h1_01->Draw();
  h1_02->Draw("sames");
  //h1_03->Draw("sames");
  h1_01->GetXaxis()->SetTitle("r_true - r_reco");
  //h1_01->GetYaxis()->SetTitle("FADC counts");
  //  
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "Chaudhuri", "apl");
  leg->AddEntry(h1_02, "Taubin", "apl");
  leg->Draw();  

  return 0;
}
