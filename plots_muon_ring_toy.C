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

Int_t plots_muon_ring_toy(){
  //
  TRandom3 *rnd = new TRandom3(1231231);
  //
  Int_t nn = 1000000;

  Double_t r_sim = 0.7;
  Double_t r_sim_sigma = 0.1;
  Double_t x;
  Double_t y;

  Double_t x_min = -2.5;
  Double_t x_max =  2.5;
  Double_t y_min = -2.5;
  Double_t y_max =  2.5;

  Double_t xc =  1.0;
  Double_t yc =  1.0;

  Double_t rpdf;

  TH2D *h2 = new TH2D("h2","h2",100,x_min-0.1,x_max+0.1,100,y_min-0.1,y_max+0.1);
  //h2->SetMaximum(1.0);
  
  for(Int_t i = 0;i<nn;i++){
    x = rnd->Uniform(x_min,x_max);
    y = rnd->Uniform(y_min,y_max);
    rpdf = TMath::Sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
    h2->Fill(x,y,TMath::Gaus(rpdf,r_sim,r_sim_sigma)*(y*y + 0.1));
  }
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,900,900);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  h2->Draw("ZCOLOR");
  //
  return 0;
}
