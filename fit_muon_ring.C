//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

#include <time.h>

#include "TMinuit.h"

using namespace std;

void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R);
void get_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Double_t &x0, Double_t &y0, Double_t &R);
void get_Rvs_theta_and_theta_dist(TGraph *gr, Double_t x0, Double_t y0, Double_t R, TGraph *gr_R, TH1D *h1_theta_deg);
void gen_Graph3D(TGraph *gr,TGraph2DErrors *gr2D);
Double_t plane_3D(Double_t *x, Double_t *par);
TGraph *_gr_to_fit = new TGraph();
double equation_of_circle(double x, double y, double *par);
void fcn(int &npar, double *gin, double &f, double *par, int iflag);

void fit_ring_with_Minuit(Double_t x0in, Double_t y0in, Double_t Rin,
			  Double_t &x0out, Double_t &y0out, Double_t &Rout,
			  Double_t &x0outerr, Double_t &y0outerr, Double_t &Routerr);

Int_t fit_muon_ring(){
  //
  TRandom3 *rnd = new TRandom3(12312312);
  //
  TString fileN01;
  fileN01 = "../scratch/simtel_data/muon/hist/hist_run1_muon.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  //TGraph *gr = (TGraph*)f01->Get("trueRing/gr_58");
  TGraph *gr = (TGraph*)f01->Get("trueRing/gr_4");
  TGraph *gr_frame = (TGraph*)f01->Get("gr_frame");
  //
  //
  Double_t xx, yy;
  for(Int_t i=0;i<gr->GetN(); i++){
    gr->GetPoint( i, xx, yy);
    _gr_to_fit->SetPoint( i, xx, yy);
  }
  //
  Double_t x0app, y0app, Rapp;
  Double_t x0app_average = 0.0;
  Double_t y0app_average = 0.0;
  Double_t Rapp_average = 0.0;
  vector<Double_t> x0app_v;
  vector<Double_t> y0app_v;
  vector<Double_t> Rapp_v;
  for(Int_t i = 0;i<50;i++){
    get_approximate_ring_parameters( gr, rnd, x0app, y0app, Rapp);
    x0app_v.push_back(x0app);
    y0app_v.push_back(y0app);
    Rapp_v.push_back(Rapp);
    x0app_average += x0app;
    y0app_average += y0app;
    Rapp_average += Rapp;
    //cout<<"x0app "<<x0app<<endl
    //    <<"y0app "<<y0app<<endl
    //    <<"Rapp  "<<Rapp<<endl;
  }
  //
  x0app_average /= x0app_v.size();
  y0app_average /= x0app_v.size();
  Rapp_average /= x0app_v.size();
  //
  Double_t x0out, y0out, Rout;
  Double_t x0outerr, y0outerr, Routerr;
  fit_ring_with_Minuit( x0app_average, y0app_average, Rapp_average,
			x0out, y0out, Rout,
			x0outerr, y0outerr, Routerr);
  //
  TGraph *gr_app_average_r0 = new TGraph();
  gr_app_average_r0->SetNameTitle("gr_app_average_r0","gr_app_average_r0");  
  gr_app_average_r0->SetPoint( 0, x0app_average, y0app_average);
  gr_app_average_r0->SetMarkerStyle(43);
  gr_app_average_r0->SetMarkerColor(kMagenta+3);
  gr_app_average_r0->SetMarkerSize(3.0);
  //
  TGraph *gr_app_average_ring = new TGraph();
  gr_app_average_ring->SetNameTitle("gr_app_average_ring","gr_app_average_ring");  
  gr_app_average_ring->SetMarkerStyle(7);
  gr_app_average_ring->SetMarkerColor(kMagenta+3);
  gr_app_average_ring->SetMarkerSize(3.0);
  gen_ring(gr_app_average_ring, 360, x0app_average, y0app_average, Rapp_average);
  //
  TGraph *gr_fit_ring = new TGraph();
  gr_fit_ring->SetNameTitle("gr_fit_ring","gr_fit_ring");
  gr_fit_ring->SetMarkerStyle(7);
  gr_fit_ring->SetMarkerColor(kRed);
  gr_fit_ring->SetMarkerSize(2.0);
  gen_ring(gr_fit_ring, 360, x0out, y0out, Rout);
  //
  TGraph *gr_app_r0 = new TGraph();
  gr_app_r0->SetNameTitle("gr_app_r0","gr_app_r0");
  for(unsigned int ii = 0;ii<x0app_v.size();ii++)
    gr_app_r0->SetPoint( ii, x0app_v.at(ii), y0app_v.at(ii));
  //
  gr->SetMarkerStyle(20);
  gr_frame->SetMarkerStyle(25);
  gr_app_r0->SetMarkerStyle(43);
  gr_app_r0->SetMarkerColor(kBlue);
  gr_app_r0->SetMarkerSize(2.0);
  //  
  TGraph *gr_R_app_average = new TGraph();
  gr_R_app_average->SetNameTitle("gr_R_app_average","gr_R_app_average");
  TH1D *h1_theta_app_average_deg = new TH1D("h1_theta_app_average_deg","h1_theta_app_average_deg", 36, 0.0, 360.0); 
  get_Rvs_theta_and_theta_dist( gr, x0app_average, y0app_average, Rapp_average, gr_R_app_average, h1_theta_app_average_deg);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1800,600);
  c1->Divide(3,1);
  c1->cd(1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_frame);
  mg->Add(gr_app_r0);  
  mg->Add(gr_app_average_r0);
  mg->Add(gr_app_average_ring);
  mg->Add(gr_fit_ring);
  mg->Draw("AP");
  //
  c1->cd(2);
  gr_R_app_average->SetMarkerStyle(20);
  gr_R_app_average->Draw("AP");  
  //
  c1->cd(3);
  h1_theta_app_average_deg->SetLineColor(kBlack);
  h1_theta_app_average_deg->SetLineWidth(2.0);
  h1_theta_app_average_deg->Draw();
  //
  //
  //
  //
  //
  TGraph2DErrors *gr2D = new TGraph2DErrors();
  TGraph2DErrors *gr2D_average = new TGraph2DErrors();
  gen_Graph3D( gr,gr2D);
  gen_Graph3D( gr_app_average_ring,gr2D_average);
  //
  //
  Double_t boxL = 5;
  TGraph2D *gr2D_box = new TGraph2D();
  gr2D_box->SetTitle("gr2D_box");
  gr2D_box->SetName("gr2D_box");
  //
  gr2D_box->SetPoint(gr2D_box->GetN(),-boxL/2.0,-boxL/2.0,-boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),-boxL/2.0, boxL/2.0,-boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),-boxL/2.0,-boxL/2.0, boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),-boxL/2.0, boxL/2.0, boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),boxL/2.0,-boxL/2.0,-boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),boxL/2.0, boxL/2.0,-boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),boxL/2.0,-boxL/2.0, boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),boxL/2.0, boxL/2.0, boxL/2.0);
  //
  //
  /*
  TCanvas *c2 = new TCanvas("c2","c2",10,10,1800,600);
  c2->Divide(3,1);
  c2->cd(1);
  gr2D->SetTitle("");
  gr2D_box->SetTitle("");
  gr2D_box->Draw("P");
  gr2D->Draw("sameP");
  //
  c2->cd(2);
  gr2D_average->SetTitle("");
  gr2D_box->Draw("P");
  gr2D_average->Draw("sameP");
  //
  c2->cd(3);
  gr2D_average->SetTitle("");
  gr2D_box->Draw("P");
  gr2D_average->Draw("sameP");
  gr2D->Draw("sameP");
  */
  //
  //
  //Double_t x_fit_min = -1.2;
  //Double_t x_fit_max = -1.2;
  //Double_t y_fit_min = -1.2;
  //Double_t y_fit_max = -1.2;
  //Int_t ndim = 2;
  //
  //TF2 TF2(const char* name, Double_t(*)(Double_t*,Double_t*) fcn, Double_t xmin = 0, Double_t xmax = 1, Double_t ymin = 0, Double_t ymax = 1, Int_t npar = 0, Int_t ndim = 2)
  //TF2 *plane_fit_f = new TF2("plane_fit_f", plane_3D, x_fit_min, x_fit_max, y_fit_min, y_fit_max, npar, ndim);
  //plane_fit_f->SetParameter(0,x0app_average);
  //plane_fit_f->SetParameter(1,y0app_average);
  //plane_fit_f->SetParameter(2,(x0app_average*x0app_average + y0app_average*y0app_average - Rapp_average*Rapp_average));
  //plane_fit_f->FixParameter(0,x0app_average);
  //plane_fit_f->FixParameter(1,y0app_average);
  //
  //gr2D_average->Fit(plane_fit_f,"");
  //gr2D->Fit(plane_fit_f);
  //
  //
  return 0;
}

void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R){
  Double_t phi = 0.0;
  TVector2 rc(x0,y0);
  for(Int_t i = 0;i<np;i++){
    TVector2 p;
    p.SetMagPhi(R,2*TMath::Pi()/(np-1)*i);
    TVector2 pt = rc + p;
    gr->SetPoint( i, pt.X(), pt.Y());
  }
}

void get_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Double_t &x0, Double_t &y0, Double_t &R){
  Double_t p1x0, p1y0;
  Double_t p2x0, p2y0; 
  Double_t maxDist = 0.0;
  Int_t p1_id = (Int_t)rnd->Uniform(0,gr->GetN());
  Int_t p2_id = (Int_t)rnd->Uniform(0,gr->GetN());
  gr->GetPoint(p1_id, p1x0, p1y0);
  TVector2 p1(p1x0, p1y0);
  gr->GetPoint(p2_id, p2x0, p2y0);
  TVector2 p2(p2x0, p2y0);
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i, p2x0, p2y0);
    p2.Set(p2x0, p2y0);
    TVector2 dp = p1 - p2;
    if(maxDist<dp.Mod()){
      maxDist = dp.Mod();
      p2_id = i;
    }
  }
  //
  gr->GetPoint(p2_id, p2x0, p2y0);
  p2.Set(p2x0, p2y0);
  //  
  TVector2 pm = (p1 + p2)/2.0;
  x0 = pm.X();
  y0 = pm.Y();
  R = maxDist/2.0;
}

void get_Rvs_theta_and_theta_dist(TGraph *gr, Double_t x0, Double_t y0, Double_t R, TGraph *gr_R, TH1D *h1_theta_deg){
  Double_t x, y;
  TVector2 pc(x0, y0);
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i, x, y);
    TVector2 p(x,y);
    TVector2 dp = p - pc;
    gr_R->SetPoint(i,dp.Phi()*180.0/TMath::Pi(),(dp.Mod() - R)/R);
    h1_theta_deg->Fill(dp.Phi()*180.0/TMath::Pi());
  }
}

void gen_Graph3D(TGraph *gr,TGraph2DErrors *gr2D){
  Double_t x, y;
  Double_t xt, yt,zt;
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i,x,y);
    xt = -2.0*x;
    yt = -2.0*y;
    zt = x*x + y*y;
    gr2D->SetPoint(i,xt,yt,zt);
    gr2D->SetPointError(i,100,100,100);
  }
}

Double_t plane_3D(Double_t *x, Double_t *par){
  //
  //Double_t X0 = par[0];
  //Double_t Y0 = par[1];
  //Double_t Z0 = 1.0;
  //Double_t D  = par[2];
  //
  //Double_t x = x[0];
  //Double_t y = x[1];
  //
  // X0x + Y0y + 1.0*z + D = 0
  return -par[0]*x[0] - par[1]*x[2] - par[2];
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag){
   double chisq = 0;
   double x, y;
   double delta;
   for (int i = 0; i<_gr_to_fit->GetN(); i++){
     _gr_to_fit->GetPoint( i, x, y);
     delta = equation_of_circle(x, y, par);
     if(delta>0.0)
       delta = delta*1.0;
     chisq += TMath::Abs(delta);
   }
   f = chisq;
}

double equation_of_circle(double x, double y, double *par){
  return par[2]*par[2] - (x-par[0])*(x-par[0]) - (y-par[1])*(y-par[1]);
}

void fit_ring_with_Minuit(Double_t x0in, Double_t y0in, Double_t Rin,
			  Double_t &x0out, Double_t &y0out, Double_t &Rout,
			  Double_t &x0outerr, Double_t &y0outerr, Double_t &Routerr){
  //
  Int_t npar = 3;
  TMinuit *gMinuit = new TMinuit(npar);
  gMinuit->SetPrintLevel(-1.0);
  gMinuit->SetFCN(fcn); 
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // 
  // Set starting values and step sizes for parameters
  gMinuit->mnparm(0, "x0", x0in, 0.001, 0,0,ierflg);
  gMinuit->mnparm(1, "y0", y0in, 0.001, 0,0,ierflg);
  gMinuit->mnparm(2, "R", Rin, 0.001, 0,0,ierflg);
  //

  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  // Print results
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);
  //
  gMinuit->GetParameter(0, x0out, x0outerr);
  gMinuit->GetParameter(1, y0out,  y0outerr);
  gMinuit->GetParameter(2, Rout, Routerr);
  //
  cout<<x0out<<endl
      <<y0out<<endl
      <<Rout<<endl;
  //
}
