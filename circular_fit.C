TRandom3 *rnd = new TRandom3(12312);
TGraph *_gr_data = new TGraph();
TGraph *_gr_weight = new TGraph();
Int_t _n_points = 0;

void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R);
void get_optimum_circular_fit_Chaudhuri( TGraph *gr_data, TGraph *gr_weight,  Double_t &cx_reco,  Double_t &cy_reco,  Double_t &r_reco);
void get_weight_data_from_TH2D( TH2D *h2, TGraph *gr_data, TGraph *gr_weight);
void generate_uniformly_filled_ring(TRandom3* rnd, Double_t cx, Double_t cy, Double_t r, Double_t &x, Double_t &y);

void fcn(int &npar, double *gin, double &f, double *par, int iflag);
double equation_of_circle(double x, double y, double *par);
double r2_of_circle(double x, double y, double *par);
void fit_ring_with_Minuit(Double_t x0in, Double_t y0in, Double_t Rin,
			  Double_t &x0out, Double_t &y0out, Double_t &Rout,
			  Double_t &x0outerr, Double_t &y0outerr, Double_t &Routerr);
void fill_gr_data_gr_weight(TGraph *gr_data, TGraph *gr_weight);

void run_circular_fit(TCanvas *c1,
		      Double_t &cx_reco,
		      Double_t &cy_reco,
		      Double_t &r_reco,
		      TString out_name = " ",
		      Int_t nn_pe = 400,
		      Int_t nn_pe_nsb = 0,
		      Int_t nn_pe_nsb_inside = 0,
		      Double_t cx_true = -2.0,
		      Double_t cy_true = 1.0,
		      Double_t r_true  = 1.2,
		      Double_t sigma   = 0.1,
		      Double_t err_lin = 0.1,
		      Double_t phi0    = TMath::Pi()-20.0/180.0*TMath::Pi(),
		      Double_t phi_sigma = 0.15*TMath::Pi());

Int_t circular_fit(){

  Int_t nn_pe = 400;
  Int_t nn_pe_nsb = 0;
  Int_t nn_pe_nsb_inside = 0;
  Double_t cx_true = -2.0;
  Double_t cy_true = 1.0;
  Double_t r_true  = 1.2;
  Double_t sigma   = 0.00;
  Double_t err_lin = 0.05;
  Double_t phi0    = TMath::Pi()-20.0/180.0*TMath::Pi();
  Double_t phi_sigma = 0.3*TMath::Pi();
  //
  Double_t cx_reco;
  Double_t cy_reco;
  Double_t r_reco;

  TH1D *h1_cx_pull = new TH1D("h1_cx_pull","h1_cx_pull",2000,-1.0,1.0);
  TH1D *h1_cy_pull = new TH1D("h1_cy_pull","h1_cy_pull",2000,-1.0,1.0);
  TH1D *h1_r_pull = new TH1D("h1_r_pull","h1_r_pull",2000,-1.0,1.0);

  Int_t nn_tests = 100;

  TCanvas *c1;
  
  for(Int_t i = 0;i<nn_tests;i++){
    TString c1_name = "c1_";
    c1_name += i;
    c1_name += "ev";
    c1 = new TCanvas(c1_name.Data(),c1_name.Data(), 10, 10, 1000, 1000);
    
    cx_true = rnd->Uniform(-3.0,3.0);
    cy_true = rnd->Uniform(-3.0,3.0);
    r_true  = rnd->Uniform(0.9,1.2);
    phi0    = TMath::Pi() + rnd->Uniform(-180.0,180.0)/180.0*TMath::Pi();
    //
    //
    TString out_name = "./out_pdf/";
    out_name += "c1_";
    out_name += i;
    out_name += "ev.pdf";
    //out_name = " ";
    run_circular_fit(c1,
		     cx_reco, cy_reco, r_reco,
		     out_name,
		     nn_pe,
		     nn_pe_nsb,
		     nn_pe_nsb_inside,
		     cx_true,
		     cy_true,
		     r_true,
		     sigma,
		     err_lin,
		     phi0,
		     phi_sigma);
    //
    h1_cx_pull->Fill(cx_true - cx_reco);
    h1_cy_pull->Fill(cy_true - cy_reco);
    h1_r_pull->Fill((r_true - r_reco));
    //
    if(i<10){
      cout<<setw(20)<<"cx_reco"<<setw(20)<<"cy_reco"<<setw(20)<<"r_reco"<<endl;
      cout<<setw(20)<<cx_reco<<setw(20)<<cy_reco<<setw(20)<<r_reco<<endl
	  <<setw(20)<<cx_true<<setw(20)<<cy_true<<setw(20)<<r_true<<endl;
    }
    delete c1;
    c1 = NULL;
  }
  
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);

  if(c1 != NULL)
    c1->Draw();

  //
  //
  
  TCanvas *c2 = new TCanvas("c2","c2", 10, 10, 1800, 600);
  c2->Divide(3,1);
  c2->cd(1);
  h1_cx_pull->Draw();
  c2->cd(2);
  h1_cy_pull->Draw();
  c2->cd(3);
  h1_r_pull->Draw();

  TString histOut = "histOut.root";
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;

  //
  //
  c2->Write();
  h1_cx_pull->Write();
  h1_cy_pull->Write();
  h1_r_pull->Write();

  rootFile->Close();

  return 0;
}

void run_circular_fit(TCanvas *c1,
		      Double_t &cx_reco,
		      Double_t &cy_reco,
		      Double_t &r_reco,
		      TString out_name,
		      Int_t nn_pe,
		      Int_t nn_pe_nsb,
		      Int_t nn_pe_nsb_inside,
		      Double_t cx_true,
		      Double_t cy_true,
		      Double_t r_true,
		      Double_t sigma,
		      Double_t err_lin,
		      Double_t phi0,
		      Double_t phi_sigma){
  
  //Double_t frame_xmin = cx_true - 1.5*r_true - 3.0*sigma;
  //Double_t frame_xmax = cx_true + 1.5*r_true + 3.0*sigma;
  //Double_t frame_ymin = cy_true - 1.5*r_true - 3.0*sigma;
  //Double_t frame_ymax = cy_true + 1.5*r_true + 3.0*sigma;

  Double_t frame_xmin = -5.0;
  Double_t frame_xmax =  5.0;
  Double_t frame_ymin = -5.0;
  Double_t frame_ymax =  5.0;  

  Double_t cx_reco_err, cy_reco_err, r_reco_err;
  Double_t cx_seed, cy_seed, r_seed;
  
  Double_t x, y;
  Double_t phi, rho;

  TGraph *gr_frame = new TGraph();
  gr_frame->SetNameTitle("gr_frame","gr_frame");
  gr_frame->SetPoint(0,frame_xmin,  frame_ymin);
  gr_frame->SetPoint(1,frame_xmin,  frame_ymax);
  gr_frame->SetPoint(2,frame_xmax,  frame_ymax);
  gr_frame->SetPoint(3,frame_xmax,  frame_ymin);
  gr_frame->SetPoint(4,frame_xmin,  frame_ymin);  
  TGraph *gr_data = new TGraph();
  gr_data->SetNameTitle("gr_data","gr_data");
  TGraph *gr_weight = new TGraph();
  gr_weight->SetNameTitle("gr_weight","gr_weight");
  //
  TGraph *gr_data_from_h2 = new TGraph();
  gr_data_from_h2->SetNameTitle("gr_data_from_h2","gr_data_from_h2");
  TGraph *gr_weight_from_h2 = new TGraph();
  gr_weight_from_h2->SetNameTitle("gr_weight_from_h2","gr_weight_from_h2");
  //
  TGraph *gr_reco = new TGraph();
  gr_reco->SetNameTitle("gr_reco","gr_reco");
  TGraph *gr_reco_c = new TGraph();
  gr_reco_c->SetNameTitle("gr_reco_c","gr_reco_c");
  TGraph *gr_true = new TGraph();
  gr_true->SetNameTitle("gr_true","gr_true");
  TGraph *gr_true_c = new TGraph();
  gr_true_c->SetNameTitle("gr_true_c","gr_true_c");

  gr_frame->SetMarkerStyle(1);
  gr_data->SetMarkerStyle(7);
  gr_reco->SetMarkerStyle(7);
  gr_reco->SetMarkerColor(kRed);
  gr_reco_c->SetMarkerStyle(65);
  gr_reco_c->SetMarkerSize(1);
  gr_reco_c->SetMarkerColor(kRed);
  gr_true->SetMarkerStyle(7);
  gr_true->SetMarkerColor(kBlue+2);
  gr_true_c->SetMarkerStyle(65);
  gr_true_c->SetMarkerSize(1);
  gr_true_c->SetMarkerColor(kBlue+2);
  
  TH2D *h2_data = new TH2D("h2_data","h2_data",70,frame_xmin,frame_xmax,70,frame_ymin,frame_ymax);
  TH1D *h1_rho = new TH1D("h1_rho","h1_rho",200,0.0,2.0*frame_xmax);
  TH1D *h1_phi = new TH1D("h1_phi","h1_phi",200,-0.05*2.0*TMath::Pi(),1.05*2.0*TMath::Pi());
  
  //
  //
  TVector2 v_tmp;
  for(Int_t i = 0; i<nn_pe; i++){
    //phi = rnd->Uniform(0,2.0*TMath::Pi());
    phi = rnd->Gaus(TMath::Pi()+phi0,phi_sigma);
    rho = rnd->Gaus(r_true,sigma);
    v_tmp.SetMagPhi(rho,phi);
    x = v_tmp.X() + cx_true + rnd->Uniform(-err_lin,err_lin);
    y = v_tmp.Y() + cy_true + rnd->Uniform(-err_lin,err_lin);
    //rho = rnd->Uniform(r_true-sigma,r_true+sigma);
    h1_rho->Fill(rho);
    h1_phi->Fill(phi);
    gr_data->SetPoint(gr_data->GetN(),x,y);
    gr_weight->SetPoint(gr_weight->GetN(),1.0,1.0);
    h2_data->Fill(x,y);
  }

  for(Int_t i = 0; i<nn_pe_nsb; i++){
    x = rnd->Uniform(frame_xmin,frame_xmax);
    y = rnd->Uniform(frame_ymin,frame_ymax);
    gr_data->SetPoint(gr_data->GetN(),x,y);
    gr_weight->SetPoint(gr_weight->GetN(),1.0,1.0);
    h2_data->Fill(x,y);
  }

  for(Int_t i = 0; i<nn_pe_nsb_inside; i++){
    generate_uniformly_filled_ring(rnd, cx_true,cy_true,r_true, x, y);
    gr_data->SetPoint(gr_data->GetN(),x,y);
    gr_weight->SetPoint(gr_weight->GetN(),1.0,1.0);
    h2_data->Fill(x,y);
  }
  
  get_weight_data_from_TH2D( h2_data, gr_data_from_h2, gr_weight_from_h2);
  //
  //fill_gr_data_gr_weight(gr_data, gr_weight);
  fill_gr_data_gr_weight(gr_data_from_h2, gr_weight_from_h2);
  
  //get_optimum_circular_fir_Chaudhuri(gr_data, gr_weight, cx_reco, cy_reco, r_reco);
  //get_optimum_circular_fit_Chaudhuri(gr_data_from_h2, gr_weight_from_h2, cx_reco, cy_reco, r_reco);

  cx_seed = 0.0;
  cy_seed = 0.0;
  r_seed = 1.0;
  fit_ring_with_Minuit(cx_seed, cy_seed, r_seed,
		       cx_reco, cy_reco, r_reco,
		       cx_reco_err, cy_reco_err, r_reco_err);

  
  gen_ring(gr_reco, 300, cx_reco, cy_reco, r_reco);

  gen_ring(gr_true, 300, cx_true, cy_true, r_true);

  gr_true_c->SetPoint(0,cx_true, cy_true);
  gr_reco_c->SetPoint(0,cx_reco, cy_reco);
  

  if(c1!=NULL){
    c1->Divide(2,2);
    gStyle->SetPalette(1);
    gStyle->SetFrameBorderMode(0);
    gROOT->ForceStyle();
    gStyle->SetStatColor(kWhite);
    gStyle->SetOptStat(kFALSE); 
    
    c1->cd(1);
    
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr_frame);
    mg->Add(gr_data);
    mg->Add(gr_reco);
    mg->Add(gr_reco_c);
    mg->Add(gr_true);
    mg->Add(gr_true_c);
    mg->Draw("ap");
    
    c1->cd(2);
    h2_data->Draw("ZCOLOR");
    gr_reco->Draw("same");
    gr_reco_c->Draw("same");
    gr_true->Draw("same");
    gr_true_c->Draw("same");
    
    c1->cd(3);
    h1_rho->Draw();
    
    c1->cd(4);
    h1_phi->Draw();

    if(out_name != " ")
      c1->SaveAs(out_name.Data());
  }

  delete h2_data;
  delete h1_rho;
  delete h1_phi;
  delete gr_data_from_h2;
  delete gr_weight_from_h2;

}

void get_optimum_circular_fit_Chaudhuri(TGraph *gr_data, TGraph *gr_weight, Double_t &cx_reco,  Double_t &cy_reco,  Double_t &r_reco){
  Double_t x, y;
  Double_t M_tot = 0.0;
  Double_t M_w_tot = 0.0;
  Double_t x_mean = 0.0;
  Double_t y_mean = 0.0;
  Double_t w = 0.0;
  Double_t A = 0.0;
  Double_t A_ = 0.0;
  Double_t B = 0.0;
  Double_t B_ = 0.0;
  Double_t C = 0.0;
  Double_t C_ = 0.0;
  for(Int_t i = 0; i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    gr_weight->GetPoint(i,w,w);
    x_mean += x*w;
    y_mean += y*w;
    M_w_tot += w;
    M_tot++;
  }
  x_mean /= M_w_tot;
  y_mean /= M_w_tot; 
  //
  for(Int_t i = 0; i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    gr_weight->GetPoint(i,w,w);
    A  += w*(x - x_mean)*x;
    A_ += w*(y - y_mean)*x;
    B  += w*(x - x_mean)*y;
    B_ += w*(y - y_mean)*y;
    C  += w*(x - x_mean)*(x*x + y*y);
    C_ += w*(y - y_mean)*(x*x + y*y);
  }
  C  /= 2;
  C_ /= 2;
  //
  cx_reco = (B_*C - B*C_) / (A*B_ - A_*B);
  cy_reco = (A_*C - A*C_) / (A_*B - A*B_);
  //
  r_reco = 0.0;
  for(Int_t i = 0; i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    gr_weight->GetPoint(i,w,w);
    r_reco += w*((x - cx_reco)*(x - cx_reco) + (y - cy_reco)*(y - cy_reco));
  }  
  r_reco /= M_w_tot;
  r_reco = TMath::Sqrt(r_reco);
  //cx_reco = x_mean;
  //cy_reco = y_mean;
  //gr_data->GetPoint(0,x,y);
  //r_reco = TMath::Sqrt((x - cx_reco)*(x - cx_reco) + (y - cy_reco)*(y - cy_reco));    
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

void generate_uniformly_filled_ring(TRandom3* rnd, Double_t cx, Double_t cy, Double_t r, Double_t &x, Double_t &y){
  bool go = false;
  Double_t cx_min = cx - r;
  Double_t cx_max = cx + r;
  Double_t cy_min = cy - r;
  Double_t cy_max = cy + r;
  while( go == false ){
    x = rnd->Uniform(cx_min, cx_max);
    y = rnd->Uniform(cy_min, cy_max);
    if(r*r>=((x - cx)*(x - cx) + (y - cy)*(y - cy)))
      return;
  }
}

void get_weight_data_from_TH2D( TH2D *h2, TGraph *gr_data, TGraph *gr_weight){
  Double_t x;
  Double_t y;
  Double_t w;
  for(Int_t i = 1; i <= h2->GetNbinsX(); i++){
    for(Int_t j = 1; j<= h2->GetNbinsY(); j++){
      x = h2->GetXaxis()->GetBinCenter(i);
      y = h2->GetYaxis()->GetBinCenter(j);
      w = h2->GetBinContent(i,j);
      if(w>0){
	gr_data->SetPoint(gr_data->GetN(),x,y);
	gr_weight->SetPoint(gr_weight->GetN(),w,w);
      }
      //cout<<"x "<<x<<endl
      //  <<"y "<<y<<endl
      //  <<"w "<<w<<endl;
    }
  }
  //cout<<"h2->GetNbinsY() "<<h2->GetNbinsY()<<endl;
}

/*
void fcn(int &npar, double *gin, double &f, double *par, int iflag){
   double chisq = 0;
   double x, y;
   double delta;
   for (int i = 0; i<_n_points; i++){
     _gr_data->GetPoint( i, x, y);
     delta = equation_of_circle(x, y, par);
     if(delta>0.0)
       delta = delta*0.8;
     chisq += TMath::Abs(delta);
   }
   f = chisq;
}
*/

void fcn(int &npar, double *gin, double &f, double *par, int iflag){
   double chisq = 0;
   double x, y;
   double delta;
   double delta_sum = 0;
   double r2_sum = 0;
   for (int i = 0; i<_n_points; i++){
     _gr_data->GetPoint( i, x, y);
     delta = equation_of_circle(x, y, par);
     delta *= delta;
     delta_sum += delta;
     r2_sum += r2_of_circle(x, y, par);
   }
   chisq = delta_sum/4.0/(r2_sum/_n_points);
   //chisq = delta_sum/r2_sum;
   f = chisq;
}

double equation_of_circle(double x, double y, double *par){
  return par[2]*par[2] - (x-par[0])*(x-par[0]) - (y-par[1])*(y-par[1]);
}

double r2_of_circle(double x, double y, double *par){
  return (x-par[0])*(x-par[0]) + (y-par[1])*(y-par[1]);
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
  gMinuit->mnparm(0, "x0", x0in, 0.01, 0,0,ierflg);
  gMinuit->mnparm(1, "y0", y0in, 0.01, 0,0,ierflg);
  gMinuit->mnparm(2, "R", Rin, 0.01, 0,0,ierflg);
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
  //cout<<x0out<<endl
  //   <<y0out<<endl
  //   <<Rout<<endl;
  //
}

void fill_gr_data_gr_weight(TGraph *gr_data, TGraph *gr_weight){
  _gr_data->Clear();
  _gr_weight->Clear();
  //delete _gr_data;
  //delete _gr_weight;
  //_gr_data = new TGraph();
  //_gr_weight = new TGraph();
  //cout<<"_gr_data->GetN()   "<<_gr_data->GetN()<<endl
  //    <<"_gr_weight->GetN() "<<_gr_weight->GetN()<<endl;
  _n_points = 0;
  Double_t x, y, w;
  for(Int_t i = 0;i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    gr_weight->GetPoint(i,w,w);
    _gr_data->SetPoint(i,x,y);
    _gr_weight->SetPoint(i,w,w);
    _n_points++;
  }  
}

