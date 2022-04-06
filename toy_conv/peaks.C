// Illustrates how to find peaks in histograms.
// This script generates a random number of gaussian peaks
// on top of a linear background.
// The position of the peaks is found via TSpectrum and injected
// as initial values of parameters to make a global fit.
// The background is computed and drawn on top of the original histogram.
//
// To execute this example, do
//  root > .x peaks.C  (generate 10 peaks by default)
//  root > .x peaks.C++ (use the compiler)
//  root > .x peaks.C++(30) (generates 30 peaks)
//
// To execute only the first part of the script (without fitting)
// specify a negative value for the number of peaks, eg
//  root > .x peaks.C(-20)
//
//Author: Rene Brun

#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <typeinfo> 
#include "TRandom.h"
//#include "TF1Convolution.h"

Int_t npeaks = 30;
Double_t fpeaks(Double_t *x, Double_t *par) {
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm  = par[3*p+2];
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}
void peaks(Int_t np=1) {
   npeaks = TMath::Abs(np);
   TH1F *h = new TH1F("h","Abs.",500,0,1000);
   Double_t par[3000];
   Int_t p;
   par[3*0+2] = 100;
   par[3*0+3] = 600;
   par[3*0+4] = 20+2;
   TF1 *f = new TF1("f",fpeaks,0,1000,2+3*npeaks);
   f->SetNpx(1000);
   f->SetParameters(par);
   TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
   c1->Divide(2,2);
   c1->cd(1);
   h->FillRandom("f",200000);
   h->Draw();
   h->Fit("gaus");

   TH1F *h2 = new TH1F("h","Empty",500,0,1000);
   TH1F *h5 = new TH1F("h","Empty",500,0,1000);
   Double_t par[3000];
   Int_t p;
   par[3*0+2] = 100;
   par[3*0+3] = 100;
   par[3*0+4] = 10+2;
   TF1 *f2 = new TF1("f2",fpeaks,0,1000,2+3*npeaks);
   f2->SetNpx(1000);
   f2->SetParameters(par);
   c1->cd(2);
   h2->FillRandom("f2",200000);
   int place = 0;
   const Int_t nbins = 500;
   int i = 0;
   while (place == 0){
	   i++;
	   if(int(h2->GetBinContent(i + 1))>int(0.0)) place = i;
   }
   for (i = 0; i < nbins; i++) h5->SetBinContent(i+1,h2->GetBinContent(place+i));
   h5->Draw();
   h5->Fit("gaus");
   
   TH1F *hcon = new TH1F("hconpic","Convolution",500,0,1000);
   for (int i=0;i<1e6;i++) {
	   Double_t x = h->GetRandom();
	   Double_t rand = gRandom->Rndm();
	   //if (rand>0.5) x += h2->GetRandom();
	   //else x -= h2->GetRandom();
	   x += h5->GetRandom()-56;
	   hcon->Fill(x);
   }
   /*
   TF1Convolution *f_conv = new TF1Convolution("gaus","gaus",1,1000,true);
   f_conv->SetRange(1.,1000.);
   f_conv->SetNofPointsFFT(1000);
   TF1   *f3 = new TF1("f",*f_conv, 0., 5., f_conv->GetNpar());
   f3->SetParameters(1.,-0.3,0.,1.);

   f3->FillRandom("f3",200000);
   f3->Draw();
   */
   c1->cd(3);
   hcon->Draw();
   //TH1F *h3 = (TH1F*)h->Clone("h3");
   
   //Use TSpectrum to find the peak candidates
   TH1F *hdecon = new TH1F("h","Gold Decon.",500,0,1000);
   Int_t i;
   Double_t xmin     = 0;
   Double_t xmax     = 2*nbins;
   Double_t source[nbins];
   Double_t response[nbins];
   float* Source_thetaX = new float[nbins];
   float* Response_thetaX = new float[nbins];
   //gROOT->ForceStyle();
   for (i = 0; i < 100; i++) Source_thetaX[i]=0;
   for (i = 100; i < nbins; i++) Source_thetaX[i]=hcon->GetBinContent(i + 1);
   for (i = 0; i < 100; i++) Response_thetaX[i]=h5->GetBinContent(i + 1);
   for (i = 100; i < nbins; i++) Response_thetaX[i]=0;
   TSpectrum *s = new TSpectrum();
   s->Deconvolution(Source_thetaX,Response_thetaX,500,10,1,0);
   for (i = 0; i < nbins; i++) hdecon->SetBinContent(i + 1,Source_thetaX[i]);
   //for (i = 0; i < nbins; i++) hdecon->SetBinContent(i + 1 - 50,Source_thetaX[i]);
 
   if (hdecon) c1->Update();
   if (np <0) return;

   c1->cd(4);
   hdecon->SetLineColor(kRed);
   hdecon->Draw("SAME L");
   hdecon->Chi2Test(h, "WW P");
   std::cout << hdecon->KolmogorovTest(h) << std::endl;
   c1->Print("Goldsimple.pdf");
}
