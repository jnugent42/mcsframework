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
void binsizepeaks(Int_t np=1) {
   Double_t *x = 0;
   TGraph *mirrorstat = new TGraph();
	   int zz = 0;
   for (int jj=23; jj<24; jj=jj+2){ 
	   std::cout << jj << std::endl;
   npeaks = TMath::Abs(np);
   const Int_t nbins = jj*5;
   const Int_t bins = 1000/nbins;
   TH1F *h = new TH1F("h","Abs.",nbins,0,1000);
   Int_t p;
   Double_t norm  = 100;
   Double_t mean  = 500;
   Double_t sigma = 30+2;
   TF1 *f = new TF1("f","gaus(0)",0,1000);
   f->SetNpx(1000);
   f->SetParameter(0,norm);
   f->SetParameter(1,mean);
   f->SetParameter(2,sigma);
   TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
   TCanvas *c3 = new TCanvas("c3","c3",10,10,1000,900);
   c3->Divide(2,2);
   c1->Divide(2,2);
   c1->cd(1);
   gPad->SetLogy();
   h->FillRandom("f",200000);
   hn = h->GetEntries();
   h->Scale(1./hn);
   h->Draw();
   TH1F *hpic = new TH1F("hpic","Abs.",jj,-100,100);
   int meanp = 0;
   int i=0;
   while (meanp==0) {
	   if (h->GetBinLowEdge(i)>h->GetMean()) meanp=i;
	   i++;
   }
   std::cout << "meanp " << meanp << std::endl;
   std::cout << "calc " << meanp-jj/2 << std::endl;
   std::cout << "mean " << h->GetMean() << std::endl;
   std::cout << "edge " << h->GetBinLowEdge(70) << std::endl;
   std::cout << "edge " << h->GetBinLowEdge(71) << std::endl;
   /* for (i = (jj+1)/2; i < jj+1; i++) hpic->SetBinContent(i,h->GetBinContent(i-((jj+1)/2)+meanp-1)); */
   /* for (i = 1; i < (jj+1)/2; i++) hpic->SetBinContent(i,h->GetBinContent(meanp-((jj+1)/2)+i-1)); */
   for (i = 1; i < (jj+1); i++) hpic->SetBinContent(i,h->GetBinContent(i-1-(jj+1)/2+meanp));
   hpic->Draw();
   int mirronbins = jj;
   TH1F *hmirrorl = new TH1F("hmirrorl","Mirror",jj,-100,100);
   TH1F *hmirrorr = new TH1F("hmirrorr","Mirror",jj,-100,100);
   for (i = 1; i < mirronbins/2+1; i++) hmirrorl->SetBinContent(i,hpic->GetBinContent(i));
   for (i = 1; i < mirronbins/2+1; i++) hmirrorr->SetBinContent(i,hpic->GetBinContent(mirronbins+1-i));
   c3->cd(1);
   hmirrorl->Draw();
   hmirrorr->SetLineColor(kRed);
   hmirrorr->Draw("SAME");
   
   TH1F *h2 = new TH1F("h","Empty",nbins,0,1000);
   TH1F *h5 = new TH1F("h","test",nbins,0,1000);
   Int_t p;
   norm  = 100;
   mean  = 100;
   sigma = 15+2;
   TF1 *f2 = new TF1("f2","gaus(0)",0,1000);
   f2->SetNpx(1000);
   f2->SetParameter(0,norm);
   f2->SetParameter(1,mean);
   f2->SetParameter(2,sigma);
   c1->cd(2);
   gPad->SetLogy();
   h2->FillRandom("f2",200000);
   //h2n = h2->GetEntries();
   //h2->Scale(1./h2n);
   h2->Draw();
   int place = 0;
   int i = 0;
   while (place == 0){
	   i++;
	   if(int(h2->GetBinContent(i + 1))>int(0.0)) place = i;
   }
   for (i = 0; i < nbins; i++) h5->SetBinContent(i+1,h2->GetBinContent(place+i));
   h5->Draw();
   TH1F *h2pic = new TH1F("h2pic","Empty",jj,-100,100);
   for (i = 1; i < nbins; i++) h2pic->SetBinContent(i,h2->GetBinContent(i));
   h2pic->Draw();
   TH1F *h2mirrorl = new TH1F("h2mirrorl","Mirror",jj,-100,100);
   TH1F *h2mirrorr = new TH1F("h2mirrorr","Mirror",jj,-100,100);
   for (i = 1; i < mirronbins/2+1; i++) h2mirrorl->SetBinContent(i,h2pic->GetBinContent(i));
   for (i = 1; i < mirronbins/2+1; i++) h2mirrorr->SetBinContent(i,h2pic->GetBinContent(mirronbins+1-i));
   c3->cd(2);
   h2mirrorl->Draw();
   h2mirrorr->SetLineColor(kRed);
   h2mirrorr->Draw("SAME");
   
   TH1F *hcon = new TH1F("hconpic","Convolution",nbins,0,1000);
   for (int i=0;i<1e6;i++) {
       Double_t x2 = h->GetRandom();
       x2 += h2->GetRandom()-100;
       hcon->Fill(x2);
   }

   c1->cd(3);
   gPad->SetLogy();
   //hconn = hcon.GetEntries();
   //hcon->Scale(1./hconn);
   hcon->Draw();
   TH1F *hconpic = new TH1F("hconpic","Convolution",jj,-100,100);
   meanp = 0;
   i=0;
   while (meanp==0) {
	   if (hcon->GetBinLowEdge(i)>hcon->GetMean()) meanp=i;
	   i++;
   }
   for (i = (jj+1)/2; i < jj+1; i++) hconpic->SetBinContent(i,hcon->GetBinContent(i-((jj+1)/2)+meanp-1));
   for (i = 1; i < (jj+1)/2; i++) hconpic->SetBinContent(i,hcon->GetBinContent(meanp-((jj+1)/2)+i-1));
   hconpic->Draw();
   TH1F *hconmirrorl = new TH1F("hconmirrorl","Mirror",jj,-100,100);
   TH1F *hconmirrorr = new TH1F("hconmirrorr","Mirror",jj,-100,100);
   for (i = 1; i < mirronbins/2+1; i++) hconmirrorl->SetBinContent(i,hconpic->GetBinContent(i));
   for (i = 1; i < mirronbins/2+1; i++) hconmirrorr->SetBinContent(i,hconpic->GetBinContent(mirronbins+1-i));
   c3->cd(3);
   hconmirrorl->Draw();
   hconmirrorr->SetLineColor(kRed);
   hconmirrorr->Draw("SAME");
   
   TH1F *hdecon = new TH1F("h","Gold Decon.",nbins,0,1000);
   Int_t i;
   Double_t xmin     = 0;
   Double_t xmax     = 2*nbins;
   Double_t source[nbins];
   Double_t response[nbins];
   float* Source_thetaX = new float[nbins];
   float* Response_thetaX = new float[nbins];
   for (i = 0; i < nbins; i++) Source_thetaX[i]=hcon->GetBinContent(i + 1);
   for (i = 0; i < nbins; i++) Response_thetaX[i]=h2->GetBinContent(i + 1);
   TSpectrum *s = new TSpectrum();
   s->DeconvolutionRL(Source_thetaX,Response_thetaX,nbins,100,1,0);
   for (i = 0; i < nbins; i++) hdecon->SetBinContent(i + 1 - jj/2,Source_thetaX[i]);
 
   if (np <0) return;

   c1->cd(4);
   gPad->SetLogy();
   hdecon->SetLineColor(kRed);
   //hdeconn = hdecon.GetEntries();
   //hdecon->Scale(1./hdeconn);
   hdecon->Draw("SAME L");
   TH1F *hdeconpic = new TH1F("hdeconpic","Gold Decon.",jj,-100,100);
   meanp = 0;
   i=0;
   while (meanp==0) {
	   if (hdecon->GetBinLowEdge(i)>hdecon->GetMean()) meanp=i;
	   i++;
   }
   std::cout << "meanp " << meanp << std::endl;
   //meanp = hdecon->GetMaximumBin();
   for (i = (jj+1)/2; i < jj+1; i++) hdeconpic->SetBinContent(i,hdecon->GetBinContent(i-((jj+1)/2)+meanp-1));
   for (i = 1; i < (jj+1)/2; i++) hdeconpic->SetBinContent(i,hdecon->GetBinContent(meanp-((jj+1)/2)+i-1));
   /* for (i = 12; i < jj+1; i++) hdeconpic->SetBinContent(i,hdecon->GetBinContent(i-((jj+1)/2)+meanp)); */
   /* for (i = 1; i < 12; i++) hdeconpic->SetBinContent(i,hdecon->GetBinContent(meanp-((jj+1)/2)+i-1)); */
   /* for (i = 1; i < nbins; i++) hdeconpic->SetBinContent(i,hdecon->GetBinContent(i+meanp-jj/2)); */
   //for (i = 1; i < nbins; i++) hdeconpic->SetBinContent(i,hdecon->GetBinContent(i+250));
   hdeconpic->SetLineColor(kRed);
   hdeconpic->DrawNormalized();
   hpic->DrawNormalized("SAME");
   c1->SetLogy();
   float residual = 0;
   /* mirronbins = 101; */
   TH1F *hdeconmirrorl = new TH1F("hmirrorl","Mirror",jj,-100,100);
   TH1F *hdeconmirrorr = new TH1F("hmirrorr","Mirror",jj,-100,100);
   for (i = 1; i < mirronbins/2+1; i++) hdeconmirrorl->SetBinContent(i,hdeconpic->GetBinContent(i));
   for (i = 1; i < mirronbins/2+1; i++) hdeconmirrorr->SetBinContent(i,hdeconpic->GetBinContent(mirronbins+1-i));
   c3->cd(4);
   hdeconmirrorl->SetMarkerStyle(20);
   hdeconmirrorl->SetMarkerSize(1);
   hdeconmirrorl->Draw();
   hdeconmirrorr->SetMarkerStyle(20);
   hdeconmirrorr->SetMarkerSize(1);
   hdeconmirrorr->SetLineColor(kRed);
   hdeconmirrorr->Draw("SAME");
   for (i = 1; i < mirronbins/2; i++) {
	   Double_t xparl;
	   Double_t yparl;
	   Double_t xparr;
	   Double_t yparr;
	   yparl = hdeconmirrorl->GetBinContent(i);
	   yparr = hdeconmirrorr->GetBinContent(i);
	   residual += sqrt(pow(yparl-yparr,2));
   }
   //if (residual>0.14e-3 && residual<0.16e-3 && jj>60 ) break;
   //if (jj>20 && hdeconpic->GetMean()>1 || hdeconpic->GetMean()<-1 ) continue;
   mirrorstat->SetPoint(zz,jj,hdeconpic->GetMean());
   residual = 0;
   zz += 1;
   hdecon->Chi2Test(h, "P");
   std::cout << hdecon->KolmogorovTest(h, "D") << std::endl;
   }
   TCanvas *c2 = new TCanvas("c2","c2",10,10,1000,900);
   mirrorstat->SetMarkerStyle(20); 
   mirrorstat->SetMarkerSize(1); 
   mirrorstat->Draw("AP"); 
   c1->Print("Goldcmean.pdf");
   c3->Print("c3Goldcmean.pdf");
}
