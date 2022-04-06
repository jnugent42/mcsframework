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

void toy_conv(Int_t np=1) {
    
    TH1F *hconpic = new TH1F("hconpic","Convolution",100,-100,100);
    TH1F *hconpictrunc = new TH1F("hconpic","Convolution",100,-100,100);
    TH1F *hconpicmoretrunc = new TH1F("hconpic","Convolution",100,-100,100);
    TH1F *difftrunc = new TH1F("difftrunc","sqrt(full-trunc)^2/full",100,-100,100);
    TH1F *diffmoretrunc = new TH1F("diffmoretrunc","sqrt(full-moretrunc)^2/full",100,-100,100);
    
    int jj=1000000;
    const Int_t nbins = 500;
    
    TH1F *h = new TH1F("h","Model.",500,0,1000);
    Double_t norm  = 100;
    Double_t mean  = 600;
    Double_t sigma = 21;
    TF1 *f = new TF1("f","gaus(0)+gaus(3)",0,1000);
    f->SetNpx(1000);
    f->SetParameter(0,norm*.98);
    f->SetParameter(1,mean);
    f->SetParameter(2,sigma);
    f->SetParameter(3,norm*.02);
    f->SetParameter(4,mean);
    f->SetParameter(5,sigma*2);
    TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
    c1->Divide(2,2);
    c1->cd(1);
    h->FillRandom("f",jj);
    h->Draw();
    TH1F *hpic = new TH1F("hpic","Model.",100,-100,100);
    for (int i = 1; i < nbins; i++) hpic->SetBinContent(i,h->GetBinContent(i+250));
    hpic->Draw();

    TH1F *h2 = new TH1F("h","Empty",500,0,1000);
    norm  = 100;
    mean  = 100;
    sigma = 12;
    double a =-60;
    double b =60;
    /* TF1* pdf = new TF1("pdf", "ROOT::Math::tdistribution_pdf(x/10,60.0)", a,b); */
    TF1 *f2 = new TF1("f2","gaus(0)+gaus(3)",0,1000);
    /* TF1 *f2 = new TF1("f2","gaus(0)",0,1000); */
    f2->SetNpx(1000);
    f2->SetParameter(0,norm*.98);
    f2->SetParameter(1,mean);
    f2->SetParameter(2,sigma);
    f2->SetParameter(3,norm*.02);
    f2->SetParameter(4,mean);
    f2->SetParameter(5,sigma*2);
    c1->cd(2);
    f2->Print();
    /* pdf->DrawCopy(); */
    h2->FillRandom("f2",jj);
    h2->Draw();
    /* for (i = 21; i < 81; i++) h2->SetBinContent(i,h2->GetBinContent(i)+1000); */
    /* for (i = 1; i < 21; i++) h2->SetBinContent(i,0); */
    /* for (i = 80; i < nbins; i++) h2->SetBinContent(i,0); */
    TH1F *h2pic = new TH1F("hpic","Empty",100,-100,100);
    for (i = 1; i < nbins; i++) h2pic->SetBinContent(i,h2->GetBinContent(i));
    h2pic->Draw();

    /* TH1F *h3 = new TH1F("h","Empty",500,0,1000); */
    /* norm  = 0.01; */
    /* mean  = 600; */
    /* sigma = 24; */
    /* TF1 *f3 = new TF1("f3","gaus(0)",0,1000); */
    /* f3->SetNpx(1000); */
    /* f3->SetParameter(0,norm); */
    /* f3->SetParameter(1,mean); */
    /* f3->SetParameter(2,sigma); */
    /* h3->FillRandom("f3",jj/1090); */
    /* TH1F *h3pic = new TH1F("hpic","Empty",100,-100,100); */
    /* TH1F *hcon4pic = new TH1F("hpic","Empty",100,-100,100); */
    /* for (i = 1; i < nbins; i++) h3pic->SetBinContent(i,h3->GetBinContent(i+250)); */
    /* /1* h3pic->Draw(); *1/ */
    /* TH1F *hcon4 = new TH1F("hconpic","Convolution",500,0,1000); */
    /* for (int i=0;i<1e6;i++) { */
    /*     Double_t x2 = h2->GetRandom(); */
    /*     Double_t x3 = h3->GetRandom(); */
    /*     Double_t x4 = x3 + x2; */
    /*     hcon4->Fill(x4); */
    /* } */
    /* for (i = 1; i < nbins; i++) hcon4pic->SetBinContent(i,hcon4->GetBinContent(i+300)); */
    /* hcon4pic->Draw(); */
    
    TH1F *hcon = new TH1F("hconpic","Convolution",500,0,1000);
    TH1F *hcon2 = new TH1F("hconpic","Convolution",500,0,1000);
    TH1F *hcon3 = new TH1F("hconpic","Convolution",500,0,1000);
    for (int i=0;i<1e6;i++) {
        Double_t x2 = h->GetRandom();
        /* Double_t x3 = pdf->GetRandom(); */
        Double_t x3 = h2->GetRandom();
        Double_t x4 = x3 + x2;
        hcon->Fill(x4);
        /* while (x3>150 || x3<50) { */
        /*     x3 = h2->GetRandom(); */ 
        /* } */
        if (x3<150 && x3>50) {
            x4 = x3 + x2;
            hcon2->Fill(x4);
        }
        /* while (x3>140 || x3<60) { */
        /*     x3 = h2->GetRandom(); */ 
        /* } */
        if (x3<140 && x3>60) {
            x4 = x3 + x2;
            hcon3->Fill(x4);

        }
    }
    c1->cd(3);
    hcon->Draw();
    for (i = 1; i < nbins; i++) hconpic->SetBinContent(i,hcon->GetBinContent(i+300));
    hconpic->Draw();

    c1->Print("fullconv.pdf");

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
    c1->Divide(2,2);
    c1->cd(1);
    hpic->Draw();

    c1->cd(2);
    for (i = 1; i < 26; i++) h2->SetBinContent(i,0);
    for (i = 76; i < nbins; i++) h2->SetBinContent(i,0);
    for (i = 1; i < nbins; i++) h2pic->SetBinContent(i,h2->GetBinContent(i));
    h2pic->Draw();

    c1->cd(3);
    hcon2->Draw();
    for (i = 1; i < nbins; i++) hconpictrunc->SetBinContent(i,hcon2->GetBinContent(i+300));
    hconpictrunc->Draw();

    c1->cd(4);
    for (i = 1; i < nbins; i++) difftrunc->SetBinContent(i,100*sqrt(pow(hconpic->GetBinContent(i)-hconpictrunc->GetBinContent(i),2))/hconpic->GetBinContent(i));
    difftrunc->GetYaxis()->SetRangeUser(0,20);
    difftrunc->GetXaxis()->SetRangeUser(-40,40);
    difftrunc->Draw();
    c1->Print("truncatedconv.pdf");

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
    c1->Divide(2,2);
    c1->cd(1);
    hpic->Draw();

    c1->cd(2);
    for (i = 1; i < 31; i++) h2->SetBinContent(i,0);
    for (i = 71; i < nbins; i++) h2->SetBinContent(i,0);
    for (i = 1; i < nbins; i++) h2pic->SetBinContent(i,h2->GetBinContent(i));
    h2pic->Draw();

    c1->cd(3);
    hcon3->Draw();
    for (i = 1; i < nbins; i++) hconpicmoretrunc->SetBinContent(i,hcon3->GetBinContent(i+300));
    hconpicmoretrunc->Draw();

    c1->cd(4);
    for (i = 1; i < nbins; i++) diffmoretrunc->SetBinContent(i,100*sqrt(pow(hconpic->GetBinContent(i)-hconpicmoretrunc->GetBinContent(i),2))/hconpic->GetBinContent(i));
    diffmoretrunc->GetYaxis()->SetRangeUser(0,20);
    diffmoretrunc->GetXaxis()->SetRangeUser(-40,40);
    diffmoretrunc->Draw();
    c1->Print("moretruncatedconv.pdf");
}
