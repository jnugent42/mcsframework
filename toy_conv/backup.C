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
    TGraph *mirrorstat = new TGraph();
    int zz = 0;
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
    Double_t sigma = 20;
    TF1 *f = new TF1("f","gaus(0)",0,1000);
    f->SetNpx(1000);
    f->SetParameter(0,norm);
    f->SetParameter(1,mean);
    f->SetParameter(2,sigma);
    TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
    c1->Divide(2,2);
    c1->cd(1);
    h->FillRandom("f",jj);
    /* for (i = 1; i < nbins; i++) h->SetBinContent(i,h->GetBinContent(i)+1e3); */
    h->Draw();
    TH1F *hpic = new TH1F("hpic","Model.",100,-100,100);
    for (int i = 1; i < nbins; i++) hpic->SetBinContent(i,h->GetBinContent(i+250));
    hpic->Draw();

    TH1F *h2 = new TH1F("h","Empty",500,0,1000);
    norm  = 100;
    mean  = 100;
    sigma = 12;
    TF1 *f2 = new TF1("f2","gaus(0)",0,1000);
    f2->SetNpx(1000);
    f2->SetParameter(0,norm);
    f2->SetParameter(1,mean);
    f2->SetParameter(2,sigma);
    c1->cd(2);
    h2->FillRandom("f2",jj);
    /* for (i = 21; i < 81; i++) h2->SetBinContent(i,h2->GetBinContent(i)+1e3); */
    TH1F *h2pic = new TH1F("hpic","Empty",100,-100,100);
    for (i = 1; i < nbins; i++) h2pic->SetBinContent(i,h2->GetBinContent(i));
    h2pic->Draw();

    TH1F *hcon = new TH1F("hconpic","Convolution",500,0,1000);
    TH1F *hcon2 = new TH1F("hconpic","Convolution",500,0,1000);
    for (int i=0;i<1e6;i++) {
        Double_t x2 = h->GetRandom();
        Double_t x3 = h2->GetRandom();
        x2 += x3;
        hcon->Fill(x2);
        if (x3<160 && x3>60) hcon2->Fill(x2);
    }
    c1->cd(3);
    hcon->Draw();
    for (i = 1; i < nbins; i++) hconpic->SetBinContent(i,hcon->GetBinContent(i+300));
    hconpic->Draw();

    c1->Print("fullconv.pdf");

    /* int jj=2000000; */
    /* const Int_t nbins = 500; */
    /* TH1F *h = new TH1F("h","Model.",500,0,1000); */
    /* Double_t norm  = 100; */
    /* Double_t mean  = 600; */
    /* Double_t sigma = 20; */
    /* TF1 *f = new TF1("f","gaus(0)",0,1000); */
    /* f->SetNpx(1000); */
    /* f->SetParameter(0,norm); */
    /* f->SetParameter(1,mean); */
    /* f->SetParameter(2,sigma); */
    TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
    c1->Divide(2,2);
    c1->cd(1);
    /* h->FillRandom("f",jj); */
    /* /1* for (i = 1; i < nbins; i++) h->SetBinContent(i,h->GetBinContent(i)+1e3); *1/ */
    /* h->Draw(); */
    /* TH1F *hpic = new TH1F("hpic","Model.",100,-100,100); */
    /* for (int i = 1; i < nbins; i++) hpic->SetBinContent(i,h->GetBinContent(i+250)); */
    hpic->Draw();

    /* TH1F *h2 = new TH1F("h","Empty",500,0,1000); */
    /* norm  = 100; */
    /* mean  = 100; */
    /* sigma = 12; */
    /* TF1 *f2 = new TF1("f2","gaus(0)",0,1000); */
    /* f2->SetNpx(1000); */
    /* f2->SetParameter(0,norm); */
    /* f2->SetParameter(1,mean); */
    /* f2->SetParameter(2,sigma); */
    c1->cd(2);
    /* h2->FillRandom("f2",jj); */
    for (i = 1; i < 12; i++) h2->SetBinContent(i,0);
    for (i = 86; i < nbins; i++) h2->SetBinContent(i,0);
    /* TH1F *h2pic = new TH1F("hpic","Empty",100,-100,100); */
    h2pic->Clear();
    for (i = 1; i < nbins; i++) h2pic->SetBinContent(i,h2->GetBinContent(i));
    /* for (i = 1; i < 12; i++) h2pic->SetBinContent(i,0); */
    /* for (i = 86; i < nbins; i++) h2pic->SetBinContent(i,0); */
    h2pic->Draw();


    /* TH1F *hcon = new TH1F("hconpic","Convolution",500,0,1000); */
    /* hcon->Clear(); */
    /* for (i = 1; i < nbins; i++) hcon->SetBinContent(i,0); */
    /* for (int i=0;i<1e6;i++) { */
    /*     Double_t x2 = h->GetRandom(); */
    /*     x2 += h2->GetRandom(); */
    /*     hcon->Fill(x2); */
    /* } */
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

    /* int jj=2000000; */
    /* const Int_t nbins = 500; */
    /* TH1F *h = new TH1F("h","Model.",500,0,1000); */
    /* Double_t norm  = 100; */
    /* Double_t mean  = 600; */
    /* Double_t sigma = 20; */
    /* TF1 *f = new TF1("f","gaus(0)",0,1000); */
    /* f->SetNpx(1000); */
    /* f->SetParameter(0,norm); */
    /* f->SetParameter(1,mean); */
    /* f->SetParameter(2,sigma); */
    /* TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900); */
    /* c1->Divide(2,2); */
    /* c1->cd(1); */
    /* h->FillRandom("f",jj); */
    /* h->Draw(); */
    /* TH1F *hpic = new TH1F("hpic","Model.",100,-100,100); */
    /* for (int i = 1; i < nbins; i++) hpic->SetBinContent(i,h->GetBinContent(i+250)); */
    /* hpic->Draw(); */

    /* TH1F *h2 = new TH1F("h","Empty",500,0,1000); */
    /* norm  = 100; */
    /* mean  = 100; */
    /* sigma = 12; */
    /* TF1 *f2 = new TF1("f2","gaus(0)",0,1000); */
    /* f2->SetNpx(1000); */
    /* f2->SetParameter(0,norm); */
    /* f2->SetParameter(1,mean); */
    /* f2->SetParameter(2,sigma); */
    /* c1->cd(2); */
    /* h2->FillRandom("f2",jj); */
    /* TH1F *h2pic = new TH1F("hpic","Empty",100,-100,100); */
    /* for (i = 1; i < nbins; i++) h2pic->SetBinContent(i,h2->GetBinContent(i)); */
    /* for (i = 1; i < 41; i++) h2pic->SetBinContent(i,0); */
    /* for (i = 61; i < nbins; i++) h2pic->SetBinContent(i,0); */
    /* h2pic->Draw(); */

    /* TH1F *hcon = new TH1F("hconpic","Convolution",500,0,1000); */
    /* for (int i=0;i<1e6;i++) { */
    /*     Double_t x2 = h->GetRandom(); */
    /*     x2 += h2->GetRandom(); */
    /*     hcon->Fill(x2); */
    /* } */
    /* c1->cd(3); */
    /* hcon->Draw(); */
    /* for (i = 1; i < nbins; i++) hconpicmoretrunc->SetBinContent(i,hcon->GetBinContent(i+300)); */
    /* hconpicmoretrunc->Draw(); */

    /* Int_t mirronbins = 101; */
    /* c1->cd(4); */
    /* for (i = 1; i < nbins; i++) diffmoretrunc->SetBinContent(i,100*sqrt(pow(hconpic->GetBinContent(i)-hconpicmoretrunc->GetBinContent(i),2))/hconpic->GetBinContent(i)); */
    /* diffmoretrunc->GetYaxis()->SetRangeUser(0,20); */
    /* diffmoretrunc->GetXaxis()->SetRangeUser(-40,40); */
    /* diffmoretrunc->Draw(); */
    /* c1->Print("moretruncatedconv.pdf"); */
}
