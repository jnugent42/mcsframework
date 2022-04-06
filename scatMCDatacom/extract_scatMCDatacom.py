import libxml2
import sys, os, subprocess, math, ROOT
from ROOT import TH1D, TCanvas, TFile, TGraphErrors, TMarker, TF1, TLegend, TText, TLine


def plotFiles(syslist):
    
    for sysfile in syslist:
        sysval = ExtractPars(sysfile)
    
def ExtractPars(xmlfile):
    print xmlfile
    doc = libxml2.parseFile(xmlfile)
    rootfile = next(x.prop("name") for x in doc.xpathEval("spec/file") \
                        if x.prop("id").find("outfile") >= 0)
    print rootfile
    rf = TFile(rootfile)
    if not rf:
        return vals
    if rf.IsZombie():
        return vals
    print rootfile
    thX_mg = rf.Get("thetaX_measdataGEANT")
    thX_mgE = rf.Get("thetaX_ref")
    
    print rootfile
    rf = TFile(rootfile)
    if not rf:
        return vals
    if rf.IsZombie():
        return vals
    print rootfile
    thX_MCmg = rf.Get("thetaX_measdataGEANT")
    thX_MCmgE = rf.Get("thetaX_ref")

    t101 = TText(0.18,0.885,"MICE Internal ISIS cycle 2015/04")
    t201 = TText(0.18,0.85,"LiH, MAUS v3.1.2")
    t101.SetNDC(1)
    t101.SetTextSize(0.04)
    t101.SetTextFont(42)
    t201.SetNDC(1)
    t201.SetTextSize(0.03)
    t201.SetTextFont(42)
    leg.AddEntry(thX_mg,"Data","p")
    leg.AddEntry(th_X_MCmg,"MC","l")
    thX_mg.SetStats(0)
    thX_mg.SetTitle("")
    thX_mg.SetMarkerStyle(20)
    th_X_MCmg.SetFillColor(ROOT.kOrange-2)
    #th_X_MCmg.SetLineStyle(2)
    th_X_MCmg.SetStats(0)
    th_X_MCmg.SetTitle("")
    th_X_MCmg.GetXaxis().SetLabelSize(0.05)
    th_X_MCmg.GetXaxis().SetTitleSize(0.05)
    th_X_MCmg.GetYaxis().SetLabelSize(0.05)
    th_X_MCmg.GetYaxis().SetTitle("No. of events")
    th_X_MCmg.GetYaxis().SetTitleSize(0.05)
    #th_X_MCmg.Sumw2()
    thX_mg.Sumw2()
    n21 = thX_mg.GetEntries()
    thX_mg.Scale(1./n21)
    l21 = th_X_MCmg.GetEntries()
    th_X_MCmg.Scale(1./l21)
    thX_mg.GetYaxis().SetRangeUser(0,0.08)
    th_X_MCmg.GetYaxis().SetRangeUser(0,0.08)
    th_X_MCmg.Draw()
    thX_mg.Draw("e1 p same")
    leg.Draw("same")
    t101.Draw("same")
    t201.Draw("same")
    c.Print(rootfile[:5]+"MCDataScatDis.eps")
    
    t101 = TText(0.18,0.885,"MICE Internal ISIS cycle 2015/04")
    t201 = TText(0.18,0.85,"LiH, MAUS v3.1.2")
    t101.SetNDC(1)
    t101.SetTextSize(0.04)
    t101.SetTextFont(42)
    t201.SetNDC(1)
    t201.SetTextSize(0.03)
    t201.SetTextFont(42)
    leg.AddEntry(thX_MCmgE,"Data","p")
    leg.AddEntry(th_X_MCmg,"MC","l")
    thX_MCmgE.SetStats(0)
    thX_MCmgE.SetTitle("")
    thX_MCmgE.SetMarkerStyle(20)
    th_X_MCmg.SetFillColor(ROOT.kOrange-2)
    #th_X_MCmg.SetLineStyle(2)
    th_X_MCmg.SetStats(0)
    th_X_MCmg.SetTitle("")
    th_X_MCmg.GetXaxis().SetLabelSize(0.05)
    th_X_MCmg.GetXaxis().SetTitleSize(0.05)
    th_X_MCmg.GetYaxis().SetLabelSize(0.05)
    th_X_MCmg.GetYaxis().SetTitle("No. of events")
    th_X_MCmg.GetYaxis().SetTitleSize(0.05)
    #th_X_MCmg.Sumw2()
    thX_MCmgE.Sumw2()
    n21 = thX_MCmgE.GetEntries()
    thX_MCmgE.Scale(1./n21)
    l21 = th_X_MCmg.GetEntries()
    th_X_MCmg.Scale(1./l21)
    thX_MCmgE.GetYaxis().SetRangeUser(0,0.08)
    th_X_MCmg.GetYaxis().SetRangeUser(0,0.08)
    th_X_MCmg.Draw()
    thX_MCmgE.Draw("e1 p same")
    leg.Draw("same")
    t101.Draw("same")
    t201.Draw("same")
    c.Print(rootfile[:5]+"MCDataScatDisE.eps")

if __name__=="__main__":
    plotFiles(sys.argv[4:])
