import libxml2
import sys, os, subprocess, math, ROOT
from ROOT import TH1D, TCanvas, TFile, TGraphErrors, TMarker, TF1, TLegend, TText, TLine, TROOT
from pprint import pprint
import re

TROOT.gROOT.SetBatch(1)
def plotFiles(syslist):

    hists  = []
    Graphs = []
    pointlist = []

    for sysfile in syslist:
        sysval = ExtractPars(sysfile)
        pointlist.append(sysval)

    pprint(pointlist)

    Graph_chi2_truthgold = TGraphErrors(20)
    for j in range(2,len(pointlist)):

	    Graph_chi2_truthgold.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][0]))
    j = 0
    c = TCanvas()
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graph_chi2_truthgold.SetMarkerStyle(20)
    Graph_chi2_truthgold.Draw('AEP')
    line3 = TLine()
    line3.DrawLine(-0.012,0,0.0115,0);

    Graph_chi2_truthgold.GetHistogram().SetTitle("")
    Graph_chi2_truthgold.GetHistogram().GetXaxis().SetTitle("Eff cut (mrad)")
    Graph_chi2_truthgold.GetHistogram().GetYaxis().SetTitle("#chi^2/NDF ")
    Graph_chi2_truthgold.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graph_chi2_truthgold.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graph_chi2_truthgold.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graph_chi2_truthgold.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graph_chi2_truthgold.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg = TLegend(0.5,0.5,0.89,0.89)
    leg.SetLineColor(10)
    leg.AddEntry(Graph_chi2_truthgold,"Data","p")

    j+=1
    t1 = TText(0.68,0.215,"MICE Internal")
    t2 = TText(0.68,0.185,"ISIS Cycle 2015/04")
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graph_chi2_truthgold.pdf")
    c.Clear()

    Graph_Kol_truthgold = TGraphErrors(20)
    for j in range(2,len(pointlist)):

	    Graph_Kol_truthgold.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][2]))
    j = 0
    c = TCanvas()
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graph_Kol_truthgold.SetMarkerStyle(20)
    Graph_Kol_truthgold.Draw('AEP')
    line3 = TLine()
    line3.DrawLine(-0.012,0,0.0115,0);

    Graph_Kol_truthgold.GetHistogram().SetTitle("")
    Graph_Kol_truthgold.GetHistogram().GetXaxis().SetTitle("Eff cut (mrad)")
    Graph_Kol_truthgold.GetHistogram().GetYaxis().SetTitle("#chi^2/NDF ")
    Graph_Kol_truthgold.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graph_Kol_truthgold.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graph_Kol_truthgold.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graph_Kol_truthgold.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graph_Kol_truthgold.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graph_Kol_truthgold,"Data","p")

    j+=1
    t1 = TText(0.68,0.215,"MICE Internal")
    t2 = TText(0.68,0.185,"ISIS Cycle 2015/04")
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graph_Kol_truthgold.pdf")
    c.Clear()

def ExtractPars(xmlfile):
    print xmlfile
    doc = libxml2.parseFile(xmlfile)
    vals = []
    rootfile = next(x.prop("name") for x in doc.xpathEval("spec/file") \
                        if x.prop("id").find("outfile") >= 0)
    print rootfile
    rf = TFile(rootfile)
    if not rf:
        return vals
    if rf.IsZombie():
        return vals
    print rootfile
    eff_cut = next(x.prop("name") for x in doc.xpathEval("spec/file") \
                        if x.prop("id").find("trkreffiname") >= 0)
    thetaY_graph = rf.Get("thetaY_graph")
    thetaY_recoGold = rf.Get("thetaY_recoGold_effi_only_gold")
    thetaY_graph.GetXaxis().SetRangeUser(-0.04,0.04)
    thetaY_recoGold.GetXaxis().SetRangeUser(-0.04,0.04)
    chi2=thetaY_graph.Chi2Test(thetaY_recoGold,"CHI2/NDF")
    vals.append(chi2)
    result = re.search('.*6_(.*).root', eff_cut)
    print(result.group(1))
    vals.append(result.group(1))

    kolmogorov=thetaY_graph.KolmogorovTest(thetaY_recoGold, "M")
    vals.append(kolmogorov)

    return vals

if __name__=="__main__":
    plotFiles(sys.argv[1:])
