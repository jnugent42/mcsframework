import libxml2
import sys, os, subprocess, math, ROOT
from ROOT import TH1D, TCanvas, TFile, TGraphErrors, TMarker, TF1, TLegend, TText, TLine
from pprint import pprint
import ROOT
from ROOT import TH1D, TH2D, TCanvas, TLegend, TFile, TText, TLine, TArrow, TROOT

def plotFiles(syslist):

    hists  = []
    Graphs = []
    pointlist = []

    for sysfile in syslist:
        sysval = ExtractPars(sysfile)
        pointlist.append(sysval)

    pprint(pointlist)

    Graph = TGraphErrors(100)
    for j in range(0,20):
        print float(pointlist[j][0])
        print float(pointlist[j][1])
        Graph.SetPoint(j,float(pointlist[j][1]),float(pointlist[j][0]))
    j = 0
    c = TCanvas()
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graph.SetMarkerStyle(20)
    #Graph.SetMarkerSize(4)
    Graph.Draw('AEP')
    line3 = TLine()
    line3.DrawLine(-0.012,0,0.0115,0);

    Graph.GetHistogram().SetTitle("")
    Graph.GetHistogram().GetXaxis().SetTitle("Offset")
    Graph.GetHistogram().GetYaxis().SetTitle("Mean")
    Graph.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graph.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graph.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graph.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graph.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    if Graph.GetName() == "integralX":
        leg = TLegend(0.5,0.15,0.79,0.39)
    else:
        leg = TLegend(0.5,0.5,0.89,0.89)
    leg.SetLineColor(10)
    leg.AddEntry(Graph,"Data","p")

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
    c.Print("Graph.pdf")
    c.Clear()

    #Graphempty = TGraphErrors(100)
    #for j in range(0,100):
	    #print 'Graphempty'
	    #print pointlist[j-1][1],pointlist[j-1][25]
	    ##Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    #Graphempty.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][2]))
	    #Graphempty.SetPointError(j,0,float(pointlist[j-1][37]))
    #j = 0
    #c.SetLeftMargin(0.125)
    #c.SetBottomMargin(0.125)
    #Graphempty.SetMarkerStyle(20)
    #Graphempty.Draw('AEP')
    #line3.DrawLine(-0.012,0,0.0115,0);

    #Graphempty.GetHistogram().SetTitle("")
    #Graphempty.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks in dXdz (rad)")
    #Graphempty.GetHistogram().GetYaxis().SetTitle("Gradient of asymmetry plot")
    #Graphempty.GetHistogram().GetXaxis().SetLabelSize(0.05)
    #Graphempty.GetHistogram().GetYaxis().SetLabelSize(0.05)
    #Graphempty.GetHistogram().GetXaxis().SetTitleSize(0.05)
    #Graphempty.GetHistogram().GetYaxis().SetTitleSize(0.05)
    #Graphempty.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    #leg.SetLineColor(10)
    #leg.AddEntry(Graphempty,"Data","p")

    #j+=1
    #t1.SetNDC(1)
    #t1.SetTextSize(0.04)
    #t1.SetTextFont(42)
    #t2.SetNDC(1)
    #t2.SetTextSize(0.03)
    #t2.SetTextFont(42)
    #t1.Draw()
    #t2.Draw()
    #c.Print("Graphempty.pdf")
    #c.Clear()

def ExtractPars(xmlfile):

    print xmlfile
    doc = libxml2.parseFile(xmlfile)
    vals = []
    rootfile = next(x.prop("name") for x in doc.xpathEval("spec/file") \
                        if x.prop("id").find("outfile") >= 0)
    rf = TFile(rootfile)
    if not rf:
        return vals
    if rf.IsZombie():
        return vals
    print rootfile
    offset = next(x.prop("value") for x in doc.xpathEval("spec/cuts") \
                        if x.prop("name").find("offset") >= 0)
    print offset

    theta = rf.Get("thetaX_recoGold")
    vals.append(theta.GetMean(1))
    vals.append(offset)
    meanabsLiH = rf.Get("thetaX_Goldmirror")

    if (int(offset)==0):
        c = TCanvas()
        hgoldmirrorl = rf.Get("hgoldmirrorl")
        hgoldmirrorr = rf.Get("hgoldmirrorr")
        hgoldmirrorl.SetStats(0)
        hgoldmirrorl.SetTitle("")
        hgoldmirrorl.SetMinimum(0)
        hgoldmirrorl.GetXaxis().SetLabelSize(0.05)
        hgoldmirrorl.GetXaxis().SetTitleSize(0.05)
        hgoldmirrorl.GetYaxis().SetLabelSize(0.05)
        hgoldmirrorl.GetYaxis().SetTitleSize(0.05)
        hgoldmirrorl.GetXaxis().SetNdivisions(508)
        hgoldmirrorl.GetYaxis().SetNdivisions(508)
        hgoldmirrorr.SetLineColor(2)
        hgoldmirrorl.Draw()
        hgoldmirrorr.Draw("SAME")
        # t1.Draw()
        c.Print("mirrorgold.eps")
        c.Print("mirrorgold_pq.jpg")

    return vals


if __name__=="__main__":
    plotFiles(sys.argv[1:])
