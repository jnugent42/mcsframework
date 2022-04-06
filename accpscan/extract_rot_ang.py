import libxml2
import sys, os, subprocess, math, ROOT
from ROOT import TH1D, TCanvas, TFile, TGraphErrors, TMarker, TF1, TLegend, TText, TLine


def plotFiles(syslist):
    
    hists  = []
    Graphs = []
    pointlist = []
    
    for sysfile in syslist:
        sysval = ExtractPars(sysfile)
        pointlist.append(sysval)
    
    Graph = TGraphErrors(20)
    for j in range(0,20):
	    print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j][0]))
	    Graph.SetPoint(j,float(pointlist[j][1]),float(pointlist[j][0]))
    j = 0
    c = TCanvas()
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graph.SetMarkerStyle(20)
    #Graph.SetMarkerSize(4)
    Graph.Draw('ap')
    line3 = TLine()
    line3.DrawLine(-0.012,0,0.0115,0);

    Graph.GetHistogram().SetTitle("")
    Graph.GetHistogram().GetXaxis().SetTitle("Acceptance cutoff (mrad)")
    Graph.GetHistogram().GetYaxis().SetTitle("Sigma of Gaussian fitter to central 80 mrad (mrad)")
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
    accpt = next(x.prop("value") for x in doc.xpathEval("spec/sys") \
                        if x.prop("name").find("accpt") >= 0)
    print accpt
    thX_rc = rf.Get("thetaX_recoCobb")
    thX_rc.Fit("gaus","Q","",-0.04,0.04)
    vals.append([1000*thX_rc.GetFunction("gaus").GetParameter("Sigma"), \
               1000*thX_rc.GetFunction("gaus").GetParameter("Sigma")/math.sqrt(thX_md.Integral())])
    return vals

if __name__=="__main__":
    plotFiles(sys.argv[4:])
