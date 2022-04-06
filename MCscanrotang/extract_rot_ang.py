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
    
    Graph = TGraphErrors(100)
    for j in range(0,100):
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
    Graph.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks in dXdz (mrad)")
    Graph.GetHistogram().GetYaxis().SetTitle("Gradient of asymmetry plot")
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

    Graphempty = TGraphErrors(100)
    for j in range(0,100):
	    print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j][0]))
	    Graphempty.SetPoint(j,float(pointlist[j][1]),float(pointlist[j][2]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphempty.SetMarkerStyle(20)
    Graphempty.Draw('ap')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphempty.GetHistogram().SetTitle("")
    Graphempty.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks in dXdz (mrad)")
    Graphempty.GetHistogram().GetYaxis().SetTitle("Gradient of asymmetry plot")
    Graphempty.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphempty.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphempty.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphempty.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphempty.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphempty,"Data","p")
    
    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphempty.pdf")
    c.Clear()

    GraphemptyX = TGraphErrors(100)
    for j in range(0,100):
	    print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j][0]))
	    GraphemptyX.SetPoint(j,float(pointlist[j][1]),float(pointlist[j][3]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphemptyX.SetMarkerStyle(20)
    GraphemptyX.Draw('ap')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphemptyX.GetHistogram().SetTitle("")
    GraphemptyX.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks in dXdz (mrad)")
    GraphemptyX.GetHistogram().GetYaxis().SetTitle("Gradient of asymmetry plot")
    GraphemptyX.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphemptyX.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphemptyX.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphemptyX.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphemptyX.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphemptyX,"Data","p")
    
    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphemptyX.pdf")
    c.Clear()

    GraphX = TGraphErrors(100)
    for j in range(0,100):
	    print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j][0]))
	    GraphX.SetPoint(j,float(pointlist[j][1]),float(pointlist[j][4]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphX.SetMarkerStyle(20)
    GraphX.Draw('ap')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphX.GetHistogram().SetTitle("")
    GraphX.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks in dXdz (mrad)")
    GraphX.GetHistogram().GetYaxis().SetTitle("Gradient of asymmetry plot")
    GraphX.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphX.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphX.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphX.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphX.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphX,"Data","p")
    
    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphX.pdf")

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
    rot_ang = next(x.prop("value") for x in doc.xpathEval("spec/sys") \
                        if x.prop("name").find("rot_ang") >= 0)
    print rot_ang
    thY_asym = rf.Get("thetaY_asymm")
    myfunc = thY_asym.GetFunction("pol1")
    grad = myfunc.GetParameter(1);
    vals.append(grad)
    vals.append(rot_ang)
    thY_asymempty = rf.Get("thetaY_emptyasymm")
    myfuncempty = thY_asymempty.GetFunction("pol1")
    gradempty = myfuncempty.GetParameter(1);
    vals.append(gradempty)
    thX_asymempty = rf.Get("thetaX_emptyasymm")
    myfuncXempty = thX_asymempty.GetFunction("pol1")
    gradXempty = myfuncXempty.GetParameter(1);
    vals.append(gradXempty)
    thX_asym = rf.Get("thetaX_asymm")
    myfuncX = thX_asym.GetFunction("pol1")
    gradX = myfuncX.GetParameter(1);
    vals.append(gradX)
    return vals

if __name__=="__main__":
    plotFiles(sys.argv[4:])
