import libxml2
import sys, os, subprocess, math, ROOT
from ROOT import TH1D, TCanvas, TFile, TGraphErrors, TMarker, TF1, TLegend, TText, TLine, TROOT
from pprint import pprint

TROOT.gROOT.SetBatch(1)
def plotFiles(syslist):

    hists  = []
    Graphs = []
    pointlist = []

    for sysfile in syslist:
        sysval = ExtractPars(sysfile)
        pointlist.append(sysval)

    pprint(pointlist)

    Graph = TGraphErrors(12)
    for j in range(0,11):

	    print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graph.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][0]))
	    Graph.SetPointError(j,0,float(pointlist[j-1][36]))
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
    Graph.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
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

    Graphempty = TGraphErrors(12)
    for j in range(0,11):
	    print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphempty.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][2]))
	    Graphempty.SetPointError(j,0,float(pointlist[j-1][37]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphempty.SetMarkerStyle(20)
    Graphempty.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphempty.GetHistogram().SetTitle("")
    Graphempty.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
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

    GraphemptyX = TGraphErrors(12)
    for j in range(0,11):
	    print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphemptyX.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][3]))
	    GraphemptyX.SetPointError(j,0,float(pointlist[j-1][38]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphemptyX.SetMarkerStyle(20)
    GraphemptyX.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0111,0);

    GraphemptyX.GetHistogram().SetTitle("")
    GraphemptyX.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
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

    GraphX = TGraphErrors(12)
    for j in range(0,11):
	    print 'fhns'
	    print pointlist[j-1][9],pointlist[j-1][4]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphX.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][4]))
	    GraphX.SetPointError(j,0,float(pointlist[j-1][39]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphX.SetMarkerStyle(20)
    GraphX.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0111,0);

    GraphX.GetHistogram().SetTitle("")
    GraphX.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
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
    c.Clear()

    Graphmean = TGraphErrors(12)
    for j in range(0,11):
	    print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphmean.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][5]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphmean.SetMarkerStyle(20)
    #Graph.SetMarkerSize(4)
    Graphmean.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphmean.GetHistogram().SetTitle("")
    Graphmean.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphmean.GetHistogram().GetYaxis().SetTitle("Mean of angle residual between US & DS tracks at absorber")
    Graphmean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphmean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphmean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphmean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphmean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    if Graphmean.GetName() == "integralX":
        leg = TLegend(0.5,0.15,0.79,0.39)
    else:
        leg = TLegend(0.5,0.5,0.89,0.89)
    leg.SetLineColor(10)
    leg.AddEntry(Graphmean,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphmean.pdf")
    c.Clear()

    GraphMCX = TGraphErrors(12)
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][6]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphMCX.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][6]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphMCX.SetMarkerStyle(20)
    GraphMCX.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphMCX.GetHistogram().SetTitle("")
    GraphMCX.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    GraphMCX.GetHistogram().GetYaxis().SetTitle("Gradient of asymmetry plot")
    GraphMCX.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphMCX.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphMCX.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphMCX.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphMCX.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphMCX,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphMCX.pdf")
    c.Clear()

    GraphMCY = TGraphErrors(12)
    for j in range(0,11):
	    print float(pointlist[j-1][9]),pointlist[j-1][7]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphMCY.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][7]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphMCY.SetMarkerStyle(20)
    GraphMCY.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphMCY.GetHistogram().SetTitle("")
    GraphMCY.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    GraphMCY.GetHistogram().GetYaxis().SetTitle("Gradient of asymmetry plot")
    GraphMCY.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphMCY.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphMCY.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphMCY.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphMCY.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphMCY,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphMCY.pdf")
    c.Clear()

    Graphchi2 = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphchi2.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][8]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphchi2.SetMarkerStyle(20)
    Graphchi2.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphchi2.GetHistogram().SetTitle("")
    Graphchi2.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphchi2.GetHistogram().GetYaxis().SetTitle("#chi^2/NDF recon & Truth for scatter plot")
    Graphchi2.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphchi2.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphchi2.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphchi2.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphchi2.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphchi2,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphchi2.pdf")
    c.Clear()

    Graphkolmogorov2 = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print float(pointlist[j-1][9]),float(pointlist[j-1][48])
	    #print pointlist[j-1][9],pointlist[j-1][36]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphkolmogorov2.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][48]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphkolmogorov2.SetMarkerStyle(20)
    Graphkolmogorov2.Draw("AEP")
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphkolmogorov2.GetHistogram().SetTitle("")
    Graphkolmogorov2.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphkolmogorov2.GetHistogram().GetYaxis().SetTitle("kolmogorov recon & Truth for scatter plot")
    Graphkolmogorov2.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphkolmogorov2.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphkolmogorov2.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphkolmogorov2.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphkolmogorov2.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphkolmogorov2,"Data","p")
    #Graphkolmogorov2.GetYaxis().SetRangeUser(0.04,0.16)
    #Graphkolmogorov2.GetXaxis().SetRangeUser(80,150)

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphkolmogorov2.pdf")
    c.Clear()

    Graphgraphmean = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphgraphmean.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][10]))
	    Graphgraphmean.SetPointError(j,0,float(pointlist[j-1][32]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphgraphmean.SetMarkerStyle(20)
    Graphgraphmean.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphgraphmean.GetHistogram().SetTitle("")
    Graphgraphmean.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphgraphmean.GetHistogram().GetYaxis().SetTitle("Mean of scatter plot")
    Graphgraphmean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphgraphmean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphgraphmean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphgraphmean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphgraphmean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphgraphmean,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphgraphmean.pdf")
    c.Clear()

    GraphrecoGoldmean = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphrecoGoldmean.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][11]))
	    GraphrecoGoldmean.SetPointError(j,0,float(pointlist[j-1][45]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphrecoGoldmean.SetMarkerStyle(20)
    GraphrecoGoldmean.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphrecoGoldmean.GetHistogram().SetTitle("")
    GraphrecoGoldmean.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    GraphrecoGoldmean.GetHistogram().GetYaxis().SetTitle("Mean of scatter plot")
    GraphrecoGoldmean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphrecoGoldmean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphrecoGoldmean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphrecoGoldmean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphrecoGoldmean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphrecoGoldmean,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphrecoGoldmean.pdf")
    c.Clear()

    Graphdatamean = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphdatamean.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][12]))
	    Graphdatamean.SetPointError(j,0,float(pointlist[j-1][46]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphdatamean.SetMarkerStyle(20)
    Graphdatamean.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphdatamean.GetHistogram().SetTitle("")
    Graphdatamean.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphdatamean.GetHistogram().GetYaxis().SetTitle("Mean of scatter plot")
    Graphdatamean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphdatamean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphdatamean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphdatamean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphdatamean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphdatamean,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphdatamean.pdf")
    c.Clear()

    Graphrefmean = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphrefmean.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][13]))
	    Graphrefmean.SetPointError(j,0,float(pointlist[j-1][46]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphrefmean.SetMarkerStyle(20)
    Graphrefmean.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphrefmean.GetHistogram().SetTitle("")
    Graphrefmean.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphrefmean.GetHistogram().GetYaxis().SetTitle("Mean of scatter plot")
    Graphrefmean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphrefmean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphrefmean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphrefmean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphrefmean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphrefmean,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphrefmean.pdf")
    c.Clear()

    Graphgraphskew = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphgraphskew.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][14]))
	    Graphgraphskew.SetPointError(j,0,float(pointlist[j-1][40]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphgraphskew.SetMarkerStyle(20)
    Graphgraphskew.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphgraphskew.GetHistogram().SetTitle("")
    Graphgraphskew.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphgraphskew.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    Graphgraphskew.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphgraphskew.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphgraphskew.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphgraphskew.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphgraphskew.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphgraphskew,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphgraphskew.pdf")
    c.Clear()

    GraphrecoGoldskew = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphrecoGoldskew.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][15]))
	    GraphrecoGoldskew.SetPointError(j,0,float(pointlist[j-1][41]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphrecoGoldskew.SetMarkerStyle(20)
    GraphrecoGoldskew.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphrecoGoldskew.GetHistogram().SetTitle("")
    GraphrecoGoldskew.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    GraphrecoGoldskew.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    GraphrecoGoldskew.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphrecoGoldskew.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphrecoGoldskew.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphrecoGoldskew.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphrecoGoldskew.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphrecoGoldskew,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphrecoGoldskew.pdf")
    c.Clear()

    Graphdataskew = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphdataskew.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][16]))
	    Graphdataskew.SetPointError(j,0,float(pointlist[j-1][42]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphdataskew.SetMarkerStyle(20)
    Graphdataskew.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphdataskew.GetHistogram().SetTitle("")
    Graphdataskew.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphdataskew.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    Graphdataskew.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphdataskew.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphdataskew.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphdataskew.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphdataskew.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphdataskew,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphdataskew.pdf")
    c.Clear()

    Graphrefskew = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphrefskew.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][17]))
	    Graphrefskew.SetPointError(j,0,float(pointlist[j-1][43]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphrefskew.SetMarkerStyle(20)
    Graphrefskew.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphrefskew.GetHistogram().SetTitle("")
    Graphrefskew.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphrefskew.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    Graphrefskew.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphrefskew.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphrefskew.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphrefskew.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphrefskew.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphrefskew,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphrefskew.pdf")
    c.Clear()

    GraphalltofdXdzres_graphmean = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphalltofdXdzres_graphmean.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][18]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphalltofdXdzres_graphmean.SetMarkerStyle(20)
    GraphalltofdXdzres_graphmean.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphalltofdXdzres_graphmean.GetHistogram().SetTitle("")
    GraphalltofdXdzres_graphmean.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    GraphalltofdXdzres_graphmean.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    GraphalltofdXdzres_graphmean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphalltofdXdzres_graphmean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphalltofdXdzres_graphmean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphalltofdXdzres_graphmean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphalltofdXdzres_graphmean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphalltofdXdzres_graphmean,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphalltofdXdzres_graphmean.pdf")
    c.Clear()

    GraphalltofdXdzres_graphskew = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphalltofdXdzres_graphskew.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][19]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphalltofdXdzres_graphskew.SetMarkerStyle(20)
    GraphalltofdXdzres_graphskew.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphalltofdXdzres_graphskew.GetHistogram().SetTitle("")
    GraphalltofdXdzres_graphskew.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    GraphalltofdXdzres_graphskew.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    GraphalltofdXdzres_graphskew.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphalltofdXdzres_graphskew.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphalltofdXdzres_graphskew.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphalltofdXdzres_graphskew.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphalltofdXdzres_graphskew.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphalltofdXdzres_graphskew,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphalltofdXdzres_graphskew.pdf")
    c.Clear()

    GraphalltofdYdzres_graphmean = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphalltofdYdzres_graphmean.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][20]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphalltofdYdzres_graphmean.SetMarkerStyle(20)
    GraphalltofdYdzres_graphmean.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphalltofdYdzres_graphmean.GetHistogram().SetTitle("")
    GraphalltofdYdzres_graphmean.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    GraphalltofdYdzres_graphmean.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    GraphalltofdYdzres_graphmean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphalltofdYdzres_graphmean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphalltofdYdzres_graphmean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphalltofdYdzres_graphmean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphalltofdYdzres_graphmean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphalltofdYdzres_graphmean,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphalltofdYdzres_graphmean.pdf")
    c.Clear()

    GraphalltofdYdzres_graphskew = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphalltofdYdzres_graphskew.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][21]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphalltofdYdzres_graphskew.SetMarkerStyle(20)
    GraphalltofdYdzres_graphskew.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    GraphalltofdYdzres_graphskew.GetHistogram().SetTitle("")
    GraphalltofdYdzres_graphskew.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    GraphalltofdYdzres_graphskew.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    GraphalltofdYdzres_graphskew.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphalltofdYdzres_graphskew.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphalltofdYdzres_graphskew.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphalltofdYdzres_graphskew.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphalltofdYdzres_graphskew.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphalltofdYdzres_graphskew,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphalltofdYdzres_graphskew.pdf")
    c.Clear()

    Graphalltofpzres_graphmean = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphalltofpzres_graphmean.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][22]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphalltofpzres_graphmean.SetMarkerStyle(20)
    Graphalltofpzres_graphmean.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphalltofpzres_graphmean.GetHistogram().SetTitle("")
    Graphalltofpzres_graphmean.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphalltofpzres_graphmean.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    Graphalltofpzres_graphmean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphalltofpzres_graphmean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphalltofpzres_graphmean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphalltofpzres_graphmean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphalltofpzres_graphmean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphalltofpzres_graphmean,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphalltofpzres_graphmean.pdf")
    c.Clear()

    Graphalltofpzres_graphskew = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphalltofpzres_graphskew.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][23]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphalltofpzres_graphskew.SetMarkerStyle(20)
    Graphalltofpzres_graphskew.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphalltofpzres_graphskew.GetHistogram().SetTitle("")
    Graphalltofpzres_graphskew.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphalltofpzres_graphskew.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    Graphalltofpzres_graphskew.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphalltofpzres_graphskew.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphalltofpzres_graphskew.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphalltofpzres_graphskew.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphalltofpzres_graphskew.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphalltofpzres_graphskew,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphalltofpzres_graphskew.pdf")
    c.Clear()

    Graphsumresiduallr = TGraphErrors(12)
    k=80
    for j in range(0,11):
	    print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphsumresiduallr.SetPoint(j,float(pointlist[j-1][9]),float(pointlist[j-1][49]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphsumresiduallr.SetMarkerStyle(20)
    Graphsumresiduallr.Draw('AEP')
    line3.DrawLine(-0.012,0,0.0115,0);

    Graphsumresiduallr.GetHistogram().SetTitle("")
    Graphsumresiduallr.GetHistogram().GetXaxis().SetTitle("Radial cut (mm)")
    Graphsumresiduallr.GetHistogram().GetYaxis().SetTitle("Sum residual l/r")
    Graphsumresiduallr.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphsumresiduallr.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphsumresiduallr.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphsumresiduallr.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphsumresiduallr.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphsumresiduallr,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphsumresiduallr.pdf")
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
    rot_ang = next(x.prop("value") for x in doc.xpathEval("spec/sys") \
                        if x.prop("name").find("rot_ang") >= 0)
    print rot_ang
    fid_rad = next(x.prop("value") for x in doc.xpathEval("spec/cuts") \
                        if x.prop("name").find("fid_rad") >= 0)
    thY_asym = rf.Get("thetaY_asymm")
    myfunc = thY_asym.GetFunction("fM2D")
    grad = myfunc.GetParameters();
    vals.append(grad[1])
    vals.append(rot_ang)
    thY_asymempty = rf.Get("thetaY_emptyasymm")
    myfuncempty = thY_asymempty.GetFunction("fM2D")
    gradempty = myfuncempty.GetParameters();
    vals.append(gradempty[1])
    thX_asymempty = rf.Get("thetaX_emptyasymm")
    myfuncXempty = thX_asymempty.GetFunction("fM2D")
    gradXempty = myfuncXempty.GetParameters();
    vals.append(gradXempty[1])
    thX_asym = rf.Get("thetaX_asymm")
    myfuncX = thX_asym.GetFunction("fM2D")
    gradX = myfuncX.GetParameters();
    vals.append(gradX[1])
    meanabsLiH = rf.Get("dataAcc_projRefDiffddz")
    meanabsLiHval = meanabsLiH.GetMean(1)
    vals.append(meanabsLiHval)
    MCthX_asym = rf.Get("MCthetaX_asymm")
    myMCfuncX = MCthX_asym.GetFunction("fM2D")
    MCgradX = myMCfuncX.GetParameters();
    vals.append(MCgradX[1])
    MCthY_asym = rf.Get("MCthetaY_asymm")
    myMCfuncY = MCthY_asym.GetFunction("fM2D")
    MCgradY = myMCfuncY.GetParameters();
    vals.append(MCgradY[1])
    thetaX_graph = rf.Get("thetaX_graph")
    thetaX_recoGold = rf.Get("thetaX_recoGold_effi_only_gold")
    thetaX_data = rf.Get("thetaX_data")
    thetaX_ref = rf.Get("thetaX_ref")
    thetaX_graph.GetXaxis().SetRangeUser(-0.04,0.04)
    thetaX_recoGold.GetXaxis().SetRangeUser(-0.04,0.04)
    chi2=thetaX_graph.Chi2Test(thetaX_recoGold,"CHI2/NDF")
    # chi2=thetaX_data.Chi2Test(thetaX_graph,"CHI2/NDF")
    vals.append(chi2)
    vals.append(fid_rad)
    #10
    thetaX_graphmean=thetaX_graph.GetMean()
    vals.append(thetaX_graphmean)
    thetaX_recoGoldmean=thetaX_recoGold.GetMean()
    vals.append(thetaX_recoGoldmean)
    thetaX_datamean=thetaX_data.GetMean()
    vals.append(thetaX_datamean)
    thetaX_refmean=thetaX_ref.GetMean()
    vals.append(thetaX_refmean)
    thetaX_graphskew=thetaX_graph.GetSkewness()
    vals.append(thetaX_graphskew)
    thetaX_recoGoldskew=thetaX_recoGold.GetSkewness()
    vals.append(thetaX_recoGoldskew)
    thetaX_dataskew=thetaX_data.GetSkewness()
    vals.append(thetaX_dataskew)
    thetaX_refskew=thetaX_ref.GetSkewness()
    vals.append(thetaX_refskew)
    alltofdXdzres = rf.Get("dataUS_alltof_dXdzres")
    alltofdYdzres = rf.Get("dataUS_alltof_dYdzres")
    alltofpzres = rf.Get("dataUS_alltof_pzres")
    alltofdXdzres_graphmean=alltofdXdzres.GetMean()
    vals.append(alltofdXdzres_graphmean)
    alltofdXdzres_graphskew=alltofdXdzres.GetSkewness()
    vals.append(alltofdXdzres_graphskew)
    #20
    alltofdYdzres_graphmean=alltofdYdzres.GetMean()
    vals.append(alltofdYdzres_graphmean)
    alltofdYdzres_graphskew=alltofdYdzres.GetSkewness()
    vals.append(alltofdYdzres_graphskew)
    alltofpzres_graphmean=alltofpzres.GetMean()
    vals.append(alltofpzres_graphmean)
    alltofpzres_graphskew=alltofpzres.GetSkewness()
    vals.append(alltofpzres_graphskew)
    preraddXdzres = rf.Get("dataUS_fid_dXdzres")
    preraddYdzres = rf.Get("dataUS_fid_dYdzres")
    preradpzres = rf.Get("dataUS_fid_pzres")
    preraddXdzres_graphmean=preraddXdzres.GetMean()
    vals.append(preraddXdzres_graphmean)
    preraddXdzres_graphskew=preraddXdzres.GetSkewness()
    vals.append(preraddXdzres_graphskew)
    preraddYdzres_graphmean=preraddYdzres.GetMean()
    vals.append(preraddYdzres_graphmean)
    preraddYdzres_graphskew=preraddYdzres.GetSkewness()
    vals.append(preraddYdzres_graphskew)
    preradpzres_graphmean=preradpzres.GetMean()
    vals.append(preradpzres_graphmean)
    preradpzres_graphskew=preradpzres.GetSkewness()
    vals.append(preradpzres_graphskew)
    #30
    dXdzres = rf.Get("dataUS_dXdzres")
    dYdzres = rf.Get("dataUS_dYdzres")
    pzres = rf.Get("dataUS_pzres")
    dXdzres_graphmean=dXdzres.GetMean()
    vals.append(dXdzres_graphmean)
    dXdzres_graphskew=dXdzres.GetSkewness()
    vals.append(dXdzres_graphskew)
    dYdzres_graphmean=dYdzres.GetMean()
    vals.append(dYdzres_graphmean)
    dYdzres_graphskew=dYdzres.GetSkewness()
    vals.append(dYdzres_graphskew)
    pzres_graphmean=pzres.GetMean()
    vals.append(pzres_graphmean)
    pzres_graphskew=pzres.GetSkewness()
    vals.append(pzres_graphskew)

    graderr = myfunc.GetParError(1);
    vals.append(graderr)
    grademptyerr = myfuncempty.GetParError(1);
    vals.append(grademptyerr)
    gradXemptyerr = myfuncXempty.GetParError(1);
    vals.append(gradXemptyerr)
    gradXerr = myfuncX.GetParError(1);
    vals.append(gradXerr)
    #40
    #thetaX_graphskewerr=thetaX_graph.GetSkewness(11)
    thetaX_graphskewerr=0
    vals.append(thetaX_graphskewerr)
    #thetaX_recoGoldskewerr=thetaX_recoGold.GetSkewness(11)
    thetaX_recoGoldskewerr=math.sqrt(6/thetaX_recoGold.GetEntries())
    vals.append(thetaX_recoGoldskewerr)
    #thetaX_dataskewerr=thetaX_data.GetSkewness(11)
    thetaX_dataskewerr=math.sqrt(6/thetaX_data.GetEntries())
    vals.append(thetaX_dataskewerr)
    #thetaX_refskewerr=thetaX_ref.GetSkewness(11)
    thetaX_refskewerr=math.sqrt(6/thetaX_ref.GetEntries())
    vals.append(thetaX_refskewerr)
    thetaX_graphmeanerr=0
    vals.append(thetaX_graphmeanerr)
    thetaX_recoGoldmeanerr=thetaX_recoGold.GetMeanError()
    vals.append(thetaX_recoGoldmeanerr)
    thetaX_datameanerr=thetaX_data.GetMeanError()
    vals.append(thetaX_datameanerr)
    thetaX_refmeanerr=thetaX_ref.GetMeanError()
    vals.append(thetaX_refmeanerr)

    kolmogorov=thetaX_graph.KolmogorovTest(thetaX_recoGold, "M")
    vals.append(kolmogorov)

    htruthmirrorl=rf.Get("htruthmirrorl")
    htruthmirrorr=rf.Get("htruthmirrorr")
    htruthl=rf.Get("thetaX_graph")
    nhtruthmirrorl = htruthl.GetEntries()
    htruthmirrorl.Scale(1/nhtruthmirrorl)
    nhtruthmirrorr = htruthl.GetEntries()
    htruthmirrorr.Scale(1/nhtruthmirrorr)
    htruthmirrorl.Draw()
    htruthmirrorr.Draw()
    residual = 0
    for i in range(1,47/2):
           yparl = htruthmirrorl.GetBinContent(i);
           yparr = htruthmirrorr.GetBinContent(i);
           residual += math.sqrt(pow(yparl-yparr,2));

    vals.append(residual)
    #50
    return vals

if __name__=="__main__":
    plotFiles(sys.argv[4:])
