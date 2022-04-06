import libxml2
import sys, os, subprocess, math, ROOT
from ROOT import TH1D, TCanvas, TFile, TGraphErrors, TMarker, TF1, TLegend, TText, TLine, TROOT, TGraph, TPad
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

    Graph = TGraphErrors(20)
    for j in range(0,20):
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
        print len(pointlist)
        print len(pointlist[j-1])
        Graph.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][0]))
        Graph.SetPointError(j,0,float(pointlist[j-1][36]))
    j = 0
    c = TCanvas()
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graph.SetMarkerStyle(20)
    #Graph.SetMarkerSize(4)
    Graph.Draw('AEP')
    line3 = TLine()
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graph.GetHistogram().SetTitle("")
    Graph.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    Graphempty = TGraphErrors(20)
    for j in range(0,20):
	    # print 'Graphempty'
	    # print pointlist[j-1][1],pointlist[j-1][25]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphempty.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][2]))
	    Graphempty.SetPointError(j,0,float(pointlist[j-1][37]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphempty.SetMarkerStyle(20)
    Graphempty.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphempty.GetHistogram().SetTitle("")
    Graphempty.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphemptyX = TGraphErrors(20)
    for j in range(0,20):
	    # print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphemptyX.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][3]))
	    GraphemptyX.SetPointError(j,0,float(pointlist[j-1][38]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphemptyX.SetMarkerStyle(20)
    GraphemptyX.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphemptyX.GetHistogram().SetTitle("")
    GraphemptyX.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphX = TGraphErrors(20)
    for j in range(0,20):
	    # print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphX.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][4]))
	    GraphX.SetPointError(j,0,float(pointlist[j-1][39]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphX.SetMarkerStyle(20)
    GraphX.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphX.GetHistogram().SetTitle("")
    GraphX.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphY = TGraphErrors(20)
    for j in range(0,20):
	    # print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphY.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][49]))
	    GraphY.SetPointError(j,0,float(pointlist[j-1][50]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphY.SetMarkerStyle(20)
    GraphY.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphY.GetHistogram().SetTitle("")
    GraphY.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
    GraphY.GetHistogram().GetYaxis().SetTitle("Gradient of asymmetry plot")
    GraphY.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphY.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphY.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphY.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphY.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphY,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphY.pdf")
    c.Clear()

    Graphmean = TGraphErrors(20)
    for j in range(0,20):
	    # print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphmean.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][5]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphmean.SetMarkerStyle(20)
    #Graph.SetMarkerSize(4)
    Graphmean.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphmean.GetHistogram().SetTitle("")
    Graphmean.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    Graphmeany = TGraphErrors(20)
    for j in range(0,20):
	    # print pointlist[j][1],pointlist[j][0]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphmeany.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][51]))
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphmeany.SetMarkerStyle(20)
    #Graph.SetMarkerSize(4)
    Graphmeany.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphmeany.GetHistogram().SetTitle("")
    Graphmeany.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
    Graphmeany.GetHistogram().GetYaxis().SetTitle("Mean of angle residual between US & DS tracks at absorber")
    Graphmeany.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphmeany.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphmeany.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphmeany.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphmeany.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    if Graphmeany.GetName() == "integralX":
        leg = TLegend(0.5,0.15,0.79,0.39)
    else:
        leg = TLegend(0.5,0.5,0.89,0.89)
    leg.SetLineColor(10)
    leg.AddEntry(Graphmeany,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphmeany.pdf")
    c.Clear()

    Graphgraphmean = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphgraphmean.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][10]))
	    Graphgraphmean.SetPointError(j,0,float(pointlist[j-1][44]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphgraphmean.SetMarkerStyle(20)
    Graphgraphmean.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphgraphmean.GetHistogram().SetTitle("")
    Graphgraphmean.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphrecoGoldmean = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphrecoGoldmean.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][11]))
	    GraphrecoGoldmean.SetPointError(j,0,float(pointlist[j-1][45]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphrecoGoldmean.SetMarkerStyle(20)
    GraphrecoGoldmean.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphrecoGoldmean.GetHistogram().SetTitle("")
    GraphrecoGoldmean.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    Graphdatamean = TGraph(18)
    k=80
    for j in range(0,20):
        # Graphdatamean.RemovePoint(0)
        print pointlist[j-1][1],pointlist[j-1][12]
        if(float(pointlist[j-1][1])!=0):
            Graphdatamean.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][12]))
            # Graphdatamean.SetPointError(j,0,float(pointlist[j-1][46]))
            k+=5
    for j in range(0,100):
            x = ROOT.Double(0)
            y = ROOT.Double(0)
            Graphdatamean.GetPoint(j,x,y)
            # print x, y
            # if (x==0.0):
            #     Graphdatamean.RemovePoint(j)
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphdatamean.SetMarkerStyle(20)
    Graphdatamean.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphdatamean.GetHistogram().SetTitle("")
    Graphdatamean.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
    Graphdatamean.GetHistogram().GetYaxis().SetTitle("Mean of scatter plot")
    Graphdatamean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphdatamean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphdatamean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphdatamean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphdatamean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphdatamean,"Data","p")
    fitmean = TF1("fitmean","pol1");
    print "about to do the thinabout to do the thingg"
    Graphdatamean.Fit("fitmean");
    mean_c = fitmean.GetParameter(0);
    mean_m = fitmean.GetParameter(1);
    fitmean.Draw("same")

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

    Graphdatameany = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
        # if (float(pointlist[j-1][51])<0.03):
            Graphdatameany.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][54]))
            Graphdatameany.SetPointError(j,0,float(pointlist[j-1][57]))
            k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphdatameany.SetMarkerStyle(20)
    Graphdatameany.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphdatameany.GetHistogram().SetTitle("")
    Graphdatameany.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
    Graphdatameany.GetHistogram().GetYaxis().SetTitle("Mean of scatter plot")
    Graphdatameany.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphdatameany.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphdatameany.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphdatameany.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphdatameany.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphdatameany,"Data","p")
    fitmeany = TF1("fitmeany","pol1");
    Graphdatameany.Fit("fitmeany");
    meany_c = fitmeany.GetParameter(0);
    meany_m = fitmeany.GetParameter(1);
    fitmeany.Draw("same")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphdatameany.pdf")
    c.Clear()

    Graphrefmean = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphrefmean.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][13]))
	    Graphrefmean.SetPointError(j,0,float(pointlist[j-1][46]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphrefmean.SetMarkerStyle(20)
    Graphrefmean.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphrefmean.GetHistogram().SetTitle("")
    Graphrefmean.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
    Graphrefmean.GetHistogram().GetYaxis().SetTitle("Mean of scatter plot")
    Graphrefmean.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphrefmean.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphrefmean.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphrefmean.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphrefmean.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphrefmean,"Data","p")
    fitrefmean = TF1("fitrefmean","pol1");
    Graphrefmean.Fit("fitrefmean");
    meanref_c = fitrefmean.GetParameter(0);
    meanref_m = fitrefmean.GetParameter(1);
    fitrefmean.Draw("same")

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

    Graphrefmeany = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphrefmeany.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][58]))
	    Graphrefmeany.SetPointError(j,0,float(pointlist[j-1][46]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphrefmeany.SetMarkerStyle(20)
    Graphrefmeany.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphrefmeany.GetHistogram().SetTitle("")
    Graphrefmeany.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
    Graphrefmeany.GetHistogram().GetYaxis().SetTitle("Mean of scatter plot")
    Graphrefmeany.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphrefmeany.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphrefmeany.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphrefmeany.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphrefmeany.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphrefmeany,"Data","p")
    fitrefmeany = TF1("fitrefmeany","pol1");
    Graphrefmeany.Fit("fitrefmeany");
    meanrefy_c = fitrefmeany.GetParameter(0);
    meanrefy_m = fitrefmeany.GetParameter(1);
    fitrefmeany.Draw("same")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphrefmeany.pdf")
    c.Clear()

    Graphgraphskew = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    #print pointlist[j-1][9],pointlist[j-1][8]

	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphgraphskew.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][14]))
	    Graphgraphskew.SetPointError(j,0,float(pointlist[j-1][40]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphgraphskew.SetMarkerStyle(20)
    Graphgraphskew.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphgraphskew.GetHistogram().SetTitle("")
    Graphgraphskew.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphrecoGoldskew = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphrecoGoldskew.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][15]))
	    GraphrecoGoldskew.SetPointError(j,0,float(pointlist[j-1][41]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphrecoGoldskew.SetMarkerStyle(20)
    GraphrecoGoldskew.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphrecoGoldskew.GetHistogram().SetTitle("")
    GraphrecoGoldskew.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphrecoGoldskewy = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphrecoGoldskewy.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][52]))
	    GraphrecoGoldskewy.SetPointError(j,0,float(pointlist[j-1][53]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphrecoGoldskewy.SetMarkerStyle(20)
    GraphrecoGoldskewy.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphrecoGoldskewy.GetHistogram().SetTitle("")
    GraphrecoGoldskewy.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
    GraphrecoGoldskewy.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    GraphrecoGoldskewy.GetHistogram().GetXaxis().SetLabelSize(0.05)
    GraphrecoGoldskewy.GetHistogram().GetYaxis().SetLabelSize(0.05)
    GraphrecoGoldskewy.GetHistogram().GetXaxis().SetTitleSize(0.05)
    GraphrecoGoldskewy.GetHistogram().GetYaxis().SetTitleSize(0.05)
    GraphrecoGoldskewy.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(GraphrecoGoldskewy,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("GraphrecoGoldskewy.pdf")
    c.Clear()

    Graphdataskew = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print 'dataskew'
	    # print j,0,float(pointlist[j-1][30])
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphdataskew.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][16]))
	    Graphdataskew.SetPointError(j,0,float(pointlist[j-1][42]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphdataskew.SetMarkerStyle(20)
    Graphdataskew.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphdataskew.GetHistogram().SetTitle("")
    Graphdataskew.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    Graphdataskewy = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print 'dataskew'
	    # print j,0,float(pointlist[j-1][30])
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphdataskewy.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][55]))
	    Graphdataskewy.SetPointError(j,0,float(pointlist[j-1][56]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphdataskewy.SetMarkerStyle(20)
    Graphdataskewy.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphdataskewy.GetHistogram().SetTitle("")
    Graphdataskewy.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
    Graphdataskewy.GetHistogram().GetYaxis().SetTitle("Skew of scatter plot")
    Graphdataskewy.GetHistogram().GetXaxis().SetLabelSize(0.05)
    Graphdataskewy.GetHistogram().GetYaxis().SetLabelSize(0.05)
    Graphdataskewy.GetHistogram().GetXaxis().SetTitleSize(0.05)
    Graphdataskewy.GetHistogram().GetYaxis().SetTitleSize(0.05)
    Graphdataskewy.GetHistogram().GetYaxis().SetTitleOffset(1.25)
    leg.SetLineColor(10)
    leg.AddEntry(Graphdataskewy,"Data","p")

    j+=1
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)
    t1.Draw()
    t2.Draw()
    c.Print("Graphdataskewy.pdf")
    c.Clear()

    Graphrefskew = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphrefskew.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][17]))
	    Graphrefskew.SetPointError(j,0,float(pointlist[j-1][43]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphrefskew.SetMarkerStyle(20)
    Graphrefskew.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphrefskew.GetHistogram().SetTitle("")
    Graphrefskew.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphalltofdXdzres_graphmean = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphalltofdXdzres_graphmean.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][18]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphalltofdXdzres_graphmean.SetMarkerStyle(20)
    GraphalltofdXdzres_graphmean.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphalltofdXdzres_graphmean.GetHistogram().SetTitle("")
    GraphalltofdXdzres_graphmean.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphalltofdXdzres_graphskew = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphalltofdXdzres_graphskew.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][19]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphalltofdXdzres_graphskew.SetMarkerStyle(20)
    GraphalltofdXdzres_graphskew.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphalltofdXdzres_graphskew.GetHistogram().SetTitle("")
    GraphalltofdXdzres_graphskew.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphalltofdYdzres_graphmean = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphalltofdYdzres_graphmean.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][20]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphalltofdYdzres_graphmean.SetMarkerStyle(20)
    GraphalltofdYdzres_graphmean.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphalltofdYdzres_graphmean.GetHistogram().SetTitle("")
    GraphalltofdYdzres_graphmean.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    GraphalltofdYdzres_graphskew = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    GraphalltofdYdzres_graphskew.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][21]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    GraphalltofdYdzres_graphskew.SetMarkerStyle(20)
    GraphalltofdYdzres_graphskew.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    GraphalltofdYdzres_graphskew.GetHistogram().SetTitle("")
    GraphalltofdYdzres_graphskew.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    Graphalltofpzres_graphmean = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphalltofpzres_graphmean.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][22]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphalltofpzres_graphmean.SetMarkerStyle(20)
    Graphalltofpzres_graphmean.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphalltofpzres_graphmean.GetHistogram().SetTitle("")
    Graphalltofpzres_graphmean.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    Graphalltofpzres_graphskew = TGraphErrors(20)
    k=80
    for j in range(0,20):
	    # print pointlist[j-1][9],pointlist[j-1][8]
	    #Graph.SetPoint(j,float(float(j-40)/2),float(pointlist[j-1][0]))
	    Graphalltofpzres_graphskew.SetPoint(j,float(pointlist[j-1][1]),float(pointlist[j-1][23]))
	    k+=5
    j = 0
    c.SetLeftMargin(0.125)
    c.SetBottomMargin(0.125)
    Graphalltofpzres_graphskew.SetMarkerStyle(20)
    Graphalltofpzres_graphskew.Draw('AEP')
    c.Update();
    line3.DrawLine(TROOT.gPad.GetUxmin(),0,TROOT.gPad.GetUxmax(),0);

    Graphalltofpzres_graphskew.GetHistogram().SetTitle("")
    Graphalltofpzres_graphskew.GetHistogram().GetXaxis().SetTitle("Rotation angle for upstream tracks")
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

    os.chdir('../')
    print os.getcwd()
    findCommand = "find . -path \"*scanrotang*\" -prune -o -name \"*.xml\""
    bashCommand = "find . -path \"*scanrotang*\" -prune -o -path \"*scanradial*\" -prune -o -name \"*.xml\" -exec sed -i s/\"rot_ang\"\ value=\".*\"/\"rot_ang\"\ value=\""
    bashCommand5 = r" -exec sed -i s/\"rot_ang\"\ value=\".*\"/\"rot_ang\"\ value=\""
    bashCommand2 = -mean_c/mean_m
    bashCommand3 = r"\"/g {} \;"""
    bashCommand4 = findCommand+bashCommand5+str(bashCommand2)+bashCommand3
    print findCommand
    print bashCommand4
    import subprocess
    processfind = subprocess.Popen(findCommand, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    process = subprocess.Popen(bashCommand4, stdout=subprocess.PIPE, shell=True)
    print process.stdout.read()
    # output, error = process.communicate()

    bashCommand = "find . -path \"*scanrotang*\" -prune -o -path \"*scanradial*\" -prune -o -name \"*.xml\" -exec sed -i s/\"rot_ang\"\ value=\".*\"/\"rot_ang\"\ value=\""
    bashCommand5 = r" -exec sed -i s/\"rot_angX\"\ value=\".*\"/\"rot_angX\"\ value=\""
    bashCommand2 = -meany_c/meany_m
    bashCommand3 = r"\"/g {} \;"""
    bashCommand4 = findCommand+bashCommand5+str(bashCommand2)+bashCommand3
    print findCommand
    print bashCommand4
    processfind = subprocess.Popen(findCommand, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    process = subprocess.Popen(bashCommand4, stdout=subprocess.PIPE, shell=True)
    print process.stdout.read()

    bashCommand = "find . -path \"*scanrotang*\" -prune -o -path \"*scanradial*\" -prune -o -name \"*.xml\" -exec sed -i s/\"rot_ang_empty\"\ value=\".*\"/\"rot_ang_empty\"\ value=\""
    bashCommand5 = r" -exec sed -i s/\"rot_ang_empty\"\ value=\".*\"/\"rot_ang_empty\"\ value=\""
    bashCommand2 = -meanref_c/meanref_m
    bashCommand3 = r"\"/g {} \;"""
    bashCommand4 = findCommand+bashCommand5+str(bashCommand2)+bashCommand3
    print findCommand
    print bashCommand4
    processfind = subprocess.Popen(findCommand, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    process = subprocess.Popen(bashCommand4, stdout=subprocess.PIPE, shell=True)
    print process.stdout.read()
    # output, error = process.communicate()

    bashCommand = "find . -path \"*scanrotang*\" -prune -o -path \"*scanradial*\" -prune -o -name \"*.xml\" -exec sed -i s/\"rot_ang_emptyX\"\ value=\".*\"/\"rot_ang_emptyX\"\ value=\""
    bashCommand5 = r" -exec sed -i s/\"rot_ang_emptyX\"\ value=\".*\"/\"rot_ang_emptyX\"\ value=\""
    bashCommand2 = -meanrefy_c/meanrefy_m
    bashCommand3 = r"\"/g {} \;"""
    bashCommand4 = findCommand+bashCommand5+str(bashCommand2)+bashCommand3
    print findCommand
    print bashCommand4
    processfind = subprocess.Popen(findCommand, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    process = subprocess.Popen(bashCommand4, stdout=subprocess.PIPE, shell=True)
    print process.stdout.read()

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
    print thY_asym.FindObject("fM2D")
    myfunc = thY_asym.GetFunction("fM2D")
    myfunc.Print()
    grad = myfunc.GetParameters();
    graderr = myfunc.GetParError(1);
    print 'grad',grad[1]
    print 'graderr',graderr
    vals.append(grad[1])
    vals.append(rot_ang)
    thY_asymempty = rf.Get("thetaY_emptyasymm")
    myfuncempty = thY_asymempty.GetFunction("fM2D")
    myfuncempty.Print()
    gradempty = myfuncempty.GetParameters();
    vals.append(gradempty[1])
    thX_asymempty = rf.Get("thetaX_emptyasymm")
    myfuncXempty = thX_asymempty.GetFunction("fM2D")
    gradXempty = myfuncXempty.GetParameters();
    vals.append(gradXempty[1])
    thX_asym = rf.Get("thetaX_asymm")
    myfuncX = thX_asym.GetFunction("fM2D")
    gradX = myfuncX.GetParameters();
    graderrX = myfuncX.GetParError(1);
    print 'here'
    print 'gradX',gradX[1]
    print 'graderrX',graderrX
    vals.append(gradX[1])
    meanabsLiH = rf.Get("dataAcc_projRefDiffddz")
    meanabsLiHval = meanabsLiH.GetMean(1)
    vals.append(meanabsLiHval)
    MCthX_asym = rf.Get("MCthetaX_asymm")
    myMCfuncX = MCthX_asym.GetFunction("pol1")
    MCgradX = 1
    vals.append(MCgradX)
    MCthY_asym = rf.Get("MCthetaY_asymm")
    myMCfuncY = MCthY_asym.GetFunction("pol1")
    MCgradY = 1
    vals.append(MCgradY)
    thetaX_graph = rf.Get("thetaX_graphmc")
    thetaX_recoGold = rf.Get("thetaX_recoGold")
    thetaX_data = rf.Get("thetaX_data")
    thetaX_ref = rf.Get("thetaX_ref")
    chi2=0
    vals.append(chi2)
    vals.append(fid_rad)
    #10 below counting from 0
    thetaX_graphmean=0
    vals.append(thetaX_graphmean)
    thetaX_recoGoldmean=thetaX_recoGold.GetMean()
    vals.append(thetaX_recoGoldmean)
    thetaX_datamean=thetaX_data.GetMean()
    vals.append(thetaX_datamean)
    thetaX_refmean=thetaX_ref.GetMean()
    vals.append(thetaX_refmean)
    thetaX_graphskew=0
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
    #20 below counting from 0
    alltofdYdzres_graphmean=alltofdYdzres.GetMean()
    vals.append(alltofdYdzres_graphmean)
    alltofdYdzres_graphskew=alltofdYdzres.GetSkewness()
    vals.append(alltofdYdzres_graphskew)
    alltofpzres_graphmean=alltofpzres.GetMean()
    vals.append(alltofpzres_graphmean)
    alltofpzres_graphskew=alltofpzres.GetSkewness()
    vals.append(alltofpzres_graphskew)
    preraddXdzres = rf.Get("dataUS_dif_dXdzres")
    preraddYdzres = rf.Get("dataUS_dif_dYdzres")
    preradpzres = rf.Get("dataUS_dif_pzres")
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
    #30 below counting from 0
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
    #40 below counting from 0
    #thetaX_graphskewerr=math.sqrt(6/thetaX_graph.GetEntries())
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

    kolmogorov=0
    vals.append(chi2)

    thY_asym = rf.Get("thetaY_asymm")
    myfuncY = thY_asym.GetFunction("fM2D")
    gradY = myfuncY.GetParameters();
    graderrY = myfuncY.GetParError(1);
    vals.append(gradY[1])
    #50 below counting from 0
    gradYerr = myfuncY.GetParError(1);
    vals.append(gradYerr)
    thetaY_recoGold = rf.Get("thetaY_recoGold")
    thetaY_data = rf.Get("thetaY_data")
    thetaY_ref = rf.Get("thetaY_ref")
    thetaY_recoGoldmean=thetaY_recoGold.GetMean()
    vals.append(thetaY_recoGoldmean)
    thetaY_recoGoldskew=thetaY_recoGold.GetSkewness()
    vals.append(thetaY_recoGoldskew)
    thetaY_recoGoldskewerr=math.sqrt(6/thetaY_recoGold.GetEntries())
    vals.append(thetaY_recoGoldskewerr)
    thetaY_datamean=thetaY_data.GetMean()
    vals.append(thetaY_datamean)
    thetaY_datamean=thetaY_data.GetMean()
    thetaY_dataskew=thetaY_data.GetSkewness()
    vals.append(thetaY_dataskew)
    thetaY_dataskewerr=math.sqrt(6/thetaY_data.GetEntries())
    vals.append(thetaY_dataskewerr)
    thetaY_datameanerr=thetaY_data.GetMeanError()
    vals.append(thetaY_datameanerr)


    thetaY_refmean=thetaY_ref.GetMean()
    vals.append(thetaY_refmean)
    # chi2 comparison to MC
    # thetaX_data.chi2(thetaX_grahmc)
    # thetaY_data.chi2(thetaY_grahmc)
    # thetaX_ref.chi2(thetaX_grahmc)
    # thetaY_ref.chi2(thetaX_grahmc)

    # vals.append(thetaX_datachi2)
    # vals.append(thetaY_datachi2)
    # vals.append(thetaY_refchi2)
    # vals.append(thetaY_refchi2)
    return vals

if __name__=="__main__":
    plotFiles(sys.argv[1:])
