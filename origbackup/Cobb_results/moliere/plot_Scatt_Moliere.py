from ROOT import TGraph, TCanvas, TH1D, TRandom3, TFile
from array import array
import math, sys

def builddistro(fname, outname):

    fout = TFile(outname, "UPDATE")
    fin = open(fname,'r')
    i = 0
    data = []
    last = 0
    for l in fin:
        ln = l.split()
        data.append([i, float(ln[0]), float(ln[1]), float(ln[0]) - last])
        last = float(ln[0])
        i+=1
    c = TCanvas()
    g = TGraph(len(data))
    xbins = array('d')
    x2bins = array('d')
    data[0][3] = data[1][3]
    for d in data:
        g.SetPoint(d[0], d[1], d[2])
        xbins.append(d[1] - d[3]/2.)
        x2bins.append(d[1] + d[3]/2.)
        x2bins.insert(0,-d[1] - d[3]/2.)
        print d[1] - d[3]/2., d[1], d[3]
        
    h = TH1D("thetaScatt_Moliere",";#theta_{Scatt};Events (arbitrary normalization)",len(data)-1,xbins)
    
    for d in data:
        h.SetBinContent(d[0] + 1, d[2])
    
    
    c.SetLogy()
    
    g.SetMarkerStyle(20)
    g.SetTitle(";#theta_{Scatt};Events (arbitrary normalization)")
    g.Draw("ap")
    h.Draw("lsame")
    fout.cd()
    h.Write()
    fout.Close()

if __name__ == '__main__':
    builddistro(sys.argv[1],sys.argv[2])
