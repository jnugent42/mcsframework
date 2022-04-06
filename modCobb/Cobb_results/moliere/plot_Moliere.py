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
        
    thetaX_Moliere = TH1D("thetaX_Moliere",";#theta_{Scatt};Events (arbitrary normalization)",len(data)-1,xbins)
    thetaY_Moliere = TH1D("thetaY_Moliere",";#theta_{Scatt};Events (arbitrary normalization)",len(data)-1,xbins)
    
    for d in data:
        thetaX_Moliere.SetBinContent(d[0] + 1, d[2])
        thetaY_Moliere.SetBinContent(d[0] + 1, d[2])
    
    
    c.SetLogy()
    
    g.SetMarkerStyle(20)
    g.SetTitle(";#theta_{Scatt};Events (arbitrary normalization)")
    g.Draw("ap")
    thetaX_Moliere.Draw("lsame")
    thetaY_Moliere.Draw("lsame")
    c.Print("Moliere_200_1.pdf")
    fout.cd()
    thetaX_Moliere.Write()
    thetaY_Moliere.Write()
    fout.Close()
    
    hthetaX = TH1D("hthetaX",";#theta_{X}; Events per 2 mrad",len(x2bins)-1,x2bins)
    hthetaY = TH1D("hthetaY",";#theta_{Y}; Events per 2 mrad",len(x2bins)-1,x2bins)
    rand = TRandom3()
    for i in range(100000):
        theta = h.GetRandom()
        # phi   = (rand.Rndm()) * math.atan(1) * 8.0
        # print phi
        X = (2*rand.Rndm() - 1)  # *math.cos(phi)
        Y = (2*rand.Rndm() - 1) #  *math.sin(phi)
        norm =  math.sqrt(X*X + Y*Y)
        X /= norm
        Y /= norm
        Z = math.cos(theta)
        thetaX = math.atan2(Y* math.sin(theta), Z)
        thetaY = math.atan2(X* math.sin(theta), Z)
        hthetaX.Fill(thetaX)
        hthetaY.Fill(thetaY)
    
    hthetaX.Draw()
    c.Print("Moliere_200_thetaX.pdf")
    
    hthetaY.Draw()
    c.Print("Moliere_200_thetaY.pdf")

if __name__ == '__main__':
    builddistro(sys.argv[1],sys.argv[2])
