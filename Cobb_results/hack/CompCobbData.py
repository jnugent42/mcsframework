from ROOT import TH1D, TFile, TCanvas, TLegend, TGraph, gDirectory
from array import array
import sys, math, os

def dataModelComps(tname, outfile):

    fout = TFile(outfile, 'RECREATE') # (tname[:-4] + '.root','RECREATE')
    histvarnames = ['thetaX','thetaY','thetaScatt']
    arrayindex   = [2, 3, 0]
    arraymin     = [-0.2, -0.2, 0.0]
    histstatenames = ['recoGEANT','GEANT']
    histstatedesc = ['Deconvolved data',
                     'GEANT4', 'Cobb-Carlisle']
    histcolors = [1, 4, 2]
    histlines = [1, 4, 8]
    histopts  = ['lp','l','l']
    lines = [line.rstrip('\n').split() for line in open(tname)]
    print lines
    q = 0
    for histvar in histvarnames:
        print [histvar+'_'+x for x in histstatenames]
        hists = []
        hists.append(TH1D(histvar+'_Cobb',histvar,800,arraymin[q],0.2))
        offset = 1
        for j in range(len(lines)):
            if histvar!="thetaScatt":
                print j, hists[-1].GetXaxis().GetBinCenter(j), float(lines[j][arrayindex[q]])
            hists[-1].SetBinContent(j+offset, float(lines[j][arrayindex[q]]))
            hists[-1].SetBinError(j+offset, math.sqrt(float(lines[j][arrayindex[q]])))
        q += 1
        hists[0].SetMarkerStyle(20)
        hists[0].SetMarkerColor(histcolors[0])
        hists[0].GetXaxis().SetLabelSize(0.05)
        hists[0].GetXaxis().SetTitleSize(0.05)
        hists[0].GetYaxis().SetLabelSize(0.05)
        hists[0].GetYaxis().SetTitleSize(0.05)
        if histvar == 'thetaScatt':
            leg = TLegend(0.65,0.5,0.89,0.89)
            hists[0].GetYaxis().SetTitle("Probability per mrad")
        else:
            leg = TLegend(0.35,0.2,0.65,0.5)
            hists[0].GetYaxis().SetTitle("Probability per 2 mrad")
        leg2 = TLegend(0.65,0.7,0.89,0.89)

        if histvar == 'thetaScatt':
            leg = TLegend(0.65,0.5,0.89,0.89)
        else:
            leg = TLegend(0.35,0.2,0.65,0.5)

        hists[-1].Write()
        for h in hists[1:]:
            #    h.Draw('same')
            fout.cd()
            h.Write()

    fout.Close()

from ROOT import TRandom3

def MoliereDist(datafile, outfile, prjdatafile):
    ofile = TFile(outfile, "UPDATE")
    sdir = '/data/neutrino03/jnugent/Unfolding/Cobb_results/hack/'
    fin = open(os.path.join(sdir, datafile),'r')
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
    g.SetName(datafile[:-4])
    xbins = array('d')
    x2bins = array('d')
    data[0][3] = data[1][3]
    for d in data:
        g.SetPoint(d[0], d[1], d[2])
        xbins.append(d[1] - d[3]/2.)
        x2bins.append(d[1] + d[3]/2.)
        x2bins.insert(0,-d[1] - d[3]/2.)
    #    print d[1] - d[3]/2., d[1], d[3]

    h = TH1D("thetaScatt_Moliere",";#theta_{Scatt};Events (arbitrary normalization)",len(data)-1,xbins)

    for d in data:
        h.SetBinContent(d[0] + 1, d[2])


    c.SetLogy()

    g.SetMarkerStyle(20)
    g.SetTitle(";#theta_{Scatt};Events (arbitrary normalization)")
    g.Draw("ap")
    h.Draw("lsame")
    c.Print("Moliere_200_1.eps")



    fin = open(os.path.join(sdir, prjdatafile),'r')
    i = 0
    data = []
    last = -0.2
    for l in fin:
        ln = l.split()
        data.append([i, float(ln[0]), float(ln[1]), float(ln[0]) - last])
        last = float(ln[0])
        i+=1
    c = TCanvas()
    g = TGraph(len(data))
    g.SetName(prjdatafile[:-4])
    xbins = array('d')
    data[0][3] = data[1][3]
    for d in data:
        g.SetPoint(d[0], d[1], d[2])
        xbins.append(d[1] - d[3]/2.)
	print xbins[d[0]],d[2]
    hthetaX = TH1D("thetaX_Moliere",";#theta_{X}; Events per 2 mrad",len(xbins)-1,xbins)
    hthetaY = TH1D("thetaY_Moliere",";#theta_{Y}; Events per 2 mrad",len(xbins)-1,xbins)
    for d in data:
        hthetaX.SetBinContent(d[0] + 1, d[2])
        hthetaY.SetBinContent(d[0] + 1, d[2])
    '''
    rand = TRandom3()
    for i in range(1000000):
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
        print thetaX
	thetaY = math.atan2(X* math.sin(theta), Z)
        hthetaX.Fill(thetaX)
        hthetaY.Fill(thetaY)
    '''
    hthetaX.Draw()
    c.Print("Moliere_200_thetaX.eps")

    hthetaY.Draw()
    c.Print("Moliere_200_thetaY.eps")

    ofile.cd()
    hthetaX.Write()
    hthetaY.Write()
    h.Write()
    g.Write()
    ofile.Close()
    fin.close()

def PickUpGEANTModel(datafile, outfile):
   modfile = TFile(outfile, "UPDATE")
   modfile.cd();
   gDirectory.Delete("thetaScatt_GEANT;1")
   gDirectory.Delete("thetaX_GEANT;1")
   gDirectory.Delete("thetaY_GEANT;1")
   MCFile = TFile(datafile, "READ")
   histx = MCFile.Get("thetaX_graph")
   histy = MCFile.Get("thetaY_graph")
   histscatt = MCFile.Get("thetaScatt_graph")
   histx.SetName("thetaX_GEANT")
   histy.SetName("thetaY_GEANT")
   histscatt.SetName("thetaScatt_GEANT")
   modfile.cd()
   histx.Write()
   histy.Write()
   histscatt.Write()
   modfile.Close()
   MCFile.Close()

   '''
   modfile.cd();
   gDirectory.Delete("thetaScat_Moliere;1")
   gDirectory.Delete("thetaX_Moliere;1")
   gDirectory.Delete("thetaY_Moliere;1")
   file.Close();
   '''

if __name__ == '__main__':
    #dataModelComps('LiH-6p5-172-B-hists.dat','LiH-6p5-172-Mol-fine.root')
    #MoliereDist('LiH-6p5-172-Mol-fine.dat','LiH-6p5-172-Mol-fine.root', 'LiH-6p5-172-Molprj-fine.dat')
    #MoliereDist('LiH-6p5-200-Mol-fine.dat','LiH-6p5-200-Mol-fine.root')
    #MoliereDist('LiH-6p5-240-Mol-fine.dat','LiH-6p5-240-Mol-fine.root')
    PickUpGEANTModel('LihMuon_03172.root','LiH-6p5-172-B-hists.root')
    PickUpGEANTModel('LihMuon_03200.root','LiH-6p5-200-B-hists.root')
    PickUpGEANTModel('LihMuon_03240.root','LiH-6p5-240-hists.root')
