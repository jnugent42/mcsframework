from ROOT import TH1D, TFile, TCanvas, TLegend, TGraph
from array import array
import sys, math, os

def dataModelComps(fname, tname, outfile):

    f = TFile(fname)
    fout = TFile(outfile, 'RECREATE') # (tname[:-4] + '.root','RECREATE')
    # will only consider histograms with the following prefixes
    histvarnames = ['thetaX','thetaY','thetaScatt']
    arrayindex   = [2, 3, 0]
    arraymin     = [-0.2, -0.2, 0.0]
    # these histograms will have the following suffixes
    histstatenames = ['recoGEANT','GEANT']
    histstatedesc = ['Deconvolved data',
                     'GEANT4', 'Cobb-Carlisle']
    histcolors = [1, 4, 2]
    histlines = [1, 4, 8]
    histopts  = ['lp','l','l']
    # read in textfile
    lines = [line.rstrip('\n').split() for line in open(tname)]
    print lines
    # create a plot for each histvarname
    q = 0
    for histvar in histvarnames:
        print [histvar+'_'+x for x in histstatenames]
        hists = [f.Get(histvar+'_'+x) for x in histstatenames]
        hists.append(TH1D(histvar+'_Cobb',histvar,800,arraymin[q],0.2))
        # if histvar=='thetaScatt' or len(lines)==800:
        offset = 1
        #else:
        #    offset = 200
        for j in range(len(lines)):
            if histvar!="thetaScatt":
                print j, hists[-1].GetXaxis().GetBinCenter(j), float(lines[j][arrayindex[q]])
            hists[-1].SetBinContent(j+offset, float(lines[j][arrayindex[q]]))
            hists[-1].SetBinError(j+offset, math.sqrt(float(lines[j][arrayindex[q]])))
        q += 1
        # c = TCanvas(fname[:-5]+'_'+histvar+'_c1')
        # print "Compatibility of Data (Deconvolved) and MC Truth",hists[0].Chi2Test(hists[1])
        # print "Compatibility of Data (Deconvolved) and Cobb-Carlisle",hists[0].Chi2Test(hists[2])
        hists[0].SetMarkerStyle(20)
        hists[0].SetMarkerColor(histcolors[0])
        hists[0].GetXaxis().SetLabelSize(0.05)
        hists[0].GetXaxis().SetTitleSize(0.05)
        hists[0].GetYaxis().SetLabelSize(0.05)
        hists[0].GetYaxis().SetTitleSize(0.05)
        resplots = [x.Clone() for x in hists]
        resplots[0].SetTitle('')
        resplots[0].GetYaxis().SetTitle("Normalized Residuals")
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
        for i in range(len(histstatedesc)):
            # if histstatedesc[i] == 'Cobb-Carlisle':
            #     hists[i].Rebin(8)
            # else:
            #    hists[i].Rebin(2)
            hists[i].SetLineColor(histcolors[i])
            hists[i].SetLineWidth(3)
            hists[i].SetLineStyle(histlines[i])
            hists[i].SetStats(0)
            resplots[i].SetLineColor(histcolors[i])
            resplots[i].SetMarkerColor(histcolors[i])
            resplots[i].SetLineWidth(1)
            resplots[i].SetMarkerStyle(20 + histlines[i])
            resplots[i].SetStats(0)
            print hists[i].GetName(),"mean =",1000*hists[i].GetMean(),"RMS =",1000*hists[i].GetRMS(),"\pm",1000*hists[i].GetRMS()/math.sqrt(hists[i].Integral()), hists[i].Integral()
            print hists[i].Fit("gaus","0","same",-0.35,0.35)
            norm = hists[i].GetMaximum()
            # norm = hists[i].Integral()
            hists[i].Scale(1./norm)
            leg.AddEntry(hists[i], histstatedesc[i], histopts[i])

            chi2 = 0
            for j in range(0, x.GetNbinsX()+1):
                res = hists[i].GetBinContent(j) - hists[0].GetBinContent(j)
                err2 = 0
                if math.fabs(res) > 0:
                    err2 = hists[i].GetBinError(j)**2 + hists[0].GetBinError(j)**2
                if err2 > 0:
                    chi2 += res*res/err2
                    resplots[i].SetBinContent(j, res/math.sqrt(err2))
                    resplots[i].SetBinError(j, 1)
                else:
                    resplots[i].SetBinContent(j, res/1)
                    resplots[i].SetBinError(j, 0)
            print "Chi-square from residuals is ",chi2, " for " , histstatedesc[i]


        # hists[0].GetXaxis().SetRangeUser(-0.06,0.06)
        #hists[0].GetYaxis().SetRangeUser(5e-4,2)
        # hists[0].Draw()
        # print "++++++++++++++ Reco Fit ++++++++++++++++++++++"
        # f0 = hists[0].Fit("gaus","","same",-0.035,0.035)
        # print "++++++++++++++++++++++++++++++++++++++++++++++"

        # print "+++++++++++++++ MC Fit +++++++++++++++++++++++"
        # f1 = hists[1].Fit("gaus","","same",-0.035,0.035)
        # print "++++++++++++++++++++++++++++++++++++++++++++++"

        # print "+++++++++++++++ MC Fit +++++++++++++++++++++++"
        # f2 = hists[2].Fit("gaus","","same",-0.035,0.035)
        # print "++++++++++++++++++++++++++++++++++++++++++++++"

        for h in hists[1:]:
            #    h.Draw('same')
            fout.cd()
            h.Write()
        # leg.Draw('same')
        # c.SetLogy()
        # c.Print(fname[:-5]+'_'+histvar+'_Cobb.pdf')

        # c.SetLogy(0)
        # resplots[0].GetYaxis().SetRangeUser(-5,5)
        # resplots[0].Draw("p")
        # for r in resplots:
        #    r.Draw('psame')
        # leg2.Draw('same')
        # c.Print(fname[:-5]+'_'+histvar+'_res.pdf')

    fout.Close()

from ROOT import TRandom3

def MoliereDist(datafile, outfile):
    ofile = TFile(outfile, "UPDATE")
    sdir = '/data/neutrino03/jnugent/Unfolding/Cobb_results/moliere/'
    # Open Scatt text file and read
    fin = open(os.path.join(sdir, datafile),'r')
    i = 0
    data = []
    last = 0
    # Append i, bin edge, prob, bin edge - last bin edge
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
    #build thetaScatt and xbin centres
    for d in data:
        g.SetPoint(d[0], d[1], d[2])
        xbins.append(d[1] - d[3]/2.)
        x2bins.append(d[1] + d[3]/2.)
        x2bins.insert(0,-d[1] - d[3]/2.)
    #    print d[1] - d[3]/2., d[1], d[3]

    h = TH1D("thetaScatt_Moliere",";#theta_{Scatt};Events (arbitrary normalization)",len(data)-1,xbins)

    # thetaScatt is just reading in text file
    for d in data:
        h.SetBinContent(d[0] + 1, d[2])


    c.SetLogy()

    g.SetMarkerStyle(20)
    g.SetTitle(";#theta_{Scatt};Events (arbitrary normalization)")
    g.Draw("ap")
    h.Draw("lsame")
    c.Print("Moliere_200_1.eps")


    hthetaX = TH1D("thetaX_Moliere",";#theta_{X}; Events per 2 mrad",len(x2bins)-1,x2bins)
    hthetaY = TH1D("thetaY_Moliere",";#theta_{Y}; Events per 2 mrad",len(x2bins)-1,x2bins)
    rand = TRandom3()
    # thetaX & Y is sampling from thetaScatt
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
   MCfile = TFile(datafile, "READ")
   histx = MCFile.Get("theta_true_x_graph")
   histy = MCFile.Get("theta_true_y_graph")
   histscatt = MCFile.Get("theta_true_scatt_graph")
   histx.SetName("thetaX_GEANT")
   histy.SetName("thetaY_GEANT")
   histscatt.SetName("thetaScatt_GEANT")
   modfile.cd()
   histx.Write()
   histy.Write()
   histscatt.Write()
   modfile.Close()
   MCfile.Close()

   '''
   modfile.cd();
   gDirectory.Delete("thetaScat_Moliere;1")
   gDirectory.Delete("thetaX_Moliere;1")
   gDirectory.Delete("thetaY_Moliere;1")
   file.Close();
   '''

if __name__ == '__main__':
    #dataModelComps(sys.argv[1],sys.argv[2],sys.argv[3])
    #MoliereDist('LiH-6p5-172-Mol-fine.dat','LiH-6p5-172-Mol-fine.root')
    #MoliereDist('LiH-6p5-200-Mol-fine.dat','LiH-6p5-200-Mol-fine.root')
    #MoliereDist('LiH-6p5-240-Mol-fine.dat','LiH-6p5-240-Mol-fine.root')
    PickUpGEANTModel('LiH-172-MCTruth.root','LiH-6p5-172-Mol-fine.root')
    PickUpGEANTModel('LiH-200-MCTruth.root','LiH-6p5-200-Mol-fine.root')
    PickUpGEANTModel('LiH-240-MCTruth.root','LiH-6p5-172-Mol-fine.root')
