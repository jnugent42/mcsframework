from ROOT import TH1D, TFile, TCanvas, TLegend
import sys, math

def dataModelComps(fname, tname, outname):

    f = TFile(fname)
    fout = TFile(outname, 'RECREATE') # (tname[:-4] + '.root','RECREATE')
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
        hists = [f.Get(histvar+'_'+x) for x in histstatenames]
        hists.append(TH1D(histvar+'_Cobb',histvar,800,arraymin[q],0.2))
        if histvar=='thetaScatt':
            offset = 0
        else:
            offset = 200
        for j in range(len(lines)):
            hists[-1].SetBinContent(j+offset, float(lines[j][arrayindex[q]]))
            hists[-1].SetBinError(j+offset, math.sqrt(float(lines[j][arrayindex[q]])))
        q += 1
        # c = TCanvas(fname[:-5]+'_'+histvar+'_c1')
        print "Compatibility of Data (Deconvolved) and MC Truth",hists[0].Chi2Test(hists[1])
        print "Compatibility of Data (Deconvolved) and Cobb-Carlisle",hists[0].Chi2Test(hists[2])
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
            # h.Draw('same')
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

if __name__ == '__main__':
    dataModelComps(sys.argv[1], sys.argv[2], sys.argv[3])
