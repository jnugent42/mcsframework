from ROOT import TCanvas, TFile, TH1D, TLegend, TText, TLatex, TMath, TROOT
import sys, math
from math import sqrt
import pprint
import scipy.special as sc

TROOT.gROOT.SetBatch(1)

class plotgen:

    def __init__(self, filename):
        self.fname = filename
        self.desc = self.fname.split("_")
        self.histvarnames = ['thetaX','thetaY','thetaScatt','theta2Scatt']
        self.histstatenames = ['measdataCobb','refconv_GEANT','refconv_Cobb', 'refconv_Moliere', 'ref', 'GEANT']
        self.histstatedesc = ['Raw Data','GEANT4 Default MCS', 'Carlisle-Cobb', 'Moliere']
        self.histcolors = [1, 4, 2, 6]
        self.histopts  = ['lp','lp','lp','lp']
        self.RMS = {}
        self.RMSErr = {}
        self.Chi2 = {}
        self.pvalue = {}
        self.RMSsysdiff = {}
        self.RMSsyserr  = {}
        self.sysFiles   = []
        self.eventnorm  = {}
        self.h0integral = {}
        self.useeventnorm = True


    def formatHists(self, hist, i):
        hist.SetMarkerColor(self.histcolors[i])
        hist.GetXaxis().SetLabelSize(0.05)
        hist.GetXaxis().SetTitleSize(0.05)
        hist.GetYaxis().SetLabelSize(0.05)
        hist.GetYaxis().SetTitleSize(0.05)
        hist.SetLineColor(self.histcolors[i])
        hist.SetMarkerColor(self.histcolors[i])
        hist.SetLineWidth(3)
        if i > 0: hist.SetMarkerStyle(21+self.histcolors[i])
        else:  hist.SetMarkerStyle(20)
        hist.SetStats(0)


    def addToSysErr(self, i, histl, histh, base, scale, histvar, sysname):
        # print sysname
        # print histl
        # print histh
        norm = 0
        if i == 0:
            norm = self.eventnorm[histvar]
        else:
            norm = histl.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
        if norm>0:
            histl.Scale(1./norm)
        norm = 0
        if i == 0:
            norm = self.eventnorm[histvar]
        else:
            norm = histh.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]

        if norm>0:
            histh.Scale(1./norm)
        if histvar == 'thetaScatt':
            self.RMSsysdiff[histvar][self.histstatenames[i]][sysname] \
                    = 1000*(histh.GetMean() - histl.GetMean())
            self.RMSsyserr[histvar][self.histstatenames[i]][sysname] \
                    = 1000*math.fabs(histh.GetMean() - histl.GetMean()) * scale
            histl.GetXaxis().SetRangeUser(-0.0,0.06)
            histh.GetXaxis().SetRangeUser(-0.0,0.06)
        elif histvar == 'theta2Scatt':
            self.RMSsysdiff[histvar][self.histstatenames[i]][sysname] \
                    = 1000*(sqrt(histh.GetMean()) - \
                    sqrt(histl.GetMean()))
            self.RMSsyserr[histvar][self.histstatenames[i]][sysname] \
                            = 1000*math.fabs(sqrt(histh.GetMean()) - \
                            sqrt(histl.GetMean())) * scale
            histl.GetXaxis().SetRangeUser(-0.0,0.0036)
            histh.GetXaxis().SetRangeUser(-0.0,0.0036)
        else:
            histl.GetXaxis().SetRangeUser(-0.06,0.06)
            histh.GetXaxis().SetRangeUser(-0.06,0.06)
            histl.Fit("gaus","Q0","same",-0.045,0.045)
            histh.Fit("gaus","Q0","same",-0.045,0.045)
            self.RMSsysdiff[histvar][self.histstatenames[i]][sysname] \
                = 1000*(histh.GetFunction('gaus').GetParameter('Sigma') - \
                histl.GetFunction('gaus').GetParameter('Sigma'))
            self.RMSsyserr[histvar][self.histstatenames[i]][sysname] \
                = math.fabs(self.RMSsysdiff[histvar][self.histstatenames[i]][sysname]) * scale

    def addToRMS(self, i, hist, base, resplot, histvar):
        norm = 0
        if i == 0:
            if self.useeventnorm:
                self.eventnorm[histvar]  = hist.GetEntries()
            self.h0integral[histvar] = hist.Integral()
            norm = self.eventnorm[histvar]
        else:
            norm = hist.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
        hist.Scale(1./norm)
        if histvar == 'thetaScatt':
            self.RMS[histvar][self.histstatenames[i]]    = (1000*hist.GetMean())
            self.RMSErr[histvar][self.histstatenames[i]] = (1000*hist.GetMeanError())
            hist.GetXaxis().SetRangeUser(-0.0,0.06)
            resplot.GetXaxis().SetRangeUser(-0.0,0.06)
        elif histvar == 'theta2Scatt':
            hist.Fit("expo","Q0","same",-0.0,0.002)
            self.RMS[histvar][self.histstatenames[i]]    = (1000/sqrt(-hist.GetFunction('expo').GetParameter('Slope')))
            self.RMSErr[histvar][self.histstatenames[i]] = (1000*hist.GetFunction('expo').GetParError(1)
                    /(-hist.GetFunction('expo').GetParameter('Slope'))**(3./2.))
            hist.GetXaxis().SetRangeUser(-0.0,0.0036)
            resplot.GetXaxis().SetRangeUser(-0.0,0.0036)
        else:
            hist.GetXaxis().SetRangeUser(-0.06,0.06)
            resplot.GetXaxis().SetRangeUser(-0.06,0.06)
            hist.Fit("gaus","Q0","same",-0.045,0.045)
            self.RMS[histvar][self.histstatenames[i]]    = (1000*hist.GetFunction('gaus').GetParameter('Sigma'))
            self.RMSErr[histvar][self.histstatenames[i]] = (1000*hist.GetFunction('gaus').GetParError(2))

    def calculateChi2(self, i, hist, base, resplot, histvar, pname):

        chi2 = 0
        cv = 0
        ndf = 0
        sysLowF  = [TFile(x[0]) for x in self.sysFiles]
        sysLowH  = [f.Get(histvar + '_' + self.histstatenames[i]) for f in sysLowF]

        sysHighF = [TFile(x[1]) for x in self.sysFiles]
        sysHighH = [f.Get(histvar + '_' + self.histstatenames[i]) for f in sysHighF]
        sysScale = [x[2] for x in self.sysFiles]
        sysName = [x[3] for x in self.sysFiles]
        syshists = []
        for q in range(len(self.sysFiles)):
            syshists.append(resplot.Clone())
            syshists[-1].SetName(resplot.GetName()+"_"+self.sysFiles[q][3])
            syshists[-1].GetYaxis().SetTitle("Upper - Lower Systematic")


        self.RMSsysdiff[histvar][self.histstatenames[i]] = {}
        self.RMSsyserr[histvar][self.histstatenames[i]]  = {}
        self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum2"] = 0
        self.RMSsyserr[histvar][self.histstatenames[i]]["Sum2"] = 0
        self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum"] = 0
        self.RMSsyserr[histvar][self.histstatenames[i]]["Sum"] = 0
        binmax = 0.045
        cloneb = TH1D("","",29,-binmax,binmax)
        cloneh = TH1D("","",29,-binmax,binmax)
        z = 1
        cv = 0
        if histvar == 'theta2Scatt':
            binmax = 0.0036
        for k in range(len(self.sysFiles)):
            if pname=="Gold" and i==1:
                continue

            self.formatHists(sysLowH[k], i)
            self.formatHists(sysHighH[k], i)
            # print self.sysFiles[k][0]
            self.addToSysErr(i, sysHighH[k], sysLowH[k], base, sysScale[k], histvar, self.sysFiles[k][3])
            self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum2"] \
                        += self.RMSsysdiff[histvar][self.histstatenames[i]][self.sysFiles[k][3]] *\
                        self.RMSsysdiff[histvar][self.histstatenames[i]][self.sysFiles[k][3]]
            self.RMSsyserr[histvar][self.histstatenames[i]]["Sum2"] \
                        += self.RMSsyserr[histvar][self.histstatenames[i]][self.sysFiles[k][3]]  *\
                        self.RMSsyserr[histvar][self.histstatenames[i]][self.sysFiles[k][3]]
            self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum"] \
                        = sqrt( self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum2"] )
            self.RMSsyserr[histvar][self.histstatenames[i]]["Sum"] \
                        = sqrt( self.RMSsyserr[histvar][self.histstatenames[i]]["Sum2"] )

        for j in range(0, base.GetNbinsX()+1):
            res = hist.GetBinContent(j) - base.GetBinContent(j)
            # resplot.SetBinContent(j, res/math.sqrt(err2hist))
            # resplot.SetBinError(j, 1)
            if pname == "Result" and i == 2:
                res = hist.GetBinContent(hist.GetMaximumBin()-24+j) - base.GetBinContent(j)
                # resplot.SetBinContent(hist.GetMaximumBin()-24+j, res/math.sqrt(err2hist))
            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax :
                cloneh.SetBinContent(z,hist.GetBinContent(j))
                cloneb.SetBinContent(z,base.GetBinContent(j))
                if hist.GetBinContent(j) > 0:
                    cv += res*res/hist.GetBinContent(j)
                z += 1
            if pname == "Result" and i == 2:
                if math.fabs(hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()-15+j)) < binmax -0.001:
                    cloneh.SetBinContent(z,hist.GetBinContent(hist.GetMaximumBin()-15+j))
                    cloneb.SetBinContent(z,base.GetBinContent(j))
                    if hist.GetBinContent(hist.GetMaximumBin()-15+j) > 0:
                        cv += res*res/hist.GetBinContent(hist.GetMaximumBin()-15+j)
                    z += 1
            err2hist = 0
            err2hist = base.GetBinError(j)**2
            # if pname == "Result" and i == 2:
            #     j = hist.GetMaximumBin()-24+j

            err2 = 0
            err2 = hist.GetBinError(j)**2

            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < 0.7:
                for k in range(len(self.sysFiles)):
                        sys = sysHighH[k].GetBinContent(j) - sysLowH[k].GetBinContent(j)
                        sys *= sysScale[k]
                        syshists[k].SetBinContent(j, sys*sys)
                        syshists[k].SetBinError(j, 0)
                        err2 += sys*sys
                        err2hist += sys*sys

                        if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax and err2hist > 0 :
                            chi2 += res*res/err2hist
                        if pname == "Result" and i == 2:
                            if math.fabs(hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()-24+j)) < binmax - 0.001 and err2hist > 0 :
                                chi2 += res*res/err2hist

                hist.SetBinError(j, sqrt(err2))

            if err2hist > 0:
                resplot.SetBinContent(j, res/math.sqrt(err2hist))
            resplot.SetBinError(j, 1)
            if pname == "Result" and i == 2:
                res = hist.GetBinContent(hist.GetMaximumBin()-24+j) - base.GetBinContent(j)
                if err2hist > 0:
                    resplot.SetBinContent(hist.GetMaximumBin()-24+j, res/math.sqrt(err2hist))
                resplot.SetBinError(hist.GetMaximumBin()-24+j, 1)
            if pname == "Result" and i == 2:
                res = hist.GetBinContent(j) - base.GetBinContent(j)
                if err2hist > 0:
                    resplot.SetBinContent(j, res/math.sqrt(err2hist))
                resplot.SetBinError(j, 1)
            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax:
                ndf += 1
            if pname == "Result" and i == 2:
                if math.fabs(hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()-24+j)) < binmax - 0.001:
                    # print "here",hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()-24+j)
                    ndf += 1
            #print i, resplot.GetXaxis().GetBinCenter(j), res, chi2, err2, ndf
        c = TCanvas(self.fname[:-5]+'_'+histvar+'_c1')
        c.SetBottomMargin(0.15)
        c.SetTopMargin(0.075)

        for k in range(len(self.sysFiles)):
            c.SetLogy()
            syshists[k].Draw("p")
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_'+self.histstatenames[i]+'_'+self.sysFiles[k][3]+'_res.pdf')
        #print "Chi-square from residuals is ",chi2, " for " ,ndf," from ", self.histstatedesc[i], binmaxa
        # PValue_root = TMath.Prob(ndf/2,chi2/2)

        # my own p-value calculation
        # ndf = 100
        # chi2 = 100
        # K = ndf * 0.5
        # X = cv * 0.5
        # PValue = sc.gammainc(0, 29)
        # if(isnan(PValue) or isinf(PValue) or PValue <= 1e-8):
        #     return 1e-14
        # PValue = PValue1/math.gamma(K)
        # PValue = (1.0 - PValue)

        histclone = hist.Clone()
        baseclone = base.Clone()
        histclone.GetXaxis().SetRangeUser(-binmax,binmax)
        baseclone.GetXaxis().SetRangeUser(-binmax,binmax)
        histclone2 = histclone.Clone()
        baseclone2 = baseclone.Clone()
        # pvalue_root =  baseclone.Chi2Test(histclone,"WW")
        # kol =  baseclone.KolmogorovTest(histclone)
        # chi2_root =  baseclone.Chi2Test(histclone,"CHI2 WW")
        print pname
        print i
        print "final cv",cv
        # pvalue_root =  baseclone2.Chi2Test(histclone2,"WW")
        # kol =  baseclone2.KolmogorovTest(histclone2)
        # chi2_root =  baseclone2.Chi2Test(histclone2,"CHI2 WW")

        pvalue_root =  cloneb.Chi2Test(cloneh,"")
        kol =  hist.KolmogorovTest(base)
        if pname == "Result" and i == 2:
            kol =  cloneb.KolmogorovTest(cloneh)
        chi2_root =  cloneb.Chi2Test(cloneh,"CHI2 ")
        # pvalue_root =  histclone.Chi2Test(histclone, "WW")
        # chi2_root =  histclone.Chi2Test(histclone,"CHI2 WW")
        # kol =  histclone.KolmogorovTest(histclone)
        ndfr = 1
        # if (histclone.Chi2Test(baseclone,"WW CHI2/NDF")!=0 ):
        #     ndfr = 1/histclone.Chi2Test(baseclone,"WW CHI2/NDF")
        # PValue1 = sc.gammainc(chi2_root/2, ndfr*chi2_root/2)
        # PValue = PValue1/math.gamma(29)


        PValue1 = sc.gammainc(chi2/2, ndf/2)
        PValue2 = sc.gammaincc(chi2/2, ndf/2)
        # PValue2 = sc.gammaincc(ndf/2, chi2/2)
        # PValue = PValue1/math.gamma(ndf)
        PValue = (1.0 - PValue2)
        self.Chi2[histvar][self.histstatenames[i]] = [chi2, ndf]
        # self.Chi2[histvar][self.histstatenames[i]] = [chi2_root, ndf]
        # self.pvalue[histvar][self.histstatenames[i]] = [kol]
        # self.pvalue[histvar][self.histstatenames[i]] = [pvalue_root]
        self.pvalue[histvar][self.histstatenames[i]] = [PValue, kol]

    def MCSPlot(self, pname):

        f = TFile(self.fname)
        self.RMS = {}
        self.RMSErr = {}
        self.Chi2 = {}
        self.RMSsysdiff = {}
        self.RMSsyserr  = {}

        # create a plot for each histvarname
        for histvar in self.histvarnames:
            self.RMS[histvar]     = {}
            self.RMSErr[histvar]  = {}
            self.RMSsysdiff[histvar] = {}
            self.RMSsyserr[histvar]  = {}
            self.Chi2[histvar]    = {}
            self.pvalue[histvar]    = {}
            names = [histvar + '_' + x for x in self.histstatenames]
            hists = [f.Get(histvar+'_'+x) for x in self.histstatenames]

            hists[0].SetTitle("")
            resplots = [x.Clone() for x in hists]
            resplots[0].SetTitle('')
            resplots[0].GetYaxis().SetTitle("Normalized Residuals")

            # leg = TLegend(0.55,0.73,0.89,0.92)
            leg = TLegend(0.25,0.23,0.59,0.42)
            leg.SetLineColor(10)
            for i in range(len(self.histstatedesc)):
                hists[i].Sumw2()
                self.formatHists(hists[i], i)
                self.formatHists(resplots[i], i)
                self.addToRMS(i, hists[i], hists[0], resplots[i], histvar)
                if histvar=='theta2Scatt':
                    hists[i].GetYaxis().SetTitle('Probability per '+str("%.2f" % round(1000*1000*hists[i].GetXaxis().GetBinWidth(4),2))+' mrad^{2}')
                else:
                    hists[i].GetYaxis().SetTitle('Probability per '+str("%.2f" % round(1000*hists[i].GetXaxis().GetBinWidth(4),2))+' mrad')
                leg.AddEntry(hists[i], self.histstatedesc[i], self.histopts[i])
                self.calculateChi2(i, hists[i], hists[0], resplots[i], histvar, pname)


            c = TCanvas(self.fname[:-5]+'_'+histvar+'_c1')
            if self.desc[0] == 'XePion':
                t1 = TText(0.18,0.885,"MICE ISIS cycle 2015/03")
                t2 = TText(0.18,0.85,"Xe, "+self.desc[1][2:5]+", MAUS v3.1.2")
            else:
                t1 = TText(0.18,0.885,"MICE ISIS cycle 2015/04")
                t2 = TText(0.18,0.85,"LiH, "+self.desc[1][2:5]+", MAUS v3.3.2")
            t1.SetNDC(1)
            t1.SetTextSize(0.04)
            t1.SetTextFont(42)
            t2.SetNDC(1)
            t2.SetTextSize(0.03)
            t2.SetTextFont(42)
            # hists[0].GetYaxis().SetRangeUser(4e-5,2.0)
            hists[0].SetTitle(";"+hists[0].GetXaxis().GetTitle()+" (radians);"+hists[0].GetYaxis().GetTitle())

            histsminus0 = TH1D("",";;asymmetry",47,-0.06,0.06)
            for j in range(0, hists[0].GetNbinsX()+1):
                histsminus0.SetBinContent(j,hists[0].GetBinContent(47-j))
            Truthasymm = hists[0].GetAsymmetry(histsminus0);
            Truthasymm.Draw('hist p')
            c.SaveAs('Truthasymm.pdf')
            c.Clear()
            histsminus1 = TH1D("",";;asymmetry",47,-0.06,0.06)
            for j in range(0, hists[1].GetNbinsX()+1):
                histsminus1.SetBinContent(j,hists[1].GetBinContent(47-j))
            Truthasymm1 = hists[1].GetAsymmetry(histsminus1);
            Truthasymm1.Draw('hist p')
            c.SaveAs('Truthasymm1.pdf')
            c.Clear()
            hists[0].Draw('ep')
            leg.SetTextSize(0.04)
            leg.Draw('same')
            c.SetBottomMargin(0.15)
            c.SetTopMargin(0.075)
            # for h in hists[1:len(self.histstatedesc)]:
            for h in hists:
                h.Draw('same')
                t1.Draw()
                t2.Draw()
                c.SetLogy()
                c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys.pdf')
                # c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys.root')
                # c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_pq.jpg')
                # leg.SetX1NDC(0.5)
                # leg.SetX2NDC(0.89)
                # leg.SetY1NDC(0.2)
                # leg.SetY2NDC(0.4)
                # leg.SetTextSize(0.04)
                # leg.Draw('same')
            c.Clear()
            c.SetLogy(0)
            resplots[0].GetYaxis().SetRangeUser(-2,2)
            resplots[0].SetTitle(";"+resplots[0].GetXaxis().GetTitle()+" (radians);"+resplots[0].GetYaxis().GetTitle())
            leg.SetX1NDC(0.5)
            leg.SetX2NDC(0.89)
            leg.SetY1NDC(0.2)
            leg.SetY2NDC(0.4)
            resplots[0].Draw("p")
            for r in resplots:
                r.GetYaxis().SetRangeUser(-5,5)
                r.Draw('psame')
            leg.SetTextSize(0.04)
            leg.Draw('same')
            t1.Draw()
            t2.Draw()
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_res_T.pdf')
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_res_pq.jpg')

        momhist = f.Get("cor_mom")
        if self.fname.find("LihMuon_03172")  >= 0:
            mom = [momhist.GetMean(), momhist.GetMeanError()]
        elif self.fname.find("LihMuon_03200")  >= 0:
            mom = [momhist.GetMean(), momhist.GetMeanError()]
        elif self.fname.find("LihMuon_03240")  >= 0:
            mom = [momhist.GetMean(), momhist.GetMeanError()]
        rms = [momhist.GetRMS(), momhist.GetRMSError()]
        summary = []
        syssummary = []
        def sigfig(x):
            if math.fabs(x) > 1e-5:
                return int(math.ceil(math.fabs(math.log(math.fabs(x),10))))
            else:
                return 1


        if pname!="Gold":

         for sys in self.sysFiles:
            stindx = 0

            difference0 = self.RMSsysdiff[self.histvarnames[0]][self.histstatenames[stindx]][sys[3]]
            difference1 = self.RMSsysdiff[self.histvarnames[1]][self.histstatenames[stindx]][sys[3]]
            # difference3 = self.RMSsysdiff[self.histvarnames[3]][self.histstatenames[stindx]][sys[3]]

            syserr0 = self.RMSsyserr[self.histvarnames[0]][self.histstatenames[stindx]][sys[3]]
            syserr1 = self.RMSsyserr[self.histvarnames[1]][self.histstatenames[stindx]][sys[3]]
            # syserr3 = self.RMSsyserr[self.histvarnames[3]][self.histstatenames[stindx]][sys[3]]

            relerr0 = syserr0/self.RMS[self.histvarnames[0]][self.histstatenames[0]]
            relerr1 = syserr1/self.RMS[self.histvarnames[1]][self.histstatenames[0]]
            # relerr3 = syserr3/self.RMS[self.histvarnames[3]][self.histstatenames[0]]

            syssummary.append("LiH & "+str("%.2f" % round(mom[0],sigfig(mom[1])))+"$\pm$"+str("%.2f" % round(mom[1],sigfig(mom[1])))+\
                    " & "+str("%.2f" % round(difference0,sigfig(difference0)))+\
                    " & "+str("%.2f" % round(syserr0,sigfig(syserr0)))+\
                    " & "+str("%.2f" % round(relerr0,sigfig(relerr0)))+\
                    " & "+str("%.2f" % round(difference1,sigfig(difference1)))+\
                    " & "+str("%.2f" % round(syserr1,sigfig(syserr1)))+\
                    " & "+str("%.2f" % round(relerr1,sigfig(relerr1)))+\
                    # " & "+str("%.2f" % round(difference3,sigfig(difference3)))+\
                    # " & "+str("%.2f" % round(syserr3,sigfig(syserr3)))+\
                    # " & "+str("%.2f" % round(relerr3,sigfig(relerr3)))
                    "\\\\")

         syssummary.append(str("%.2f" % round(mom[0],2))+"$\pm$"+str("%.2f" % round(mom[1],2))+\
                " & "+ str("%.2f" % round(rms[0],2))+"$\pm$"+str("%.2f" % round(rms[1],2))+\
                " & "+str("%.2f" % round(self.RMSsysdiff[self.histvarnames[0]][self.histstatenames[0]]['Sum'],2))+\
                " & "+str("%.2f" % round(self.RMSsyserr[self.histvarnames[0]][self.histstatenames[0]]['Sum'],2))+\
                " & "+str("%.2f" % round(self.RMSsyserr[self.histvarnames[0]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[0]][self.histstatenames[0]],2))+\
                " & "+str("%.2f" % round(self.RMSsysdiff[self.histvarnames[1]][self.histstatenames[0]]['Sum'],2))+\
                " & "+str("%.2f" % round(self.RMSsyserr[self.histvarnames[1]][self.histstatenames[0]]['Sum'],2))+\
                " & "+str("%.2f" % round(self.RMSsyserr[self.histvarnames[1]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[1]][self.histstatenames[0]],2))+\
                # " & "+str("%.2f" % round(self.RMSsysdiff[self.histvarnames[3]][self.histstatenames[0]]['Sum'],2))+
                # " & "+str("%.2f" % round(self.RMSsyserr[self.histvarnames[3]][self.histstatenames[0]]['Sum'],2))+
                # " & "+str("%.2f" % round(self.RMSsyserr[self.histvarnames[3]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[2]][self.histstatenames[0]],2))
                "\\\\")


        # if pname == 'Result' or pname == 'con' or pname == 'Gold':
        # if pname == 'Gold':
        #  for histvar in self.histvarnames:
        #      summary.append(str("%.2f" % round(mom[0],2))+"$\pm$"+str("%.2f" % round(mom[1],2))+\
        #              "& $\ "+histvar+"$ & "+str("%.2f" % round(self.RMS[histvar][self.histstatenames[0]],2))+ \
        #              "$\pm$"+str("%.2f" % round(self.RMSErr[histvar][self.histstatenames[0]],2))+ \
        #              "$\pm$"+str("%.2f" % round(self.RMSsyserr[histvar][self.histstatenames[0]]["Sum"],2))+ \
        #              " & "+str("%.2f" % round(self.RMS[histvar][self.histstatenames[1]],2))+ \
        #              "$\pm$"+str("%.2f" % round(self.RMSErr[histvar][self.histstatenames[1]],2))+ \
        #              " & "+str("%.2f" % round(self.Chi2[histvar][self.histstatenames[1]][0],1))+ \
        #              " / "+ str(self.Chi2[histvar][self.histstatenames[1]][1]) \
        #              +"\\\\")
        # else:
        #      for histvar in self.histvarnames:
        #          summary.append(str("%.2f" % round(mom[0],2))+"$\pm$"+str("%.2f" % round(mom[1],2))+\
        #                  "& $\\"+histvar+"$ & "+str("%.2f" % round(self.RMS[histvar][self.histstatenames[0]],2))+ \
        #                   "$\pm$"+str("%.2f" % round(self.RMSErr[histvar][self.histstatenames[0]],2))+ \
        #                   "$\pm$"+str("%.2f" % round(self.RMSsyserr[histvar][self.histstatenames[0]]["Sum"],2))+ \
        #                   " & "+str("%.2f" % round(self.RMS[histvar][self.histstatenames[1]],2))+ \
        #                   "$\pm$"+str("%.2f" % round(self.RMSErr[histvar][self.histstatenames[1]],2))+ \
        #                   # " & "+str("%.2f" % round(self.Chi2[histvar][self.histstatenames[1]][0],1))+ \
        #                   " & "+str(self.Chi2[histvar][self.histstatenames[1]][0])+ \
        #                   " / "+ str(self.Chi2[histvar][self.histstatenames[1]][1])+ \
        #                   # " & "+str("%.2f" % round(self.pvalue[histvar][self.histstatenames[1]][0],1))+ \
        #                   " & "+str(self.pvalue[histvar][self.histstatenames[1]][0])+ \
        #                   " & "+str(self.pvalue[histvar][self.histstatenames[1]][1])+ \
        #                   " & "+str("%.2f" % round(self.RMS[histvar][self.histstatenames[2]],2))+ \
        #                  "$\pm$"+str("%.2f" % round(self.RMSErr[histvar][self.histstatenames[2]],2))+ \
        #                  # " & "+str("%.2f" % round(self.Chi2[histvar][self.histstatenames[2]][0],1))+ \
        #                  " & "+str(self.Chi2[histvar][self.histstatenames[2]][0])+ \
        #                  " / "+ str(self.Chi2[histvar][self.histstatenames[2]][1])+ \
        #                  # " & "+str("{:e}".format(round(self.pvalue[histvar][self.histstatenames[2]][0],1)))+ \
        #                  " & "+str(self.pvalue[histvar][self.histstatenames[2]][0])+ \
        #                  " & "+str(self.pvalue[histvar][self.histstatenames[2]][1])+ \
        #                  "\\\\")
        if pname == 'Gold':
         for histvar in self.histvarnames:
             summary.append(str("%.2f" % mom[0])+"$\pm$"+str("%.2f" % mom[1])+\
                     "& $\ "+histvar+"$ & "+str("%.2f" % self.RMS[histvar][self.histstatenames[0]])+ \
                     "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[0]])+ \
                     "$\pm$"+str("%.2f" % self.RMSsyserr[histvar][self.histstatenames[0]]["Sum"])+ \
                     " & "+str("%.2f" % self.RMS[histvar][self.histstatenames[1]])+ \
                     "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[1]])+ \
                     " & "+str("%.2f" % self.Chi2[histvar][self.histstatenames[1]][0],1)+ \
                     " / "+ str(self.Chi2[histvar][self.histstatenames[1]][1]) \
                     +"\\\\")
        else:
             for histvar in self.histvarnames:
                 summary.append(str("%.2f" % mom[0])+"$\pm$"+str("%.2f" % mom[1])+\
                         "& $\\"+histvar+"$ & "+str("%.2f" % self.RMS[histvar][self.histstatenames[0]])+ \
                          "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[0]])+ \
                          "$\pm$"+str("%.2f" % self.RMSsyserr[histvar][self.histstatenames[0]]["Sum"])+ \
                          " & "+str("%.2f" % self.RMS[histvar][self.histstatenames[1]])+ \
                          "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[1]])+ \
                          " & "+str("%.2f" % self.Chi2[histvar][self.histstatenames[1]][0])+ \
                          " / "+ str(self.Chi2[histvar][self.histstatenames[1]][1])+ \
                          # " & "+str("%.2f" % self.pvalue[histvar][self.histstatenames[1]][0],1))+ \
                          " & "+str("{:.2e}".format(self.pvalue[histvar][self.histstatenames[1]][0]))+ \
                          " & "+str("%.2f" % self.pvalue[histvar][self.histstatenames[1]][1])+ \
                          " & "+str("%.2f" % self.RMS[histvar][self.histstatenames[2]])+ \
                         "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[2]])+ \
                         " & "+str("%.2f" % self.Chi2[histvar][self.histstatenames[2]][0])+ \
                         " / "+ str(self.Chi2[histvar][self.histstatenames[2]][1])+ \
                         # " & "+str("{:e}".format(self.pvalue[histvar][self.histstatenames[2]][0],1)))+ \
                         " & "+str("{:.2e}".format(self.pvalue[histvar][self.histstatenames[2]][0]))+ \
                         " & "+str("%.2f" % self.pvalue[histvar][self.histstatenames[2]][1])+ \
                         "\\\\")

        f.Close()
        return [summary, syssummary]

if __name__ == '__main__':

    raw172 = plotgen("LihMuon_03172.root")
    raw172.sysFiles = [
                      ["../TOFsys/LiHMu_3200_tof_lim28.5537703996_u28.7537703996.root",
                        "../TOFsys/LiHMu_3200_tof_lim28.7537703996_u28.9537703996.root",
                        70./200.,"TOF"],
                      ["../angdef/angdef_0/LiHMu_3172_15.root",
                       "../angdef/angdef_0/LiHMu_3172_100.root",
                       1,"angdef"],
                      ["../MClih172/LihMuon_03172.root",
                       "../MClih172/nopionLihMuon_03172.root",
                       1,"picon"],
                      ["../fiducial/LiHMu_3172_R_-1_G_1.root",
                       "../fiducial/LiHMu_3172_R_1_G_1.root",
                      .478/20.,"Fid. Radius"],
                      ["../alignment/LiHMu_3172_-1.root",
                       "../alignment/LiHMu_3172_1.root",
                       8.9e-5/8.1e-4,"Alignment"]
                      ]

    raw200 = plotgen("LihMuon_03200.root")
    raw200.sysFiles = [
                       ["../TOFsys/LiHMu_3200_tof_lim27.7537703996_u27.9537703996.root",
                        "../TOFsys/LiHMu_3200_tof_lim28.1537703996_u28.3537703996.root",
                        70./400.,"TOF"],
                       ["../angdef/angdef_0/LiHMu_3200_15.root",
                        "../angdef/angdef_0/LiHMu_3200_100.root",
                        1,"angdef"],
                        ["../MClih200/LihMuon_03200.root",
                        "../MClih200/nopionLihMuon_03200.root",
                        1,"picon"],
                       ["../fiducial/LiHMu_3200_R_-1_G_1.root",
                        "../fiducial/LiHMu_3200_R_1_G_1.root",
                        0.478/20.,"Fid. Radius"],
                       ["../alignment/LiHMu_3200_-1.root",
                        "../alignment/LiHMu_3200_1.root",
                        2.1e-5/9.1e-5,"Alignment"]
                       ]

    raw240 = plotgen("LihMuon_03240.root")
    raw240.sysFiles = [
                       ["../TOFsys/LiHMu_3200_tof_lim27.1537703996_u27.3537703996.root",
                        "../TOFsys/LiHMu_3200_tof_lim27.3537703996_u27.5537703996.root",
                        70./200.,"TOF"],
                       ["../angdef/angdef_0/LiHMu_3240_15.root",
                        "../angdef/angdef_0/LiHMu_3240_100.root",
                        1,"angdef"],
                       ["../MClih240/LihMuon_03240.root",
                        "../MClih240/nopionLihMuon_03240.root",
                        1,"picon"],
                       ["../fiducial/LiHMu_3240_R_-1_G_1.root",
                        "../fiducial/LiHMu_3240_R_1_G_1.root",
                        0.478/20.,"Fid. Radius"],
                       ["../alignment/LiHMu_3240_-1.root",
                        "../alignment/LiHMu_3240_1.root",
                        1.2e-5/6.5e-5,"Alignment"]
                       ]

    [out172, sys172] = raw172.MCSPlot("raw")
    [out200, sys200] = raw200.MCSPlot("raw")
    [out240, sys240] = raw240.MCSPlot("raw")

    print "\n"
    raw172.histstatenames = ['measdataGEANT','ref','refconv_GEANT']
    raw172.histstatedesc = ['Raw LiH Data','Empty AFC', 'GEANT Model conv. Empty']
    raw200.histstatenames = ['measdataGEANT','ref','refconv_GEANT']
    raw200.histstatedesc = ['Raw LiH Data','Empty AFC', 'GEANT Model conv Empty']
    raw240.histstatenames = ['measdataGEANT','ref','refconv_GEANT']
    raw240.histstatedesc = ['Raw LiH Data','Empty AFC', 'GEANT Model conv Empty']
    [ref172, refsys172] = raw172.MCSPlot("refG4")
    [ref200, refsys200] = raw200.MCSPlot("refG4")
    [ref240, refsys240] = raw240.MCSPlot("refG4")

    print "\n"
    raw172.histstatenames = ['measdataGEANT','refconv_GEANT','refconv_Moliere']
    raw172.histstatedesc = ['Raw LiH Data', 'GEANT Model conv. Empty','Moliere Model conv. Empty']
    raw200.histstatenames = ['measdataGEANT','refconv_GEANT','refconv_Moliere']
    raw200.histstatedesc = ['Raw LiH Data', 'GEANT Model conv. Empty','Moliere Model conv. Empty']
    raw240.histstatenames = ['measdataGEANT','refconv_GEANT','refconv_Moliere']
    raw240.histstatedesc = ['Raw LiH Data', 'GEANT Model conv. Empty','Moliere Model conv. Empty']
    [con172, consys172] = raw172.MCSPlot("con")
    [con200, consys200] = raw200.MCSPlot("con")
    [con240, consys240] = raw240.MCSPlot("con")

    raw172.histvarnames = ['thetaX','thetaY']
    raw200.histvarnames = ['thetaX','thetaY']
    raw240.histvarnames = ['thetaX','thetaY']
    raw172.histstatenames = ['recoGold_effi_only_gold', 'GEANT', 'Moliere']
    raw172.histstatedesc = ['Deconvolved Gold', 'GEANT', 'Moliere Model']
    raw200.histstatenames = ['recoGold_effi_only_gold', 'GEANT', 'Moliere']
    raw200.histstatedesc = ['Deconvolved Gold', 'GEANT', 'Moliere Model']
    raw240.histstatenames = ['recoGold_effi_only_gold', 'GEANT', 'Moliere']
    raw240.histstatedesc = ['Deconvolved Gold', 'GEANT', 'Moliere Model']

    [Gold172, Goldsys172] = raw172.MCSPlot("Result")
    [Gold200, Goldsys200] = raw200.MCSPlot("Result")
    [Gold240, Goldsys240] = raw240.MCSPlot("Result")

    print 'Gold sys'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3]
        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
        print Goldsys172[i]
        print Goldsys200[i]
        print Goldsys240[i]

    print 'Gold sys summary'
    print '\hline\n'
    print '\hline'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3], Goldsys200[i]
    print '\hline\n'
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print Goldsys200[-1]
    print '\hline\n'
    print '\hline'

    print '\hline'
    print "Gold deco table"
    print '\hline'
    print '\hline'
    print Gold172[0]
    print Gold172[1]
    print '\hline'
    print Gold200[0]
    print Gold200[1]
    print '\hline'
    print Gold240[0]
    print Gold240[1]

    print '\hline'
    print '\hline'
    # print Gold172[3]
    # print Gold200[3]
    # print Gold240[3]

    print '\hline'
    print '\hline\n'
    print "Raw deco table"
    print '\hline'
    print '\hline'
    print out172[0]
    print out172[1]
    print '\hline'
    print out200[0]
    print out200[1]
    print '\hline'
    print out240[0]
    print out240[1]

    print '\hline'
    print '\hline'
    print out172[3]
    print out200[3]
    print out240[3]
    print '\hline'
    print '\hline\n'

    print 'con'
    print '\hline'
    print '\hline'
    print con172[0]
    print con172[1]
    print '\hline'
    print con200[0]
    print con200[1]
    print '\hline'
    print con240[0]
    print con240[1]

    print '\hline'
    print '\hline'
    print con172[3]
    print con200[3]
    print con240[3]
    print '\hline'
    print '\hline\n'

    print 'raw sys'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3]
        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
        print sys172[i]
        print sys200[i]
        print sys240[i]
        print '\hline\n'
    print '\hline'

    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print sys172[-1]
    print sys200[-1]
    print sys240[-1]
    print '\hline\n'
    print '\hline'

    Goldfile172 = plotgen("MCLihMuon_03172.root")
    Goldfile172.sysFiles = [["../TOFsys/LiHMu_3200_tof_lim28.5537703996_u28.7537703996.root",
                        "../TOFsys/LiHMu_3200_tof_lim28.7537703996_u28.9537703996.root",
                        70./200.,"TOF"],
                        ["../angdef/angdef_0/LiHMu_3172_15.root",
                        "../angdef/angdef_0/LiHMu_3172_100.root",
                        1,"angdef"],
                        ["../MClih172/LihMuon_03172.root",
                        "../MClih172/nopionLihMuon_03172.root",
                        1,"picon"],
                        ["../fiducial/LiHMu_3172_R_-1_G_1.root",
                        "../fiducial/LiHMu_3172_R_1_G_1.root",
                        .478/20.,"Fid. Radius"],
                        ["../alignment/LiHMu_3172_-1.root",
                         "../alignment/LiHMu_3172_1.root",
                         8.9e-5/8.1e-4,"Alignment"]
                        ]

    Goldfile200 = plotgen("MCLihMuon_03200.root")
    Goldfile200.sysFiles = [
                        ["../TOFsys/LiHMu_3200_tof_lim27.7537703996_u27.9537703996.root",
                        "../TOFsys/LiHMu_3200_tof_lim28.1537703996_u28.3537703996.root",
                        70./400.,"TOF"],
                        ["../angdef/angdef_0/LiHMu_3200_15.root",
                        "../angdef/angdef_0/LiHMu_3200_100.root",
                        1,"angdef"],
                        ["../MClih200/LihMuon_03200.root",
                        "../MClih200/nopionLihMuon_03200.root",
                        1,"picon"],
                        ["../fiducial/LiHMu_3200_R_-1_G_1.root",
                        "../fiducial/LiHMu_3200_R_1_G_1.root",
                        0.478/20.,"Fid. Radius"],
                        ["../alignment/LiHMu_3200_-1.root",
                        "../alignment/LiHMu_3200_1.root",
                        2.1e-5/9.1e-5,"Alignment"]
                        ]

    Goldfile240 = plotgen("MCLihMuon_03240.root")
    Goldfile240.sysFiles = [["../TOFsys/LiHMu_3200_tof_lim27.1537703996_u27.3537703996.root",
                        "../TOFsys/LiHMu_3200_tof_lim27.3537703996_u27.5537703996.root",
                        70./200.,"TOF"],
                        ["../angdef/angdef_0/LiHMu_3240_15.root",
                        "../angdef/angdef_0/LiHMu_3240_100.root",
                        1,"angdef"],
                        ["../MClih240/LihMuon_03240.root",
                        "../MClih240/nopionLihMuon_03240.root",
                        1,"picon"],
                        ["../fiducial/LiHMu_3240_R_-1_G_1.root",
                        "../fiducial/LiHMu_3240_R_1_G_1.root",
                        0.478/20.,"Fid. Radius"],
                        ["../alignment/LiHMu_3240_-1.root",
                        "../alignment/LiHMu_3240_1.root",
                        1.2e-5/6.5e-5,"Alignment"]]

    # [Goldout172, Goldoutsys172] = Goldfile172.MCSPlot("raw")
    # [Goldout200, Goldoutsys200] = Goldfile200.MCSPlot("raw")
    # [Goldout240, Goldoutsys240] = Goldfile240.MCSPlot("raw")

    Goldfile172.useeventnorm = False
    Goldfile200.useeventnorm = False
    Goldfile240.useeventnorm = False
    Goldfile200.histvarnames = ['thetaX','thetaY']
    # Goldfile200.histstatenames = ['graph', 'recoGold_noeffi']
    # Goldfile200.histstatedesc = ['Truth', 'Deconvolved Gold']
    Goldfile200.histstatenames = ['recoGold_effi_only_gold', 'graph']
    Goldfile200.histstatedesc = ['Deconvolved Gold', 'Truth']
    # [Goldout200, Goldoutsys200] = Goldfile200.MCSPlot("Gold")



