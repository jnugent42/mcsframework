from math import sqrt
import libxml2
import sys, os, subprocess, math, ROOT
from ROOT import TH1D, TCanvas, TFile, TGraphErrors, TMarker, TF1, TLegend, TText, TLatex

class plotgen:
    
    def __init__(self, filename):
        self.fname = filename
        self.desc = self.fname.split("_")
        self.histvarnames = ['thetaX','thetaY','thetaScatt','theta2Scatt']
        self.histstatenames = ['measdataCobb','refconv_GEANT','refconv_Cobb', 'ref', 'GEANT', 'Moliere']
        self.histstatedesc = ['Raw Data',
                              'GEANT4 Default MCS', 'Carlisle-Cobb', 'Moliere']
        self.histcolors = [1, 4, 2, 6]
        self.histopts  = ['lp','lp','lp','lp']
        self.RMS = {}
        self.RMSErr = {}
        self.Chi2 = {}
        self.RMSsysdiff = {}
        self.RMSsyserr  = {}
        self.sysFiles   = []
        self.eventnorm  = {}
        self.h0integral = {}
        self.useeventnorm = True

        
    def formatHists(self, hist, i):
        # hist.Sumw2()
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


    def addToSysErr(self, i, histl, histh, base, scale, histvar, sysname, sysangdef):
        norm = 0
        if i == 0:
            norm = self.eventnorm[histvar]
        else:
            norm = histl.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
            # if histvar == 'thetaScatt' or histvar == 'theta2Scatt':
            #     norm /= base.GetXaxis().GetBinWidth(10)/histl.GetXaxis().GetBinWidth(10)
        # print norm
        histl.Scale(1./norm)
        norm = 0
        if i == 0:
            norm = self.eventnorm[histvar]
        else:
            norm = histh.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
            # if histvar == 'thetaScatt' or histvar == 'theta2Scatt':
            #     norm /= base.GetXaxis().GetBinWidth(10)/histh.GetXaxis().GetBinWidth(10)
        # print norm
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
            print sysname
	    print histvar
	    print histh
	    print histl
	    self.RMSsysdiff[histvar][self.histstatenames[i]][sysname] \
                = 1000*(histh.GetFunction('gaus').GetParameter('Sigma') - \
                        histl.GetFunction('gaus').GetParameter('Sigma'))
            self.RMSsyserr[histvar][self.histstatenames[i]][sysname] \
                = math.fabs(self.RMSsysdiff[histvar][self.histstatenames[i]][sysname]) * scale 
	'''
	if sysname == 'angdef':
   	       self.RMSsysdiff[histvar][self.histstatenames[i]][sysname] \
                   = 1000*(histangdef.GetRMS(2)
               self.RMSsyserr[histvar][self.histstatenames[i]][sysname] \
                   = math.fabs(self.RMSsysdiff[histvar][self.histstatenames[i]][sysname]) * scale 
        '''
            
    def addToRMS(self, i, hist, base, resplot, histvar):
        norm = 0
        if i == 0:
            if self.useeventnorm:
                self.eventnorm[histvar]  = hist.GetEntries()
            self.h0integral[histvar] = hist.Integral()
            norm = self.eventnorm[histvar]
            print histvar, i, norm, hist.Integral(), hist.Integral()/norm
        else:
            norm = hist.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
            # if histvar == 'thetaScatt' or histvar == 'theta2Scatt':
            #     norm /= base.GetXaxis().GetBinWidth(10)/hist.GetXaxis().GetBinWidth(10)
        # print norm
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
	    print self.RMS

    def calculateChi2(self, i, hist, base, resplot, histvar, pname):
        
        chi2 = 0
        ndf = 0
        # if i%2 == 0 and i > 0:
        sysLowF  = [TFile(x[0]) for x in self.sysFiles]
        sysLowH  = [f.Get(histvar + '_' + self.histstatenames[i]) for f in sysLowF]
        sysangdef = [f.Get(histvar + '_' + self.histstatenames[i] + '_' self.histangdef[0]) for f in self.fname]

        sysHighF = [TFile(x[1]) for x in self.sysFiles]
        sysHighH = [f.Get(histvar + '_' + self.histstatenames[i]) for f in sysHighF]
        sysScale = [x[2] for x in self.sysFiles]
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
        if histvar == 'theta2Scatt':
            binmax = 0.0036
        for k in range(len(self.sysFiles)):
            if histvar == 'theta2Scatt':
                sysLowH[k].Rebin(8)
                sysHighH[k].Rebin(8)
            elif histvar == 'thetaScatt':
                sysLowH[k].Rebin(1)
                sysHighH[k].Rebin(1)
            elif self.sysFiles[k][-1] != "TOF":
		print self.sysFiles[k] 
		print sysLowH[k] 
		sysLowH[k].Rebin(4)
                sysHighH[k].Rebin(4)
            self.formatHists(sysLowH[k], i)
            self.formatHists(sysHighH[k], i)
            self.addToSysErr(i, sysHighH[k], sysLowH[k], base, sysScale[k], histvar, self.sysFiles[k][3], sysangdef)
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
        for j in range(0, hist.GetNbinsX()+1):
            res = hist.GetBinContent(j) - base.GetBinContent(j)                
            #sys = hsysup[i].GetBinContent(j) - hsysdn[i].GetBinContent(j)
            #sys *= 0.06/0.8
            err2 = 0
            if math.fabs(res) > 0:
                err2 = hist.GetBinError(j)**2  # + base.GetBinError(j)**2
            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < 0.6 and err2 > 0:
                for k in range(len(self.sysFiles)):
                    sys = sysHighH[k].GetBinContent(j) - sysLowH[k].GetBinContent(j)
                    sys *= sysScale[k]
                    syshists[k].SetBinContent(j, sys*sys/err2)
                    syshists[k].SetBinError(j, 0.)
                    err2 += sys*sys
                if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax:
                    chi2 += res*res/err2
                resplot.SetBinContent(j, res/math.sqrt(err2))
                resplot.SetBinError(j, 1)
            else:
                resplot.SetBinContent(j, res/1)
                resplot.SetBinError(j, 0)
                
            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax:
                ndf += 1
                # print i, resplot.GetXaxis().GetBinCenter(j), res, chi2
        c = TCanvas(self.fname[:-5]+'_'+histvar+'_c1')
        c.SetBottomMargin(0.15)
        c.SetTopMargin(0.075)
        
        if pname == 'rawG4' or pname == 'decoG4':
            for k in range(len(self.sysFiles)):
                syshists[k].Draw("p")
                c.Print(pname+'_'+self.fname[:-5]+'_'+histvar+'_'+self.sysFiles[k][3]+'_res.eps')
        # print "Chi-square from residuals is ",chi2, " for " ,ndf," from ", self.histstatedesc[i]
        self.Chi2[histvar][self.histstatenames[i]] = [chi2, ndf]
        
    def MCSPlot(self, pname):

        print self.fname
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
            names = [histvar + '_' + x for x in self.histstatenames]
            # print names
            hists = [f.Get(histvar+'_'+x) for x in self.histstatenames]
            
            hists[0].SetTitle("")
            # hists[3].Scale(norm)
            # self.formatHist(hist[0], 0)
            resplots = [x.Clone() for x in hists]
            resplots[0].SetTitle('')
            resplots[0].GetYaxis().SetTitle("Normalized Residuals")
        
            # if histvar == 'thetaScatt':
            leg = TLegend(0.55,0.73,0.89,0.92)
            leg.SetLineColor(10)
            # else:
            #    leg = TLegend(0.35,0.2,0.65,0.5)
            for i in range(len(self.histstatedesc)):
                hists[i].Sumw2()
                if histvar=='theta2Scatt':
                    hists[i].Rebin(8)
                    resplots[i].Rebin(8)
                elif histvar=='thetaScatt':
                    hists[i].Rebin(1)
                    resplots[i].Rebin(1)
                else:
                    hists[i].Rebin(1)
                    resplots[i].Rebin(1)
                self.formatHists(hists[i], i)
                self.formatHists(resplots[i], i)
                self.addToRMS(i, hists[i], hists[0], resplots[i], histvar)
                if histvar=='theta2Scatt':
                    hists[i].GetYaxis().SetTitle('Probability per '+str(round(1000*1000*hists[i].GetXaxis().GetBinWidth(4),2))+' mrad^{2}')
                else:
                    hists[i].GetYaxis().SetTitle('Probability per '+str(round(1000*hists[i].GetXaxis().GetBinWidth(4),2))+' mrad')
                leg.AddEntry(hists[i], self.histstatedesc[i], self.histopts[i])
		print hists[0]
                self.calculateChi2(i, hists[i], hists[0], resplots[i], histvar, pname)

            
            c = TCanvas(self.fname[:-5]+'_'+histvar+'_c1')
            if self.desc[0] == 'XePion':
                t1 = TText(0.18,0.885,"MICE ISIS cycle 2015/03")
                t2 = TText(0.18,0.85,"Xe, "+self.desc[1][2:5]+", MAUS v2.9.1")
            else:
                t1 = TText(0.18,0.885,"MICE ISIS cycle 2015/04")
                t2 = TText(0.18,0.85,"LiH, "+self.desc[1][2:5]+", MAUS v2.9.1")
            t1.SetNDC(1)
            t1.SetTextSize(0.04)
            t1.SetTextFont(42)
            t2.SetNDC(1)
            t2.SetTextSize(0.03)
            t2.SetTextFont(42)
            
            hists[0].GetYaxis().SetRangeUser(4e-4,2.0)
        
            hists[0].SetTitle(";"+hists[0].GetXaxis().GetTitle()+" (radians);"+hists[0].GetYaxis().GetTitle())
            hists[0].Draw('ep')
            
            c.SetBottomMargin(0.15)
            c.SetTopMargin(0.075)

            for h in hists[1:len(self.histstatedesc)]:
                h.Draw('epsame')
            leg.SetTextSize(0.04)
            leg.Draw('same')
            t1.Draw()
            t2.Draw()
            c.SetLogy()
            c.Print(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys.eps')
            c.Print(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_pq.jpg')
            c.Clear()
            
            c.SetLogy(0)
            resplots[0].GetYaxis().SetRangeUser(-10,8)
            resplots[0].SetTitle(";"+resplots[0].GetXaxis().GetTitle()+" (radians);"+resplots[0].GetYaxis().GetTitle())
            leg.SetX1NDC(0.5)
            leg.SetX2NDC(0.89)
            leg.SetY1NDC(0.2)
            leg.SetY2NDC(0.4)
            resplots[0].Draw("p")
            for r in resplots:
                r.Draw('psame')
            leg.SetTextSize(0.04)
            leg.Draw('same')
            t1.Draw()
            t2.Draw()
            # pblock.Draw()
            c.Print(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_res_T.eps')
            c.Print(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_res_pq.jpg')

        momhist = f.Get("thetaScatt_measdata_vpGEANT")    
        #mom = [momhist.GetMean() + 19.468, momhist.GetMeanError()]
        #if self.fname.find("LiHMuon_03172")  >= 0:
        #    mom = [momhist.GetMean()*1.107 + 1.05, momhist.GetMeanError()]
        #elif self.fname.find("LiHMuon_03200")  >= 0:
        #    mom = [momhist.GetMean()*1.104 + 1.139, momhist.GetMeanError()]
        #elif self.fname.find("LiHMuon_03240")  >= 0:
        #    mom = [momhist.GetMean()*1.17 - 9.41, momhist.GetMeanError()]
        if self.fname.find("LiHMuon_03172")  >= 0:
            mom = [momhist.GetMean(), momhist.GetMeanError()]
        elif self.fname.find("LiHMuon_03200")  >= 0:
            mom = [momhist.GetMean(), momhist.GetMeanError()]
        elif self.fname.find("LiHMuon_03240")  >= 0:
            mom = [momhist.GetMean(), momhist.GetMeanError()]
        rms = [momhist.GetRMS(), momhist.GetRMSError()]
        summary = []
        syssummary = []
        def sigfig(x):
            if math.fabs(x) > 1e-5:
                return int(math.ceil(math.fabs(math.log(math.fabs(x),10))))
            else:
                return 1
        
        # syssummary.append("p (MeV/c) & "+self.histvarnames[0]+"&"+self.histvarnames[1]+"&"+self.histvarnames[3]+"\\\\")

        
        for sys in self.sysFiles:
            # if sys[3] == 'Material': stindx = 1
            # else:
            stindx = 0
            # print sys[3], self.histstatenames[stindx]
            difference0 = self.RMSsysdiff[self.histvarnames[0]][self.histstatenames[stindx]][sys[3]]
            difference1 = self.RMSsysdiff[self.histvarnames[1]][self.histstatenames[stindx]][sys[3]]
            difference3 = self.RMSsysdiff[self.histvarnames[3]][self.histstatenames[stindx]][sys[3]]

            syserr0 = self.RMSsyserr[self.histvarnames[0]][self.histstatenames[stindx]][sys[3]]
            syserr1 = self.RMSsyserr[self.histvarnames[1]][self.histstatenames[stindx]][sys[3]]
            syserr3 = self.RMSsyserr[self.histvarnames[3]][self.histstatenames[stindx]][sys[3]]

            relerr0 = syserr0/self.RMS[self.histvarnames[0]][self.histstatenames[0]]
            relerr1 = syserr1/self.RMS[self.histvarnames[1]][self.histstatenames[0]]
            relerr3 = syserr3/self.RMS[self.histvarnames[3]][self.histstatenames[0]]
            
            syssummary.append(str(round(mom[0],sigfig(mom[1])))+"$\pm$"+str(round(mom[1],sigfig(mom[1])))+\
                              " & "+str(round(difference0,sigfig(difference0)))+\
                              " & "+str(round(syserr0,sigfig(syserr0)))+\
                              " & "+str(round(relerr0,sigfig(relerr0)))+\
                              " & "+str(round(difference1,sigfig(difference1)))+\
                              " & "+str(round(syserr1,sigfig(syserr1)))+\
                              " & "+str(round(relerr1,sigfig(relerr1)))+\
                              " & "+str(round(difference3,sigfig(difference3)))+\
                              " & "+str(round(syserr3,sigfig(syserr3)))+\
                              " & "+str(round(relerr3,sigfig(relerr3)))+"\\\\")
            # print syssummary[-1]
            
        syssummary.append(str(round(mom[0],2))+"$\pm$"+str(round(mom[1],2))+\
                          " & "+ str(round(rms[0],2))+"$\pm$"+str(round(rms[1],2))+\
                          " & "+str(round(self.RMSsysdiff[self.histvarnames[0]][self.histstatenames[0]]['Sum'],2))+\
                          " & "+str(round(self.RMSsyserr[self.histvarnames[0]][self.histstatenames[0]]['Sum'],2))+\
                          " & "+str(round(self.RMSsyserr[self.histvarnames[0]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[0]][self.histstatenames[0]],2))+\
                          " & "+str(round(self.RMSsysdiff[self.histvarnames[1]][self.histstatenames[0]]['Sum'],2))+\
                          " & "+str(round(self.RMSsyserr[self.histvarnames[1]][self.histstatenames[0]]['Sum'],2))+\
                          " & "+str(round(self.RMSsyserr[self.histvarnames[1]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[1]][self.histstatenames[0]],2))+\
                          " & "+str(round(self.RMSsysdiff[self.histvarnames[3]][self.histstatenames[0]]['Sum'],2))+
                          " & "+str(round(self.RMSsyserr[self.histvarnames[3]][self.histstatenames[0]]['Sum'],2))+
                          " & "+str(round(self.RMSsyserr[self.histvarnames[3]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[2]][self.histstatenames[0]],2))+"\\\\")
            
            
        # summary.append("p (MeV/c) &  &"+str(self.histstatenames[0])+" & "+str(self.histstatenames[1])+" & $\chi^{2}$/ndf & "+
        #                +str(self.histstatenames[2])+" & $\chi^{2}$/ndf \\\\")
        # print mom, self.RMS, self.RMSErr, self.Chi2
        for histvar in self.histvarnames:
            summary.append(str(round(mom[0],2))+"$\pm$"+str(round(mom[1],2))+\
                           "& "+histvar+" & "+str(round(self.RMS[histvar][self.histstatenames[0]],2))+ \
                           "$\pm$"+str(round(self.RMSErr[histvar][self.histstatenames[0]],2))+ \
                           "$\pm$"+str(round(self.RMSsyserr[histvar][self.histstatenames[0]]["Sum"],2))+ \
                           " & "+str(round(self.RMS[histvar][self.histstatenames[1]],2))+ \
                           "$\pm$"+str(round(self.RMSErr[histvar][self.histstatenames[1]],2))+ \
                           # "$\pm$"+str(round(self.RMSsyserr[histvar][self.histstatenames[1]]["Sum"],2))+\
                           " & "+str(round(self.Chi2[histvar][self.histstatenames[1]][0],1))+ \
                           " / "+ str(self.Chi2[histvar][self.histstatenames[1]][1])+ \
                           " & "+str(round(self.RMS[histvar][self.histstatenames[2]],2))+ \
                           "$\pm$"+str(round(self.RMSErr[histvar][self.histstatenames[2]],2))+ \
                           #"$\pm$"+str(round(self.RMSsyserr[histvar][self.histstatenames[2]]["Sum"],2))+\
                           " & "+str(round(self.Chi2[histvar][self.histstatenames[2]][0],1))+ \
                           " / "+ str(self.Chi2[histvar][self.histstatenames[2]][1])
                           +"\\\\")
            # print summary[-1]
            
        f.Close()
        return [summary, syssummary]
    #  print "+++++++++++++++ systematics +++++++++++++++++++++++"
    # for i in range(len(RMS)):
    #     print " & ", histvarnames[i]," & ",round(RMSsysdiff[i][0],2)," & ",round(RMSsyserr[i][0],3)," & ",round(RMSsysdiff[i][1],2)," & ",round(RMSsyserr[i][1],3)," & ",round(RMSsysdiff[i][2],2)," & ",round(RMSsyserr[i][2],3),"\\\\"
        
if __name__ == '__main__':
    # baseSet = sys.argv[1]
    raw172 = plotgen("../lih172/LiHMuon_03172.root")
    raw172.sysFiles = [["../TOFsys/LiHMu_3172_tof_lim28.7374090796_u28.9374090796.root",
                        "../TOFsys/LiHMu_3172_tof_lim29.5374090796_u29.7374090796.root",
                        129./800.,"TOF"],
                       ["../materials/material_0/LiHMu_3172_res_lim0.891.root",
                        "../materials/material_0/LiHMu_3172_res_lim1.109.root",
                        sqrt(0.004/0.218),"Material"],
                       ["../angdef/LiHMu_3172_45.root",
                        "../angdef/LiHMu_3172_135.root",
                        1,"angdef"],
                       ["../MClih172/LiHMu_3172.root",
                        "../MClih172/acceptLiHMu_3172.root",
                        1,"accept"],
                       ["../MClih172/LiHMu_3172.root",
                        "../MClih172/nopionLiHMu_3172.root",
                        1,"picon"],
                       ["../fiducial/LiHMu_3172_R_5_G_-17.root",
                        "../fiducial/LiHMu_3172_R_5_G_-15.root",
                        0.5/10,"Fid. Pitch"],
                       ["../fiducial/LiHMu_3172_R_4_G_-17.root",
                        "../fiducial/LiHMu_3172_R_6_G_-17.root",
                       .478/20.,"Fid. Radius"],
                       ["../alignment/LiHMu_3172_13.root",
                        "../alignment/LiHMu_3172_10.root",
                        8.9e-5/8.1e-4,"Alignment"]]

    raw200 = plotgen("../lih200/LiHMuon_03200.root")
    raw200.sysFiles = [["../TOFsys/LiHMu_3200_tof_lim27.8059285384_u28.0059285384.root",
                        "../TOFsys/LiHMu_3200_tof_lim28.6059285384_u28.8059285384.root",
                        129./800.,"TOF"],
                       ["../materials/material_0/LiHMu_3200_res_lim0.891.root",
                        "../materials/material_0/LiHMu_3200_res_lim1.109.root",
                        sqrt(0.004/0.218),"Material"],
                       ["../angdef/LiHMu_3200_45.root",
                        "../angdef/LiHMu_3200_135.root",
                        1,"angdef"],
                       ["../MClih200/LiHMu_3200.root",
                        "../MClih200/nopionLiHMu_3200.root",
                        1,"picon"],
                       ["../MClih200/LiHMu_3200.root",
                        "../MClih200/acceptLiHMu_3200.root",
                        1,"accept"],
                       ["../fiducial/LiHMu_3200_R_5_G_-17.root",
                        "../fiducial/LiHMu_3200_R_5_G_-15.root",
                        0.5/10,"Fid. Pitch"],
                       ["../fiducial/LiHMu_3200_R_4_G_-17.root",
                        "../fiducial/LiHMu_3200_R_6_G_-17.root",
                        0.478/20.,"Fid. Radius"],
                       ["../alignment/LiHMu_3200_1.root",
                        "../alignment/LiHMu_3200_7.root",
                        2.1e-5/9.1e-5,"Alignment"]]
    
    raw240 = plotgen("../lih240/LiHMuon_03240.root")
    raw240.sysFiles = [["../TOFsys/LiHMu_3240_tof_lim27.0570317938_u27.2570317938.root",
                        "../TOFsys/LiHMu_3240_tof_lim27.8570317938_u28.0570317938.root",
                        129./800.,"TOF"],
                       ["../materials/material_0/LiHMu_3240_res_lim0.891.root",
                        "../materials/material_0/LiHMu_3240_res_lim1.109.root",
                        sqrt(0.004/0.218),"Material"],
                       ["../angdef/LiHMu_3240_45.root",
                        "../angdef/LiHMu_3240_135.root",
                        1,"angdef"],
                       ["../MClih240/LiHMu_3240.root",
                        "../MClih240/nopionLiHMu_3240.root",
                        1,"picon"],
                       ["../MClih240/LiHMu_3240.root",
                        "../MClih240/acceptLiHMu_3240.root",
                        1,"accept"],
                       ["../fiducial/LiHMu_3240_R_5_G_-17.root",
                        "../fiducial/LiHMu_3240_R_5_G_-15.root",
                        0.5/10,"Fid. Pitch"],
                       ["../fiducial/LiHMu_3240_R_4_G_-17.root",
                        "../fiducial/LiHMu_3240_R_6_G_-17.root",
                        0.478/20.,"Fid. Radius"],
                       ["../alignment/LiHMu_3240_59.root",
                        "../alignment/LiHMu_3240_19.root",
                        1.2e-5/6.5e-5,"Alignment"]]
    
    [out172, sys172] = raw172.MCSPlot("raw")
    [out200, sys200] = raw200.MCSPlot("raw")
    [out240, sys240] = raw240.MCSPlot("raw")
    '''
    piraw240 = plotgen("XePion_3240_0.root")
    piraw240.sysFiles = [["TOFsys/XePion_3240_tof_lim26.6_u28.0.root",
                          "TOFsys/XePion_3240_tof_lim27.4_u28.8.root",
                          128./800.,"TOF"],
                         ["material_0/XePion_3240_res_lim0.891.root",
                          "material_0/XePion_3240_res_lim1.109.root",
                          sqrt(0.010/0.218),"Material"],
                       ["../angdef/LiHMu_3172_45.root",
                        "../angdef/LiHMu_3172_135.root",
                        1,"angdef"],
                         ["fiducial_0/XePion_3240_R_5_G_0.root",
                          "fiducial_0/XePion_3240_R_5_G_2.root",
                          0.5/10,"Fid. Pitch"],
                         ["fiducial_0/XePion_3240_R_4_G_1.root",
                          "fiducial_0/XePion_3240_R_6_G_1.root",
                          0.478/20.,"Fid. Radius"]]

    piraw240.histstatenames = ['measdataCobb','refconv_truthsim','refconv_Cobb', 'ref', 'truthsim']
    [piout240, pisys240] = piraw240.MCSPlot("raw")
    '''
    
    print "\n"
    raw172.histstatenames = ['measdataGEANT','recoGEANT','GEANT']
    raw172.histstatedesc = ['Raw Data','Deconvolved, GEANT', 'GEANT Model']
    raw200.histstatenames = ['measdataGEANT','recoGEANT','GEANT']
    raw200.histstatedesc = ['Raw Data','Deconvolved, GEANT', 'GEANT Model']
    raw240.histstatenames = ['measdataGEANT','recoGEANT','GEANT']
    raw240.histstatedesc = ['Raw Data','Deconvolved, GEANT', 'GEANT Model']
    raw172.useeventnorm = True
    raw200.useeventnorm = True
    raw240.useeventnorm = True
    raw172.MCSPlot("rawG4")
    raw200.MCSPlot("rawG4")
    raw240.MCSPlot("rawG4")
    
    print "\n"
    raw172.histstatenames = ['measdataCobb','recoCobb','Cobb']
    raw172.histstatedesc = ['Raw Data','Deconv., Carlisle-Cobb', 'Carlisle-Cobb Model']
    raw200.histstatenames = ['measdataCobb','recoCobb','Cobb']
    raw200.histstatedesc = ['Raw Data','Deconv., Carlisle-Cobb', 'Carlisle-Cobb Model']
    raw240.histstatenames = ['measdataCobb','recoCobb','Cobb']
    raw240.histstatedesc = ['Raw Data','Deconv., Carlisle-Cobb', 'Carlisle-Cobb Model']
    raw172.MCSPlot("rawCC")
    raw200.MCSPlot("rawCC")
    raw240.MCSPlot("rawCC")

    print "\n"
    raw172.histstatenames = ['measdataGEANT','ref','GEANT']
    raw172.histstatedesc = ['Raw Data','Empty AFC', 'GEANT Model']
    raw200.histstatenames = ['measdataGEANT','ref','GEANT']
    raw200.histstatedesc = ['Raw Data','Empty AFC', 'GEANT Model']
    raw240.histstatenames = ['measdataGEANT','ref','GEANT']
    raw240.histstatedesc = ['Raw Data','Empty AFC', 'GEANT Model']
    # piraw240.histstatenames = ['measdatatruthsim','ref','truthsim']
    # piraw240.histstatedesc = ['Raw Data','Empty AFC', 'GEANT Model']
    [ref172, refsys172] = raw172.MCSPlot("refG4")
    [ref200, refsys200] = raw200.MCSPlot("refG4")
    [ref240, refsys240] = raw240.MCSPlot("refG4")
    # [piref240, pirefsys240] = piraw240.MCSPlot("refG4")
    # piraw240.histstatenames = ['measdatatruthsim','ref','Cobb']
    # piraw240.histstatedesc = ['Raw Data','Empty AFC', 'Carlisle-Cobb Model']
    # piraw240.MCSPlot("refCC")
    
    print "\n"
    raw172.histstatenames = ['recoGEANT','GEANT','Cobb']
    raw172.histstatedesc = ['Deconvolved, GEANT', 'GEANT Model', 'Carlisle-Cobb Model']
    raw200.histstatenames = ['recoGEANT','GEANT','Cobb']
    raw200.histstatedesc = ['Deconvolved, GEANT', 'GEANT Model', 'Carlisle-Cobb Model']
    raw240.histstatenames = ['recoGEANT','GEANT','Cobb']
    raw240.histstatedesc = ['Deconvolved, GEANT', 'GEANT Model', 'Carlisle-Cobb Model']
    # piraw240.histstatenames = ['recotruthsim','truthsim','Cobb']
    # piraw240.histstatedesc = ['Deconvolved, GEANT', 'GEANT Model', 'Carlisle-Cobb Model']
    raw172.useeventnorm = False
    raw200.useeventnorm = False
    raw240.useeventnorm = False
    # piraw240.useeventnorm = False
    [G4out172, G4sys172] = raw172.MCSPlot("decoG4")
    [G4out200, G4sys200] = raw200.MCSPlot("decoG4")
    [G4out240, G4sys240] = raw240.MCSPlot("decoG4")
    # [piG4out240, piG4sys240] = piraw240.MCSPlot("decoG4")
    raw172.histstatenames = ['recoMoliere','Moliere','Cobb']
    raw172.histstatedesc = ['Deconvolved, Moliere', 'Moliere Model', 'Carlisle-Cobb Model']
    raw200.histstatenames = ['recoMoliere','Moliere','Cobb']
    raw200.histstatedesc = ['Deconvolved, Moliere', 'Moliere Model', 'Carlisle-Cobb Model']
    raw240.histstatenames = ['recoMoliere','Moliere','Cobb']
    raw240.histstatedesc = ['Deconvolved, Moliere', 'Moliere Model', 'Carlisle-Cobb Model']
    # piraw240.histstatenames = ['recotruthsim','truthsim','Cobb']
    # piraw240.histstatedesc = ['Deconvolved, Moliere', 'Moliere Model', 'Carlisle-Cobb Model']
    raw172.useeventnorm = False
    raw200.useeventnorm = False
    raw240.useeventnorm = False
    # piraw240.useeventnorm = False
    [Moliereout172, Molieresys172] = raw172.MCSPlot("decoG4")
    [Moliereout200, Molieresys200] = raw200.MCSPlot("decoG4")
    [Moliereout240, Molieresys240] = raw240.MCSPlot("decoG4")
    print "\n"
    raw172.histstatenames = ['recoCobb','GEANT','Cobb']
    raw172.histstatedesc = ['Deconv., Carlisle-Cobb','GEANT Model', 'Carlisle-Cobb Model']
    raw200.histstatenames = ['recoCobb','GEANT','Cobb']
    raw200.histstatedesc = ['Deconv., Carlisle-Cobb','GEANT Model', 'Carlisle-Cobb Model']
    raw240.histstatenames = ['recoCobb','GEANT','Cobb']
    raw240.histstatedesc = ['Deconv., Carlisle-Cobb','GEANT Model', 'Carlisle-Cobb Model']
    # piraw240.histstatenames = ['recoCobb','truthsim','Cobb']
    # piraw240.histstatedesc = ['Deconv., Carlisle-Cobb','GEANT Model', 'Carlisle-Cobb Model']
    [CCout172, CCsys172] = raw172.MCSPlot("decoCC")
    [CCout200, CCsys200] = raw200.MCSPlot("decoCC")
    [CCout240, CCsys240] = raw240.MCSPlot("decoCC")
    plotFiles(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
    # [piCCout240, piCCsys240] = piraw240.MCSPlot("decoCC")

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
    # print '\hline'
    # print piout240[0]
    # print piout240[1]

    print '\hline'
    print '\hline'
    print out172[3]
    print out200[3]
    print out240[3]
    # print piout240[3]
    print '\hline'
    print '\hline\n'

    print '\hline'
    print '\hline'
    print ref172[0]
    print ref172[1]
    print '\hline'
    print ref200[0]
    print ref200[1]
    print '\hline'
    print ref240[0]
    print ref240[1]
    # print '\hline'
    # print piref240[0]
    # print piref240[1]

    print '\hline'
    print '\hline'
    print ref172[3]
    print ref200[3]
    print ref240[3]
    # print piref240[3]
    print '\hline'
    print '\hline\n'

    
    print "GEANT4 deco table" 
    print '\hline'
    print '\hline'
    print G4out172[0]
    print G4out172[1]
    print '\hline'
    print G4out200[0]
    print G4out200[1]
    print '\hline'
    print G4out240[0]
    print G4out240[1]
    # print '\hline'
    # print piG4out240[0]
    # print piG4out240[1]

    print '\hline'
    print '\hline'
    print G4out172[3]
    print G4out200[3]
    print G4out240[3]
    # print piG4out240[3]
    print '\hline'
    print '\hline\n'

    print "Moliere deco table" 
    print '\hline'
    print '\hline'
    print Moliereout172[0]
    print Moliereout172[1]
    print '\hline'
    print Moliereout200[0]
    print Moliereout200[1]
    print '\hline'
    print Moliereout240[0]
    print Moliereout240[1]
    # print '\hline'
    # print piMoliereout240[0]
    # print piMoliereout240[1]

    print '\hline'
    print '\hline'
    print Moliereout172[3]
    print Moliereout200[3]
    print Moliereout240[3]
    # print piG4out240[3]
    print '\hline'
    print '\hline\n'
    
    print "Cobb-Carlisle deco table" 
    print '\hline'
    print '\hline'
    print CCout172[0]
    print CCout172[1]
    print '\hline'
    print CCout200[0]
    print CCout200[1]
    print '\hline'
    print CCout240[0]
    print CCout240[1]
    # print '\hline'
    # print piCCout240[0]
    # print piCCout240[1]

    print '\hline'
    print '\hline'
    print CCout172[3]
    print CCout200[3]
    print CCout240[3]
    # print piCCout240[3]
    print '\hline'
    print '\hline'

    
    print '\hline'
    print '\hline'
   
    print 'raw sys'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3]
        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
        print sys172[i]
        print sys200[i]
        print sys240[i]
        # print pisys240[i]
        print '\hline\n'
    print '\hline'
    
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print sys172[-1]
    print sys200[-1]
    print sys240[-1]
    # print pisys240[-1]
    print '\hline\n'
    print '\hline'
    print 'G4 sys'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3]
        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
        print G4sys172[i]
        print G4sys200[i]
        print G4sys240[i]
        # print piG4sys240[i]
        print '\hline\n'
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print sys172[-1]
    print sys200[-1]
    print sys240[-1]
    # print pisys240[-1]
    print '\hline\n'
    print '\hline'
    print 'Moliere sys'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3]
        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
        #print Molieresys172[i]
        #print Molieresys200[i]
        #print Molieresys240[i]
        # print piG4sys240[i]
        print '\hline\n'
    print '\hline'
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print G4sys172[-1]
    print G4sys200[-1]
    print G4sys240[-1]
    # print piG4sys240[-1]
    print '\hline\n'
    print '\hline'
    print 'CC sys'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3]
        print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
        print CCsys172[i]
        print CCsys200[i]
        print CCsys240[i]
        print '\hline\n'
    #print raw172.RMSsyserr
    #print raw200.RMSsyserr
    #print raw240.RMSsyserr

    raw200.histvarnames = ['thetaX']
    raw200.histstatenames = ['theta_true_x_graph']
    raw200.histstatedesc = ['Truth Data']
    raw200.MCSPlot('')

'''
def plotFiles(meas172, measfile, meas240, syslist):
    
    hists  = []
    graphs = []

    m172 = ExtractPars(meas172)
    mval = ExtractPars(measfile)
    m240 = ExtractPars(meas240)
    m172G4 = ExtractParsG4(meas172)
    m200G4 = ExtractParsG4(measfile)
    m240G4 = ExtractParsG4(meas240)
    print "mval",mval
    print "m172",m172
    print "m240",m240
    pars = [["TOF_ll",-0.8,0.8], ["TOF_ul",0.5,2.1], 
            ["integralX",0,2*mval[2][0]], ["amplitudeX", 0, 2*mval[3][0]], 
            ["rawmeanX",mval[4][0]-3*mval[4][1], mval[4][0]+3*mval[4][1]],
            ["rawsigmaX",mval[5][0]-3*mval[5][1], mval[5][0]+3*mval[5][1]],
            ["meanX",mval[6][0]-3*mval[6][1], mval[6][0]+3*mval[6][1]], 
            ["sigmaX",mval[7][0]-3*mval[7][1], mval[7][0]+3*mval[7][1]],
            ["integralY",0,2*mval[8][0]], ["amplitudeY", 0, 2*mval[9][0]], 
            ["rawmeanY",mval[10][0]-3*mval[10][1], mval[10][0]+3*mval[10][1]],
            ["rawsigmaY",mval[11][0]-3*mval[11][1], mval[11][0]+3*mval[11][1]],
            ["meanY",mval[12][0]-3*mval[12][1], mval[12][0]+3*mval[12][1]], 
            ["sigmaY",mval[13][0]-3*mval[13][1], mval[13][0]+3*mval[13][1]],
            ["rawmeantheta2",mval[14][0]-3*mval[12][1], mval[14][0]+3*mval[12][1]], 
            ["meantheta2",mval[15][0]-3*mval[12][1], mval[15][0]+3*mval[12][1]], 
            ["momentum",mval[16][0]-3*mval[14][1], mval[16][0]+3*mval[14][1]]
            ]
    # for par in pars:
    #     hists.append(TH1D(par[0], ";"+par[0],100,par[1],par[2]))

    pointlist = []
    
    scale =  1. # 1.18 # 1.13 # 0.894
    offset =  0. # -11.28 # -5.95 # 46.4
    for sysfile in syslist:
        sysval = ExtractPars(sysfile)
        pcomp = 200
        
        if sysfile.find("LiHMu_3172") >= 0:
            pcomp = 174
            '''
p0                        =      1.05205   +/-   0.273774    
p1                        =      1.10656   +/-   0.0016539
            '''
            #sysval[16][0] *= 1.107
            #sysval[16][0] += 1.05
            #sysval[16][1] *= 1.107
        elif sysfile.find("LiHMu_3200") >= 0:
            pcomp = 200
            '''
p0                        =      1.13876   +/-   0.671571    
p1                        =      1.10475   +/-   0.00407707
            '''
            #sysval[16][0] *= 1.104
            #sysval[16][1] *= 1.104
            #sysval[16][0] += 1.139
        elif sysfile.find("LiHMu_3240") >= 0:
            pcomp = 240
            '''
p0                        =     -9.41424   +/-   0.293949    
p1                        =      1.17466   +/-   0.00142602  
            '''
            #sysval[16][0] *= 1.17
            #sysval[16][1] *= 1.17
            #sysval[16][0] -= 9.41
            # sysval[16][1] *= 0.780
        
        if sysval[2][0] > 5000.: # and pcomp < sysval[-1][0] + 10:
           print "sysval",sysval
           pointlist.append(sysval)
            
    for par in pars:
        graphs.append(TGraphErrors(len(pointlist)))
        graphs[-1].SetName(par[0])
        # graphs[-1].SetTitle("; #Delta t_{12}(ns);"
        if par[0] == "rawmeantheta2" or par[0] == "meantheta2":
            graphs[-1].SetTitle("; Momentum (MeV/c);  #sqrt{<#theta^{2}>/2} (milliradians)")
        elif par[0] == "rawsigmaX" or par[0] == "sigmaX":
            graphs[-1].SetTitle("; Momentum (MeV/c);  #Theta_{X} (milliradians)")
        elif par[0] == "rawsigmaY" or par[0] == "sigmaY":
            graphs[-1].SetTitle("; Momentum (MeV/c);  #Theta_{Y} (milliradians)")
        elif par[0] == "integralX":
            graphs[-1].SetTitle("; Momentum (MeV/c);  Integral of Events in Bin")
        else:
            graphs[-1].SetTitle("; Momentum (MeV/c);#Delta t_{01} (ns)")
    j = 0
    FillGraphs(graphs, pointlist, offset, scale)
    # func = TF1("func",'[0] + [1]*13.6*sqrt(x*x + 105.65 * 105.65)/x/x + [2]*13.6*13.6*(x*x + 105.65 * 105.65)/x/x/x/x',150,280)
    # func = TF1("func",'[0] + [1]*13.6*sqrt(1 + 105.65 * 105.65/x/x)/x',150,280)
    # func = TF1("func",'[0] + [1]*13.6*sqrt(1 + 105.65 * 105.65/(1.21*x - 14.9)/(1.21*x - 14.9))/(1.21*x - 14.9)',150,280)
    func = TF1("func",'[0]*13.6*sqrt(1 + 105.65 * 105.65/x/x)/(x)',140,280)
    # func = TF1("func",'[0]*13.6*105.65*105.65/sqrt(x*x + 105.65 * 105.65)/x/x',150,280)
    func0 = TF1("func0",'[0]*13.6*sqrt(1 + 105.65 * 105.65/x/x)/(x)',140,280)
    #ftof = TF1("ftof",'[0]/x + [1]', 26.8, 29.4)
    etof = (12929.2608098 - 5287.24720607) / 0.299792458 / 1000.
    ftof = TF1("ftof",'105*105/(y*y/etof*etof-1)', 26.8, 29.4)
    i = 0
    for graph in graphs:
        c = TCanvas()
        c.SetLeftMargin(0.125)
        c.SetBottomMargin(0.125)
        graph.SetMarkerStyle(21)
        graph.Draw('ap')
        graph.GetHistogram().GetXaxis().SetLabelSize(0.05)
        graph.GetHistogram().GetYaxis().SetLabelSize(0.05)
        graph.GetHistogram().GetXaxis().SetTitleSize(0.05)
        graph.GetHistogram().GetYaxis().SetTitleSize(0.05)
        graph.GetHistogram().GetYaxis().SetTitleOffset(1.25)
        if graph.GetName() == "integralX":
            leg = TLegend(0.5,0.15,0.89,0.39)
        else: 
            leg = TLegend(0.5,0.5,0.89,0.89)
        leg.SetLineColor(10)
        leg.AddEntry(graph,"Data","p")
        
        if graph.GetName() == "rawmeantheta2" or graph.GetName() == "meantheta2" \
           or graph.GetName() == "rawsigmaX" or graph.GetName() == "sigmaX" \
           or graph.GetName() == "rawsigmaY" or graph.GetName() == "sigmaY":
            func.SetLineWidth(1)
	    graph.Fit('func')
            x = ROOT.Double()
            y = ROOT.Double()
            for k in range(len(pointlist)):
                graph.GetPoint(k, x, y)
                ex = graph.GetErrorX(k)
                ey = graph.GetErrorY(k)
                fup = graph.GetFunction('func').Eval(x + 4)
                fdn = graph.GetFunction('func').Eval(x - 4)
                time = (float(pointlist[k][1][0]) + float(pointlist[k][0][0]))/2.
                dpdt = (1.04e3/(time - 23.2)**2)
                dsdp = (fup - fdn)/8.
                errsys = dpdt * 0.129  
		print k, j, pointlist[k][1][0], pointlist[k][0][0], x, y, fup, fdn, errsys, graph.GetErrorY(k), math.sqrt(pow(float(pointlist[k][i][1]),2)), \
                    math.sqrt(pow(float(pointlist[k][i][1]),2) + errsys*errsys)
                graph.SetPointError(k, math.sqrt(ex*ex + errsys*errsys), math.sqrt(ey*ey + 0.0*errsys*errsys))
            graph.Fit('func')
            
            leg.AddEntry(func, "Fit to #frac{13.6 a}{p#beta} (a="+str(round(func.GetParameter(0),2))+"#pm"+str(round(func.GetParError(0),2))+")", "l")
            fup = graph.GetFunction('func').Clone()
            fup.SetName('fup')
            fup.SetParameter(0, graph.GetFunction('func').GetParameter('p0') + \
                             graph.GetFunction('func').GetParError(0))
            # fup.SetParameter(1, graph.GetFunction('func').GetParameter('p1') + \
            #                 graph.GetFunction('func').GetParError(1))
            fup.SetLineColor(4)
            fup.SetLineWidth(1)
            leg.AddEntry(fup, "Fit plus/minus error", "l")
            fdn = graph.GetFunction('func').Clone()
            fdn.SetName('fdn')
            fdn.SetParameter(0, graph.GetFunction('func').GetParameter('p0') - \
                             graph.GetFunction('func').GetParError(0))
            # fdn.SetParameter(1, graph.GetFunction('func').GetParameter('p1') - \
            #                  graph.GetFunction('func').GetParError(1))
            fdn.SetLineColor(4)
            fdn.SetLineWidth(1)
            # leg.AddEntry(fdn, "Fit minus stat. error", "l")
            fup.Draw('lsame')
            fdn.Draw('lsame')
            func0.SetParameter(0, 253*(1 + 0.036*math.log(0.253*0.253)))
            func0.SetLineColor(8)
            #func0.Draw("lsame")
            #leg.AddEntry(func0,"PDG prediction","l")
            
        if graph.GetName() == "momentum":
            graph.Fit('ftof')
            print "TOF for 172 MeV/c is ",graph.GetFunction('ftof').Eval(172)
            print "TOF for 200 MeV/c is ",graph.GetFunction('ftof').Eval(200)
            print "TOF for 240 MeV/c is ",graph.GetFunction('ftof').Eval(240)
            leg.AddEntry(ftof, "Fit to #frac{m}{p} + c", "l")
            

        if j!=2 and j!=3 and j!=4 and j!=6 and j!=8 and j!=9 \
           and j!=10 and j!=12 and j!=len(graphs)-1:
            graph.GetHistogram().GetYaxis().SetRangeUser(10.1,28)
        elif j==4 or j==6 or j==10 or j==12:
            graph.GetHistogram().GetYaxis().SetRangeUser(-5,10)
            # graph.GetHistogram().GetYaxis().SetRangeUser(float(mval[j][0]) - 30*float(mval[j][1]),
            #                                             float(mval[j][0]) + 30*float(mval[j][1]))
        # m = TMarker((float(mval[1][0]) + float(mval[0][0]), float(mval[j][0]), 25)
        elif j==2:
            graph.GetHistogram().GetYaxis().SetRangeUser(0,10000)
        #graph.GetHistogram().GetXaxis().SetRangeUser(-125,300)
        #graph.GetHistogram().GetYaxis().SetRangeUser(0.1,50)
            
        if graph.GetName() == "momentum":
            #m = TMarker(1.104 * float(mval[j][0]) + 1.139, (float(mval[1][0]) + float(mval[0][0]))/2., 20)
            m = TMarker(float(mval[j][0]), (float(mval[1][0]) + float(mval[0][0]))/2., 20)
            m_172 = TMarker(1.107 * float(m172[j][0]) + 1.105, (float(m172[1][0]) + float(m172[0][0]))/2., 20)
            m_240 = TMarker(1.175 * float(m240[j][0]) - 9.41, (float(m240[1][0]) + float(m240[0][0]))/2., 20)
            # mG4_172 = TMarker(scale * float(m172G4[j][0]) + offset, (float(m172G4[1][0]) + float(m172G4[0][0]))/2., 20)
            # mG4_200 = TMarker(scale * float(m200G4[j][0]) + offset, (float(m200G4[1][0]) + float(m200G4[0][0]))/2., 20)
            # mG4_240 = TMarker(scale * float(m240G4[j][0]) + offset, (float(m240G4[1][0]) + float(m240G4[0][0]))/2., 20)
        else:
            #m = TMarker((1.104 * float(mval[-1][0]) + 1.139), float(mval[j][0]), 20)
            m = TMarker((float(mval[-1][0])), float(mval[j][0]), 20)
            m_172 = TMarker((float(m172[-1][0])), float(m172[j][0]), 20)
            m_240 = TMarker((float(m240[-1][0])), float(m240[j][0]), 20)
            #m_172 = TMarker((1.107 * float(m172[-1][0]) + 1.105), float(m172[j][0]), 20)
            #m_240 = TMarker((1.175 * float(m240[-1][0]) - 9.41), float(m240[j][0]), 20)
            
            mG4_172 = TMarker((1.107 * float(m172G4[-1][0]) + 1.105), float(m172G4[j][0]), 21)
            mG4_200 = TMarker((1.104 * float(m200G4[-1][0]) + 1.139), float(m200G4[j][0]), 21)
            mG4_240 = TMarker((1.175 * float(m240G4[-1][0]) - 9.41), float(m240G4[j][0]), 21)
        m.SetMarkerColor(2)
        m.SetMarkerSize(2)
        #m.Draw("psame")
        m_172.SetMarkerColor(2)
        m_172.SetMarkerSize(2)
        #m_172.Draw("psame")
        m_240.SetMarkerColor(2)
        m_240.SetMarkerSize(2)
        #m_240.Draw("psame")
        #leg.AddEntry(m,"Nominal Momenta","p")
        if graph.GetName() != "momentum" and graph.GetName() != "integralX":
            mG4_172.SetMarkerColor(4)
            mG4_172.SetMarkerSize(2)
            # mG4_172.Draw("psame")
            mG4_200.SetMarkerColor(4)
            mG4_200.SetMarkerSize(2)
            #mG4_200.Draw("psame")
            mG4_240.SetMarkerColor(4)
            mG4_240.SetMarkerSize(2)
            # mG4_240.Draw("psame")
            #leg.AddEntry(mG4_200,"Cobb Predictions","p")
        leg.Draw("same")
        j+=1
        i+=1
        t1 = TText(0.18,0.215,"MICE Preliminary")
        t2 = TText(0.18,0.185,"ISIS Cycle 2015/04")
        # t1 = TText(0.15,0.875,"MICE preliminary")
        # t2 = TText(0.15,0.825,self.desc[0] + ", " + self.desc[1][2:5]+" MeV/c, March 2016, MAUS v2.5")
        t1.SetNDC(1)
        t1.SetTextSize(0.04)
        t1.SetTextFont(42)
        t2.SetNDC(1)
        t2.SetTextSize(0.03)
        t2.SetTextFont(42)
        t1.Draw()
        t2.Draw()
        c.Print(graph.GetName()+".pdf")
        # print h.GetName(), "Mean Value = ", hist.GetMean(), "+/-", hist.GetMeanError(), \
        #     "RMS Value = ", hist.GetRMS(),"+/-",hist.GetRMSError()


def FillGraphs(graphs, pointlist, offset, scale):
    
    j = 0
    for point in pointlist:
        for i in range(2,len(point)):
            if graphs[i].GetName() == "momentum":
                graphs[i].SetPoint(j, 
                                   scale * float(point[i][0]) + offset,
                                   (float(point[1][0]) + float(point[0][0]))/2.)
                graphs[i].SetPointError(j, 
                                        scale * float(point[i][1]),
                                        (float(point[1][0]) - float(point[0][0]))/2.)
            else:
                graphs[i].SetPoint(j, (scale * float(point[16][0]) + offset),
                                   float(point[i][0]))
                # 
                graphs[i].SetPointError(j, math.sqrt(pow(scale * float(point[16][1]),2)*(1 + 0.2*0.2)),
                                        math.sqrt(pow(float(point[i][1]),2)))
        j+=1


        
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
    vals.append([next(x.prop("value") for x in doc.xpathEval("spec/cuts") \
                          if x.prop("name").find("TOF_ll") >= 0), 0.1])
    vals.append([next(x.prop("value") for x in doc.xpathEval("spec/cuts") \
                          if x.prop("name").find("TOF_ul") >= 0), 0.1])
    
    # try:
    thX_md = rf.Get("thetaX_measdataCobb")
    # except Exception as inst:
    #     print inst
    #     return vals
    # if thX_md.ClassName() != 'TH1D':
    #    return vals
    thX_md.Rebin(2)
    thY_md = rf.Get("thetaY_measdataCobb")
    thY_md.Rebin(2)
    th2_md = rf.Get("theta2Scatt_measdataCobb")
    thX_rc = rf.Get("thetaX_recoCobb")
    thX_rc.Rebin(2)
    thY_rc = rf.Get("thetaY_recoCobb")
    thY_md.Rebin(2)
    th2_rc = rf.Get("theta2Scatt_recoCobb")
    pcalc  = rf.Get("cor_mom")
    
    vals.append([thX_md.Integral(), math.sqrt(thX_md.Integral())])
    if thX_md.Integral() < 0:
        return vals
    else:
        limits = math.sqrt(th2_md.GetMean())
        thX_md.Fit("gaus","Q","",-limits,limits)
        vals.append([thX_md.GetFunction("gaus").GetParameter("Constant"), \
                     thX_md.GetFunction("gaus").GetParError(0)])
        vals.append([1000*thX_md.GetFunction("gaus").GetParameter("Mean"), \
                     1000*thX_md.GetFunction("gaus").GetParameter("Mean")/math.sqrt(thX_md.Integral())])
        vals.append([1000*thX_md.GetFunction("gaus").GetParameter("Sigma"), \
                     1000*thX_md.GetFunction("gaus").GetParameter("Sigma")/math.sqrt(thX_md.Integral())])
        thX_rc.Fit("gaus","Q","",-limits,limits)
        vals.append([1000*thX_rc.GetFunction("gaus").GetParameter("Mean"), \
                     1000*thX_rc.GetFunction("gaus").GetParameter("Mean")/math.sqrt(thX_md.Integral())])
        vals.append([1000*thX_rc.GetFunction("gaus").GetParameter("Sigma"), \
                     1000*thX_rc.GetFunction("gaus").GetParameter("Sigma")/math.sqrt(thX_md.Integral())])
        thY_md.Fit("gaus","Q","",-limits,limits)
        vals.append([thY_md.Integral(), math.sqrt(thY_md.Integral())])
        vals.append([thY_md.GetFunction("gaus").GetParameter("Constant"), \
                     thY_md.GetFunction("gaus").GetParError(0)])
        vals.append([1000*thY_md.GetFunction("gaus").GetParameter("Mean"), \
                     1000*thY_md.GetFunction("gaus").GetParameter("Mean")/math.sqrt(thX_md.Integral())])
        vals.append([1000*thY_md.GetFunction("gaus").GetParameter("Sigma"), \
                     1000*thY_md.GetFunction("gaus").GetParameter("Sigma")/math.sqrt(thX_md.Integral())])
        thY_rc.Fit("gaus","Q","",-limits,limits)
        vals.append([1000*thY_rc.GetFunction("gaus").GetParameter("Mean"), \
                     1000*thY_rc.GetFunction("gaus").GetParameter("Mean")/math.sqrt(thX_md.Integral())])
        vals.append([1000*thY_rc.GetFunction("gaus").GetParameter("Sigma"), \
                     1000*thY_rc.GetFunction("gaus").GetParameter("Sigma")/math.sqrt(thX_md.Integral())])
        # th2_md.GetXaxis.SetRangeUser
        vals.append([1000 * math.sqrt(th2_md.GetMean()/2.), \
                     1000 * th2_md.GetMeanError()/math.sqrt(th2_md.GetMean()/2.)])
        vals.append([1000 * math.sqrt(th2_rc.GetMean()/2.), \
                     1000 * th2_rc.GetMeanError()/math.sqrt(th2_rc.GetMean()/2.)])
        
        vals.append([pcalc.GetMean(), pcalc.GetMeanError()])
        #vals.append([pcalc.GetMean(), pcalc.GetRMS()])
        return vals


def ExtractParsG4(xmlfile):
    # print xmlfile
    doc = libxml2.parseFile(xmlfile)
    vals = []
    rootfile = next(x.prop("name") for x in doc.xpathEval("spec/file") \
                        if x.prop("id").find("outfile") >= 0)
    print "Check MC in ",rootfile
    vals.append([next(x.prop("value") for x in doc.xpathEval("spec/cuts") \
                          if x.prop("name").find("TOF_ll") >= 0), 0.1])
    vals.append([next(x.prop("value") for x in doc.xpathEval("spec/cuts") \
                          if x.prop("name").find("TOF_ul") >= 0), 0.1])
    
    rf = TFile(rootfile)
    thX_md = rf.Get("thetaX_measuredCobb")
    thX_md.Rebin(2)
    thY_md = rf.Get("thetaY_measuredCobb")
    thY_md.Rebin(2)
    th2_md = rf.Get("theta2Scatt_measuredCobb")
    thX_rc = rf.Get("thetaX_Cobb")
    thX_rc.Rebin(2)
    thY_rc = rf.Get("thetaY_Cobb")
    thY_md.Rebin(2)
    th2_rc = rf.Get("theta2Scatt_Cobb")
    # try:
    pcalc  = rf.Get("thetaScatt_measdata_vpCobb")
    # except:
    #    pcalc  = rf.Get("thetaScatt_measdata_vpCobb")

    vals.append([thX_md.Integral(), math.sqrt(thX_md.Integral())])
    if thX_md.Integral() < 2000:
        return vals
    else:
        limits = math.sqrt(th2_md.GetMean())
        thX_md.Fit("gaus","Q","",-limits,limits)
        vals.append([thX_md.GetFunction("gaus").GetParameter("Constant"), \
                     thX_md.GetFunction("gaus").GetParError(0)])
        vals.append([1000*thX_md.GetFunction("gaus").GetParameter("Mean"), \
                     1000*thX_md.GetFunction("gaus").GetParameter("Mean")/math.sqrt(thX_md.Integral())])
        vals.append([1000*thX_md.GetFunction("gaus").GetParameter("Sigma"), \
                     1000*thX_md.GetFunction("gaus").GetParameter("Sigma")/math.sqrt(thX_md.Integral())])
        thX_rc.Fit("gaus","Q","",-limits,limits)
        vals.append([1000*thX_rc.GetFunction("gaus").GetParameter("Mean"), \
                     1000*thX_rc.GetFunction("gaus").GetParameter("Mean")/math.sqrt(thX_md.Integral())])
        vals.append([1000*thX_rc.GetFunction("gaus").GetParameter("Sigma"), \
                     1000*thX_rc.GetFunction("gaus").GetParameter("Sigma")/math.sqrt(thX_md.Integral())])
        thY_md.Fit("gaus","Q","",-limits,limits)
        vals.append([thY_md.Integral(), math.sqrt(thY_md.Integral())])
        vals.append([thY_md.GetFunction("gaus").GetParameter("Constant"), \
                     thY_md.GetFunction("gaus").GetParError(0)])
        vals.append([1000*thY_md.GetFunction("gaus").GetParameter("Mean"), \
                     1000*thY_md.GetFunction("gaus").GetParameter("Mean")/math.sqrt(thX_md.Integral())])
        vals.append([1000*thY_md.GetFunction("gaus").GetParameter("Sigma"), \
                     1000*thY_md.GetFunction("gaus").GetParameter("Sigma")/math.sqrt(thX_md.Integral())])
        thY_rc.Fit("gaus","Q","",-limits,limits)
        vals.append([1000*thY_rc.GetFunction("gaus").GetParameter("Mean"), \
                     1000*thY_rc.GetFunction("gaus").GetParameter("Mean")/math.sqrt(thX_md.Integral())])
        vals.append([1000*thY_rc.GetFunction("gaus").GetParameter("Sigma"), \
                     1000*thY_rc.GetFunction("gaus").GetParameter("Sigma")/math.sqrt(thX_md.Integral())])
        # th2_md.GetXaxis.SetRangeUser
        vals.append([1000 * math.sqrt(th2_md.GetMean()/2.), \
                     1000 * th2_md.GetMeanError()/math.sqrt(th2_md.GetMean()/2.)])
        vals.append([1000 * math.sqrt(th2_rc.GetMean()/2.), \
                     1000 * th2_rc.GetMeanError()/math.sqrt(th2_rc.GetMean()/2.)])
        
        vals.append([pcalc.GetMean(), pcalc.GetMeanError()])
        return vals

                     

if __name__=="__main__":
'''



