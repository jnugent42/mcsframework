from ROOT import TCanvas, TFile, TH1D, TLegend, TText, TLatex, TMath, TROOT, TMatrix, TMatrixDSym, TMatrixD
from ROOT import TH1D, TCanvas, TFile, TGraphErrors, TMarker, TF1, TLegend, TText, TLine, TROOT, TGraph, TPad
import sys, math
from math import sqrt
import pprint
import scipy.special as sc
import ROOT
import numpy

stops = [0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000]
red   = [0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764]
green = [0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832]
blue  = [0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539]
s = numpy.array(stops)
r = numpy.array(red)
g = numpy.array(green)
b = numpy.array(blue)
ncontours = 255
npoints = len(s)
ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
ROOT.gStyle.SetNumberContours(ncontours)
# axes and labels
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.15)
ROOT.gStyle.SetPadRightMargin(0.15)
for axis in "X", "Y":
    ROOT.gStyle.SetNdivisions(505, axis)
    ROOT.gStyle.SetLabelSize(0.05, axis)
    ROOT.gStyle.SetTitleSize(0.06, axis)
    ROOT.gStyle.SetTitleOffset(1.10, axis)
level=0.9
ROOT.TColor.GetColor(0.2082*level, 0.1664*level, 0.5293*level)

TROOT.gROOT.SetBatch(1)

class plotgen:

    def __init__(self, filename):
        self.fname = filename
        self.desc = self.fname.split("_")
        self.histvarnames = ['thetaX','thetaY','thetaScatt','theta2Scatt']
        # self.histstatenames = ['measdataCobb','refconv_GEANT','refconv_Cobb', 'refconv_Moliere', 'ref', 'GEANT']
        self.histstatenames = ['measdataCobb','refconv_GEANT','refconv_Cobb', 'refconv_Moliere']
        self.normnames = ['measdata', 'GEANT', 'recoGold','conv', 'ref', 'graph']
        self.histstatedesc = ['Raw Data','GEANT4 Model Default MCS', 'Carlisle-Cobb', 'Moliere']
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
#        hist.GetXaxis().SetLabelSize(0.05)
#        hist.GetXaxis().SetTitleSize(0.05)
#        hist.GetYaxis().SetLabelSize(0.05)
#        hist.GetYaxis().SetTitleSize(0.05)
        hist.SetLineColor(self.histcolors[i])
        hist.SetMarkerColor(self.histcolors[i])
        hist.SetLineWidth(3)
        if i > 0: hist.SetMarkerStyle(21+self.histcolors[i])
        else:  hist.SetMarkerStyle(20)
        hist.SetStats(0)


    def addToSysErr(self, i, histl, histh, base, scale, histvar, sysname, norms):
        # print sysname
        # print norms
        # if i == 0:
        #     norm = self.eventnorm[histvar]
        # else:
        #     norm = histl.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
        # if norm>0:
        #     histl.Scale(1./norm)
        # if i == 0:
        #     norm = self.eventnorm[histvar]
        # else:
        #     norm = histh.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
        # if norm>0:
        #     histh.Scale(1./norm)


        if sysname != "decon":
            norm = 0
            if norms[0]>0:
                histl.Scale(1./norms[0])
            norm = 0
            if norms[0]>0:
                histh.Scale(1./norms[0])


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
            # histl.GetXaxis().SetRangeUser(-0.06,0.06)
            # histh.GetXaxis().SetRangeUser(-0.06,0.06)
            histl.GetXaxis().SetRangeUser(-0.0465,0.045)
            histh.GetXaxis().SetRangeUser(-0.0465,0.045)
            histl.Fit("gaus","Q0","same",-0.045,0.045)
            histh.Fit("gaus","Q0","same",-0.045,0.045)
            self.RMSsysdiff[histvar][self.histstatenames[i]][sysname] \
                = 1000*(histh.GetFunction('gaus').GetParameter('Sigma') - \
                histl.GetFunction('gaus').GetParameter('Sigma'))
            self.RMSsyserr[histvar][self.histstatenames[i]][sysname] \
                = math.fabs(self.RMSsysdiff[histvar][self.histstatenames[i]][sysname]) * scale


    def addToRMS(self, i, hist, base, resplot, histvar, norms, pname):
        norm = 0
        # if i == 0:
        #     if self.useeventnorm:
        #         self.eventnorm[histvar]  = hist.GetEntries()
        #     self.h0integral[histvar] = hist.Integral()
        #     norm = self.eventnorm[histvar]
        # else:
        #     norm = hist.Integral() * self.eventnorm[histvar] / self.h0integral[histvar]
        # hist.Scale(1./norm)
        # print hist
        # print norms[i]
        # if (i==2):
        #     norms[i]=1e6
        if norms[i]:
            hist.Scale(1./norms[i])
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
            # hist.GetXaxis().SetRangeUser(-0.06,0.06)
            # resplot.GetXaxis().SetRangeUser(-0.06,0.06)
            hist.GetXaxis().SetRangeUser(-0.0465,0.045)
            resplot.GetXaxis().SetRangeUser(-0.0465,0.045)
            hist.Fit("gaus","Q0","same",-0.045,0.045)
            self.RMS[histvar][self.histstatenames[i]]    = (1000*hist.GetFunction('gaus').GetParameter('Sigma'))
            self.RMSErr[histvar][self.histstatenames[i]] = (1000*hist.GetFunction('gaus').GetParError(2))
            if (i == 2 and pname == "Result"):
                # self.RMSErr[histvar][self.histstatenames[i]] = (self.RMS[histvar][self.histstatenames[i]]/sqrt(norms[i]))
                self.RMSErr[histvar][self.histstatenames[i]] = 0

    def calculateChi2(self, i, hist, base, resplot, histvar, pname, norms):

        chi2 = 0
        cv = 0
        ndf = 0
        sysLowF = []
        sysLowH = []
        sysHighF = []
        sysHighH = []
        Icovar = TFile(self.fname).Get("cov_matrix")
        print self.fname
        Ucovar = TFile(self.fname).Get("cov_matrix")
        epsilon = TFile(self.fname).Get("epsilon")
        A = TFile(self.fname).Get("A")
        # A = TFile(self.fname).Get("rogersA")
        # TMatrixD(A)
        # epsilon = TFile(self.fname).Get("epsilon")
        # epsilon = TMatrixD(29,29)
        # epsilon = TMatrixD(A.Zero())
        # for lll in range(0,29):
        #     print epsilon
        #     print epsilon(lll,lll)
        #     # epsilon[lll][lll] = 0.01/1.1
            # epsilon[lll][lll] = 1
        # TMatrixD h(5,5);
        # TArrayD data(25);
        # for lll in range(0, 25):
        #      data[i] = 1./(ir+ic);
        # h.SetMatrixArray(data.GetArray());
        B = A

        # covar = TMatrix.UnitMatrix(A)

        # covar = Ucovar
        # Ucovar = Ucovar*0.03448275862

        Atran = TMatrixD(29,29)
        # Atran = A.Invert()
        Atran.Transpose(A)
        covar = A*epsilon*Atran
        # covar = A*epsilon
        Icovar = covar
        Icovar = A.Invert()
        # Icovar = epsilon.Invert()
        # if pname == "Result" and i == 2:
            # print Icovar.Determinant()
            # print "A"
            # A.Print()
            # print A.GetNrows()
            # print A.GetNcols()
            # print "epsilon"
            # epsilon.Print()
            # print "Atran"
            # Atran.Print()
            # Icovar.Print()
        # sysLowF  = [TFile(x[0]) for x in self.sysFiles]
        # sysLowH  = [f.Get(histvar + '_' + self.histstatenames[0]) for f in sysLowF]
        # owH  = [f.Get(histvar + '_' + self.histstatenames[i])
        for x in self.sysFiles:
            sysLowF.append(TFile(x[0]))
            f = TFile(x[0])
            if x[3]=="TOF" and histvar=="thetaX":
               plot = sysLowF[-1].Get("thetaY" + '_' + self.histstatenames[0])
               sysLowH.append(sysLowF[-1].Get("thetaY" + '_' + self.histstatenames[0]))
            else:
               plot = sysLowF[-1].Get(histvar + '_' + self.histstatenames[0])
               sysLowH.append(sysLowF[-1].Get(histvar + '_' + self.histstatenames[0]))

        # sysHighF = [TFile(x[1]) for x in self.sysFiles]
        # sysHighH = [f.Get(histvar + '_' + self.histstatenames[0]) for f in sysHighF]
        # sysHighH = [f.Get(histvar + '_' + self.histstatenames[i]) for f in sysHighF]
        for x in self.sysFiles:
            sysHighF.append(TFile(x[1]))
            f = TFile(x[1])
            if x[3]=="TOF" and histvar=="thetaX":
               plot = sysHighF[-1].Get("thetaY" + '_' + self.histstatenames[0])
               sysHighH.append(sysHighF[-1].Get("thetaY" + '_' + self.histstatenames[0]))
            else:
               plot = sysHighF[-1].Get(histvar + '_' + self.histstatenames[0])
               sysHighH.append(sysHighF[-1].Get(histvar + '_' + self.histstatenames[0]))
        sysScale = [x[2] for x in self.sysFiles]
        sysName = [x[3] for x in self.sysFiles]
        if (len(sysName)==6):
            if (sysName[5] == "decon"):
                sysLowH[5]  = sysLowF[5].Get(histvar + '_' + "graph")
                sysLowH[5].Scale(1./float(sysLowH[5].GetEntries()))
                sysHighH[5]  = sysHighF[5].Get(histvar + '_' + "recoGoldhold3")
                nsysHighH  = sysHighF[5].Get(histvar + '_' + "data").GetEntries()
                sysHighH[5].Scale(1./float(nsysHighH))
        syshists = []
        for q in range(len(self.sysFiles)):
            syshists.append(resplot.Clone())
            syshists[-1].SetName(resplot.GetName()+"_"+self.sysFiles[q][3])
            syshists[-1].GetYaxis().SetTitle("Upper - Lower Systematic Squared")


        self.RMSsysdiff[histvar][self.histstatenames[i]] = {}
        self.RMSsyserr[histvar][self.histstatenames[i]]  = {}
        self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum2"] = 0
        self.RMSsyserr[histvar][self.histstatenames[i]]["Sum2"] = 0
        self.RMSsysdiff[histvar][self.histstatenames[i]]["Sum"] = 0
        self.RMSsyserr[histvar][self.histstatenames[i]]["Sum"] = 0
        binmax = 0.045
        cloneb = TH1D("","",29,-binmax,binmax)
        cloneh = TH1D("","",29,-binmax,binmax)
        mcloneb = TH1D("","",47,-0.0705,0.0705)
        mcloneh = TH1D("","",47,-0.0705,0.0705)
        z = 1
        cv = 0
        firstbin = 99
        if histvar == 'theta2Scatt':
            binmax = 0.0036
        for k in range(len(self.sysFiles)):
            if pname=="Gold" and i==1:
                continue

            print sysLowF[k]
            self.formatHists(sysLowH[k], i)
            self.formatHists(sysHighH[k], i)
            # print self.sysFiles[k][0]
            self.addToSysErr(i, sysHighH[k], sysLowH[k], base, sysScale[k], histvar, self.sysFiles[k][3], norms)

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
            # print self.sysFiles[k][3]
            # print "Sum",self.RMSsyserr[histvar][self.histstatenames[i]]["Sum"]

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
                if math.fabs(hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()-23+j)) < binmax-0.003 :
                    mcloneh.SetBinContent(z,hist.GetBinContent(hist.GetMaximumBin()-23+j))
                    mcloneb.SetBinContent(z,base.GetBinContent(j))
                    if hist.GetBinContent(hist.GetMaximumBin()-23+j) > 0:
                        cv += res*res/hist.GetBinContent(hist.GetMaximumBin()-23+j)
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
                        # print "sys",sys
                        sys *= sysScale[k]

                        # if (sys*sys>1e-8):
                        syshists[k].SetBinContent(j, sys*sys)
                        # if pname == "Result" and k == 5 and i == 0 :
                        #      print j," & ",sys," \\"
                        #      if j == 1:
                        #         print "error from decon"
                        # syshists[k].SetBinContent(j, sysLowH[k].GetBinContent(j))
                        syshists[k].SetBinError(j, syshists[k].GetBinContent(j)*1e-3)
                        # err2 += sys*sys
                        if pname == "Result" and k == 5 :
                            err2hist += 0
                        else:
                            err2hist += sys*sys

            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax and err2hist > 0 and firstbin > j:
                 firstbin = j
            if i == 0:
                self.total = base.GetEntries()
            flag = 0
            if pname == "Result" and i == 2:
                flag = 1
            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax and err2hist > 0 and pname != "Result" :
                chi2 += res*res/err2hist
                flag = 0
                if i == 0:
                    self.central += base.GetBinContent(j)

                # for z in range(1,29):
                #     chi2 += res*covar(j,z)*(hist.GetBinContent(firstbin+z-1) - base.GetBinContent(firstbin+z-1))/err2hist
                # chi2 += res*res
                # print "chi2",chi2
                # print "res*res",res*res
                # print "err2hist",err2hist
            hist.SetBinError(j, sqrt(err2hist))

            if pname == "Gold" and i == 1:
                hist.SetBinError(j, sqrt(err2))
            if pname == "Result" and i == 1:
                if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax - 0.001 and err2hist > 0 :
                    chi2 += res*res/err2hist
                    # for z in range(0,29):
                        # z = j-10
                        # chi2 += res*Icovar(j-10,z)*(hist.GetBinContent(z+10) - base.GetBinContent(z+10))/err2hist
            if pname == "Result" and i == 2:
                # chi2 = 0
                if math.fabs(hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()-24+j)) < binmax - 0.001 and err2hist > 0 :
                    chi2 += res*res/err2hist
                        # print j
                    # for z in range(0,29):
                        # chi2 += res*Icovar(j-10,z)*res/err2hist
                        # chi2 += res*Icovar(j-10,z)*(hist.GetBinContent(hist.GetMaximumBin()-14+z) - base.GetBinContent(z+10))/err2hist

                        # chi2 += res*covar(j-10,z)*(hist.GetBinContent(hist.GetMaximumBin()-24+firstbin+z-1) - base.GetBinContent(firstbin+z-1))/err2hist
                        # z = j-10

                        # chi2 += res*covar(j-10,z)*res/err2hist
                    # chi2 += res*res
                    # print "chi2",chi2
                    # print "res",res
                    # print "err2hist",err2hist

                hist.SetBinError(hist.GetMaximumBin()-24+j, 0)

            if err2hist > 0:
                resplot.SetBinContent(j, res/sqrt(err2hist))
                # resplot.SetBinContent(j, 1e3*res/sqrt(err2hist))
                # resplot.SetBinContent(j, math.sqrt(err2hist))
                # resplot.SetBinContent(j, sqrt(err2hist))
            resplot.SetBinError(j, 1)

            if pname == "Gold" and i == 1:
                resplot.SetBinContent(j, res/(base.GetBinError(j)))
                resplot.SetBinError(j, 1)

            if pname == "Result" and i == 2:
                res = hist.GetBinContent(hist.GetMaximumBin()-24+j) - base.GetBinContent(j)
                if err2hist > 0:
                    resplot.SetBinContent(hist.GetMaximumBin()-24+j, res/math.sqrt(err2hist))
                    # resplot.SetBinContent(hist.GetMaximumBin()-24+j, math.sqrt(err2hist))
                if base.GetBinError(j) > 0:
                    resplot.SetBinError(hist.GetMaximumBin()-24+j, err2/base.GetBinError(j))
            # if pname == "Result" and i == 2:
            #     res = hist.GetBinContent(j) - base.GetBinContent(j)
            #     if err2hist > 0:
            #         resplot.SetBinContent(j, res/math.sqrt(err2hist))
                # if base.GetBinError(i) > 0:
                #     resplot.SetBinError(j, err2hist/base.GetBinError(i))
            if i != 0 and base.GetBinError(j) > 0:
                resplot.SetBinError(j, sqrt(err2)/base.GetBinError(j))
            if pname == "con" and i == 2:
                # if base.GetBinError(j) > 0:
                    resplot.SetBinError(j, err2)
            if math.fabs(hist.GetXaxis().GetBinCenter(j)) < binmax:
                ndf += 1
            if pname == "Result" and i == 2:
                if math.fabs(hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()-24+j)) < binmax - 0.001:
                    # print "here",hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()-24+j)
                    ndf += 1
            #print i, resplot.GetXaxis().GetBinCenter(j), res, chi2, err2, ndf
            #print i, resplot.GetXaxis().GetBinCenter(j), res, chi2, err2, ndf

        ndf += 2
        c = TCanvas(self.fname[:-5]+'s_'+histvar+'_c1')
        c.SetBottomMargin(0.15)
        c.SetTopMargin(0.075)

        if pname == "Result" and i == 0:
            # syshists[5].GetYaxis().SetRangeUser(1e-9,1e-4)
            syshists[5].Rebin(6)
            syshists[5].Draw()
            f1 = TF1("f1", "gaus")
            # def fline(x, par):
            #     if (x < 2.5):
            #         TF1::RejectPoint()
            #         return 0
            #     return par[0] + par[1]*x[0];
            # fline = TF1("fline", fline);
            # f1.SetParLimits(0,1e-5,1e-5)
            # f1.SetParLimits(1,1e-14,1e-12)
            # f1.SetParLimits(2,-4e-3,-2e-3)
            # f1.SetParameter(0,1e-5)
            # f1.FixParameter(1,0)
            # f1.FixParameter(2,-3e-3)
            syshists[5].Fit(f1,"M")
            # syshists[5].Fit(fline)
            xaxis = syshists[5].GetXaxis()
            for j in range(0, base.GetNbinsX()+1):
                binCenter = xaxis.GetBinCenter(j)
                fup = syshists[5].GetFunction("f1").Eval(binCenter)/6
                # print base.GetBinError(j)
                # print fup
                hist.SetBinError(j, sqrt(hist.GetBinError(j)**2+sqrt(fup**2)))
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
        histclone2.Scale(1/norms[i])
        baseclone2 = baseclone.Clone()
        baseclone2.Scale(1/norms[0])
        # pvalue_root =  baseclone.Chi2Test(histclone,"WW")
        # kol =  baseclone.KolmogorovTest(histclone)
        # chi2_root =  baseclone.Chi2Test(histclone,"CHI2 WW")
        # print pname
        # print i
        # print "final cv",cv
        # pvalue_root =  baseclone2.Chi2Test(histclone2,"WW")
        # kol =  baseclone2.KolmogorovTest(histclone2)
        # chi2_root =  baseclone2.Chi2Test(histclone2,"CHI2 WW")

        pvalue_root =  cloneb.Chi2Test(cloneh,"")
        kol =  hist.KolmogorovTest(base, "M")
        # mcloneh.Scale(norms[i])
        # mcloneb.Scale(norms[0])
        if pname == "Result" and i == 2:
            kol =  mcloneh.KolmogorovTest(mcloneb, "M")
            # print kol
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
        # PValue2 = sc.gammaincc(chi2/2, ndf/2)

        PValue = sc.gammaincc(ndf/2, chi2/2)

        # PValue = PValue1/math.gamma(ndf)
        # PValue = (1.0 - PValue2)
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
            # hists = [f.Get(histvar+'_'+x) for x in self.histstatenames]
            hists = []
            kkk = 0
            for x in self.histstatenames:
               # print self.histstatedesc
               # print self.histstatenames
               # print kkk
               if self.histstatedesc[kkk]=="LiH MC":
                  fMC = TFile("MC"+self.fname)
                  x = 'measdataGEANT'
                  hists.append(fMC.Get(histvar+'_'+x))
                  # print "MC"
                  # print hists
               elif self.histstatedesc[kkk]=="No Absorber MC":
                  fMC = TFile("MC"+self.fname)
                  x = 'ref'
                  hists.append(fMC.Get(histvar+'_'+x))
               else:
                  hists.append(f.Get(histvar+'_'+x))
                  # print "data"
                  # print hists
               kkk += 1
            # print hists
            count = 0
            norms = []
            jjj = 0
            for n in self.histstatenames:
                for x in self.normnames:
                 # print 'n',n
                 # print 'x',x
                 if self.histstatedesc[jjj]=="LiH MC":
                        fMCnorm = TFile("MC"+self.fname)
                        norms.append(fMCnorm.Get(x)[0]*3)
                        jjj += 1
                        # print "sdgbsfgbsfgbsfbgs"
                        break

                 if self.histstatedesc[jjj]=="No Absorber MC":
                        fMCnorm = TFile("MC"+self.fname)
                        norms.append(fMCnorm.Get(x)[0]*3)
                        jjj += 1
                        # print "sdgbsfgbsfgbsfbgs"
                        break

                 if ("conv" in n):
                        norms.append(1e6*3)
                        break
                 if (n in "Moliere" and len(norms)>1):
                        # print hists[2].Integral()
                        # norms.append((norms[0]*3e-3))
                        norms.append(hists[2].Integral()*3)
                        # norms.append(10000)
                        # break
                 if (x in n):
                    # print f.Get(x)[0]
                    # print 'n',n
                    # print 'x',x
                    jjj += 1

                    if (f.Get(x)[0]>0):
                        norms.append(f.Get(x)[0]*3)
                    else:
                        norms.append(1)
                    if (x=="conv"):
                        norms.append(f.Get("refconv")[0]*3)
                    break


            # print hists[0]
            # print norms
            hists[0].SetTitle("")
            # print hists
            resplots = [x.Clone() for x in hists]
            resplots[0].SetTitle('')
            resplots[0].GetYaxis().SetTitle("Normalized Residuals")

            leg = TLegend(0.18,0.6,0.75,0.8)
            leg.SetLineColor(10)
            self.central = 0
            self.total = 0
            c = TCanvas(self.fname[:-5]+'_'+histvar+'_c1')
            for i in range(len(self.histstatedesc)):
                hists[i].Sumw2()
                hists[i]*3
                print hists[i]
                self.formatHists(hists[i], i)
                self.formatHists(resplots[i], i)
                # print norms
                self.addToRMS(i, hists[i], hists[0], resplots[i], histvar, norms, pname)
                if histvar=='theta2Scatt':
                    # hists[i].GetYaxis().SetTitle('Probability per '+str("%.2f" % round(1000*1000*hists[i].GetXaxis().GetBinWidth(4),2))+' mrad^{2}')
                    hists[i].GetYaxis().SetTitle('Probability per mrad^{2}')
                else:
                    # hists[i].GetYaxis().SetTitle('Probability per '+str("%.2f" % round(1000*hists[i].GetXaxis().GetBinWidth(4),2))+' mrad')
                    hists[i].GetYaxis().SetTitle('Probability per mrad')
                    c.SetLeftMargin(0.15);
#                    hists[i].GetYaxis().SetTitleOffset(1.4)
                leg.AddEntry(hists[i], self.histstatedesc[i], self.histopts[i])
                # print hists[i]
                # print hists[0]
                self.calculateChi2(i, hists[i], hists[0], resplots[i], histvar, pname, norms)



            if self.desc[0] == 'XePion':
                t1 = TText(0.18,0.885,"MICE ISIS cycle 2015/03")
                t2 = TText(0.18,0.85,"Xe, "+self.desc[1][2:5]+", MAUS v3.3.2")
            elif "LiH" in self.histstatedesc[i]:
                t1 = TText(0.18,0.885,"MICE Internal")
                t6 = TText(0.18,0.84,"ISIS cycle 2015/04")
                t2 = TText(0.18,0.795,"LiH, "+self.desc[1][2:5]+" MeV/c, MAUS v3.3.2")
                t3 = TText(0.18,0.84,"Central, "+str((self.central/self.total)*100)+",%")
            elif "No Absorber" in self.histstatedesc[i]:
                t1 = TText(0.18,0.885,"MICE Internal")
                t6 = TText(0.18,0.84,"ISIS cycle 2015/04")
                t2 = TText(0.18,0.795,"No Absorber, "+self.desc[1][2:5]+" MeV/c, MAUS v3.3.2")
                t3 = TText(0.18,0.84,"Central, "+str((self.central/self.total)*100)+",%")
            elif "conv" in self.histstatedesc[i]:
                t1 = TText(0.18,0.885,"MICE Internal")
                t6 = TText(0.18,0.84,"ISIS cycle 2015/04")
                t2 = TText(0.18,0.795,"LiH, "+self.desc[1][2:5]+" MeV/c, MAUS v3.3.2")
                t3 = TText(0.18,0.84,"Central, "+str((self.central/self.total)*100)+",%")
            else:
                t1 = TText(0.18,0.885,"MICE Internal")
                t6 = TText(0.18,0.84,"ISIS cycle 2015/04")
                t2 = TText(0.18,0.795,"LiH, "+self.desc[1][2:5]+" MeV/c, MAUS v3.3.2")
                t3 = TText(0.18,0.84,"Central, "+str((self.central/self.total)*100)+",%")
            print self.histstatedesc[i]
            t1.SetNDC(1)
            t1.SetTextSize(0.04)
            t1.SetTextFont(42)
            t2.SetNDC(1)
            t2.SetTextSize(0.04)
            t2.SetTextFont(42)
            t6.SetNDC(1)
            t6.SetTextSize(0.04)
            t6.SetTextFont(42)
            hists[0].GetYaxis().SetRangeUser(4e-5,0.06)
            if "No Absorber" in self.histstatedesc[0]:
                hists[0].GetYaxis().SetRangeUser(4e-5,0.1)
            if self.histstatedesc[1]=="ref MC":
                hists[0].GetYaxis().SetRangeUser(4e-5,0.1)
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
                t6.Draw()
                t2.Draw()
                integral = hists[0].Integral()
                sumofw = hists[0].GetSumOfWeights()
                t4 = TText(0.18,0.8,"Sum of Weights "+str(sumofw))
                # t4 = TText(0.18,0.8,"Integral "+str(integral))
                t4.SetNDC(1)
                t4.SetTextSize(0.03)
                t4.SetTextFont(42)
                # t4.Draw()
                if histvar == 'theta2Scatt':
                   c.SetLogy()
                   # h.Fit("pol1","","",0.0,0.002)
                   # hists[0].SetTitle(";"+hists[0].GetXaxis().GetTitle()+";"+hists[0].GetYaxis().GetTitle())
                   hists[0].SetTitle("; #theta_{Scatt}^{2} (radians^{2});"+hists[0].GetYaxis().GetTitle())
                   # TF1  *flin = new TF1("flin","sin(x)/x",0,0.002)
                   hists[0].Fit("expo","","",0.0,0.002)
                   g = hists[0].GetListOfFunctions().FindObject("expo")
                   par =  g.GetParameters()
                   t5 = TText(0.7,0.885,"m="+str("%.2f" % par[1]))
                   t5.SetNDC(1)
                   t5.SetTextSize(0.04)
                   t5.SetTextFont(42)
                   t5.Draw()
                   h.GetYaxis().SetRangeUser(1e-4,1)
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
            # c.SetLogy(0)
            # resplots[0].GetYaxis().SetRangeUser(-2,2)
            resplots[0].SetTitle(";"+resplots[0].GetXaxis().GetTitle()+" (radians);"+resplots[0].GetYaxis().GetTitle())
            leg.SetX1NDC(0.5)
            leg.SetX2NDC(0.89)
            leg.SetY1NDC(0.2)
            leg.SetY2NDC(0.4)
            resplots[0].Draw("p")
            for r in resplots:
                r.GetYaxis().SetRangeUser(-10,10)
                # r.GetYaxis().SetRangeUser(-1e-2,1e-2)
                r.Draw('psame')
            leg.SetTextSize(0.04)
            leg.Draw('same')
            t1.Draw()
            t6.Draw()
            t2.Draw()
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_res_T.pdf')
            c.SaveAs(pname+'_'+self.fname[:-5]+'_'+histvar+'_sys_res_pq.jpg')

        # if pname=='Result':
        #     self.RMSsysdiff['thetaX']['recoGold_effi_only_gold']['TOF'] \
        #         = self.RMSsysdiff['thetaY']['recoGold_effi_only_gold']['TOF']
        #     self.RMSsyserr["thetaX"]['recoGold_effi_only_gold']['TOF'] \
        #         = self.RMSsyserr['thetaY'][self.histstatenames[i]]['TOF']
        #     self.RMSsysdiff['thetaX']['recoGold_effi_only_gold']["Sum2"] \
        #                 += self.RMSsysdiff['thetaY']['recoGold_effi_only_gold']['TOF'] *\
        #                 self.RMSsysdiff['thetaY']['recoGold_effi_only_gold']['TOF']
        #     self.RMSsyserr['thetaX']['recoGold_effi_only_gold']["Sum2"] \
        #                 += self.RMSsyserr['thetaY']['recoGold_effi_only_gold']['TOF']  *\
        #                 self.RMSsyserr['thetaY']['recoGold_effi_only_gold']['TOF']
        #     self.RMSsysdiff['thetaX']['recoGold_effi_only_gold']["Sum"] \
        #                 = sqrt( self.RMSsysdiff['thetaX']['recoGold_effi_only_gold']["Sum2"] )
        #     self.RMSsyserr['thetaX']['recoGold_effi_only_gold']["Sum"] \
        #                 = sqrt( self.RMSsyserr['thetaX']['recoGold_effi_only_gold']["Sum2"] )
        # if pname=='raw':
        #     self.RMSsysdiff['thetaX']['measdataCobb']['TOF'] \
        #         = self.RMSsysdiff['thetaY']['measdataCobb']['TOF']
        #     self.RMSsyserr["thetaX"]['measdataCobb']['TOF'] \
        #         = self.RMSsyserr['thetaY'][self.histstatenames[i]]['TOF']
        #     self.RMSsysdiff['thetaX']['measdataCobb']["Sum2"] \
        #                 += self.RMSsysdiff['thetaY']['measdataCobb']['TOF'] *\
        #                 self.RMSsysdiff['thetaY']['measdataCobb']['TOF']
        #     self.RMSsyserr['thetaX']['measdataCobb']["Sum2"] \
        #                 += self.RMSsyserr['thetaY']['measdataCobb']['TOF']  *\
        #                 self.RMSsyserr['thetaY']['measdataCobb']['TOF']
        #     self.RMSsysdiff['thetaX']['measdataCobb']["Sum"] \
        #                 = sqrt( self.RMSsysdiff['thetaX']['measdataCobb']["Sum2"] )
        #     self.RMSsyserr['thetaX']['measdataCobb']["Sum"] \
        #                 = sqrt( self.RMSsyserr['thetaX']['measdataCobb']["Sum2"] )
        momhist = f.Get("LiH_mom")
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

            # syssummary.append("LiH & "+str("%.2f" % round(mom[0],sigfig(mom[1])))+"$\pm$"+str("%.2f" % round(mom[1],sigfig(mom[1])))+\
            syssummary.append("LiH & "+str("%.2f" % mom[0])+\
                    " & "+str("%.2f" % difference0)+\
                    " & "+str("%.2f" % syserr0)+\
                    " & "+str("%.2f" % relerr0)+\
                    " & "+str("%.2f" % difference1)+\
                    " & "+str("%.2f" % syserr1)+\
                    " & "+str("%.2f" % relerr1)+\
                    # " & "+str("%.2f" % round(difference3,sigfig(difference3)))+\
                    # " & "+str("%.2f" % round(syserr3,sigfig(syserr3)))+\
                    # " & "+str("%.2f" % round(relerr3,sigfig(relerr3)))
                    "\\\\")

         # syssummary.append(str("%.2f" % round(mom[0],2))+"$\pm$"+str("%.2f" % round(mom[1],2))+\
         syssummary.append(str("%.2f" % mom[0])+\
                " & "+ str("%.2f" % rms[0])+"$\pm$"+str("%.2f" % rms[1])+\
                " & "+str("%.2f" % self.RMSsysdiff[self.histvarnames[0]][self.histstatenames[0]]['Sum'])+\
                " & "+str("%.2f" % self.RMSsyserr[self.histvarnames[0]][self.histstatenames[0]]['Sum'])+\
                " & "+str("%.2f" % (self.RMSsyserr[self.histvarnames[0]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[0]][self.histstatenames[0]]))+\
                " & "+str("%.2f" % self.RMSsysdiff[self.histvarnames[1]][self.histstatenames[0]]['Sum'])+\
                " & "+str("%.2f" % self.RMSsyserr[self.histvarnames[1]][self.histstatenames[0]]['Sum'])+\
                " & "+str("%.2f" % (self.RMSsyserr[self.histvarnames[1]][self.histstatenames[0]]['Sum']/self.RMS[self.histvarnames[1]][self.histstatenames[0]]))+\
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
        if pname == 'Gold' or pname == "MCdata"  or pname == "MCref" :
         for histvar in self.histvarnames:
             # summary.append(str("%.2f" % mom[0])+"$\pm$"+str("%.2f" % mom[1])+\
             summary.append(str("%.2f" % mom[0])+\
                         "& $\\"+histvar+"$ & "+str("%.2f" % self.RMS[histvar][self.histstatenames[0]])+ \
                          "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[0]])+ \
                          "$\pm$"+str("%.2f" % self.RMSsyserr[histvar][self.histstatenames[0]]["Sum"])+ \
                          " & "+str("%.2f" % self.RMS[histvar][self.histstatenames[1]])+ \
                          "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[1]])+ \
                          " & "+str("%.2f" % self.Chi2[histvar][self.histstatenames[1]][0])+ \
                          " / "+ str(self.Chi2[histvar][self.histstatenames[1]][1])+ \
                          # " & "+str("%.2f" % self.pvalue[histvar][self.histstatenames[1]][0],1))+ \
                          " & "+str("{:.2f}".format(self.pvalue[histvar][self.histstatenames[1]][0]))+ \
                     # "& $\ "+histvar+"$ & "+str("%.2f" % self.RMS[histvar][self.histstatenames[0]])+ \
                     # "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[0]])+ \
                     # "$\pm$"+str("%.2f" % self.RMSsyserr[histvar][self.histstatenames[0]]["Sum"])+ \
                     # " & "+str("%.2f" % self.RMS[histvar][self.histstatenames[1]])+ \
                     # "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[1]])+ \
                     # " & "+str("%.2f" % self.Chi2[histvar][self.histstatenames[1]][0],1)+ \
                     # # " / "+ str(self.Chi2[histvar][self.histstatenames[1]][1])+ \
                     "\\\\")
        else:
             for histvar in self.histvarnames:
                 # summary.append(str("%.2f" % mom[0])+"$\pm$"+str("%.2f" % mom[1])+\
                 summary.append(str("%.2f" % mom[0])+\
                         "& $\\"+histvar+"$ & "+str("%.2f" % self.RMS[histvar][self.histstatenames[0]])+ \
                          "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[0]])+ \
                          "$\pm$"+str("%.2f" % self.RMSsyserr[histvar][self.histstatenames[0]]["Sum"])+ \
                          " & "+str("%.2f" % self.RMS[histvar][self.histstatenames[1]])+ \
                          "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[1]])+ \
                          "$\pm$"+str("%.2f" % self.RMSsyserr[histvar][self.histstatenames[1]]["Sum"])+ \
                          " & "+str("%.2f" % self.Chi2[histvar][self.histstatenames[1]][0])+ \
                          " / "+ str(self.Chi2[histvar][self.histstatenames[1]][1])+ \
                          # " & "+str("%.2f" % self.pvalue[histvar][self.histstatenames[1]][0],1))+ \
                          " & "+str("{:.2f}".format(self.pvalue[histvar][self.histstatenames[1]][0]))+ \
                          # " & "+str("%.2f" % self.pvalue[histvar][self.histstatenames[1]][1])+ \
                          " & "+str("%.2f" % self.RMS[histvar][self.histstatenames[2]])+ \
                         "$\pm$"+str("%.2f" % self.RMSErr[histvar][self.histstatenames[2]])+ \
                         " & "+str("%.2f" % self.Chi2[histvar][self.histstatenames[2]][0])+ \
                         " / "+ str(self.Chi2[histvar][self.histstatenames[2]][1])+ \
                         # " & "+str("{:e}".format(self.pvalue[histvar][self.histstatenames[2]][0],1)))+ \
                         " & "+str("{:.2f}".format(self.pvalue[histvar][self.histstatenames[2]][0]))+ \
                         # " & "+str("%.2f" % self.pvalue[histvar][self.histstatenames[2]][1])+ \
                         "\\\\")

        f.Close()
        return [summary, syssummary]

if __name__ == '__main__':

    raw172 = plotgen("LihMuon_03172.root")
    raw172.sysFiles = [
                      ["../TOFsys/LiHMu_3200_tof_lim28.4918072213_u28.6918072213.root",
                        "../TOFsys/LiHMu_3200_tof_lim28.6918072213_u28.8918072213.root",
                        70./200.,"TOF"],
                      ["../angdef/angdef_0/LiHMu_3172_15.root",
                       "../angdef/angdef_0/LiHMu_3172_100.root",
                       1/1483,"angdef"],
                      ["MCLihMuon_03172.root",
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
                       ["../TOFsys/LiHMu_3200_tof_lim27.6918072213_u27.8918072213.root",
                        "../TOFsys/LiHMu_3200_tof_lim28.0918072213_u28.2918072213.root",
                        70./400.,"TOF"],
                       ["../angdef/angdef_0/LiHMu_3200_15.root",
                        "../angdef/angdef_0/LiHMu_3200_100.root",
                        1/1483,"angdef"],
                        ["MCLihMuon_03200.root",
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
                       ["../TOFsys/LiHMu_3200_tof_lim27.0918072213_u27.2918072213.root",
                        "../TOFsys/LiHMu_3200_tof_lim27.2918072213_u27.4918072213.root",
                        70./200.,"TOF"],
                       ["../angdef/angdef_0/LiHMu_3240_15.root",
                        "../angdef/angdef_0/LiHMu_3240_100.root",
                        1/1483,"angdef"],
                       ["MCLihMuon_03240.root",
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
    raw172.normnames = ['measdata', 'ref', 'conv']
    raw200.normnames = ['measdata', 'ref', 'conv']
    raw240.normnames = ['measdata', 'ref', 'conv']
    raw172.histstatenames = ['measdataGEANT','ref','refconv_GEANT']
    raw172.histstatedesc = ['LiH Data','No Absorber Data', 'GEANT4 conv. No Absorber Data']
    raw200.histstatenames = ['measdataGEANT','ref','refconv_GEANT']
    raw200.histstatedesc = ['LiH Data','No Absorber Data', 'GEANT4 conv. No Absorber Data']
    raw240.histstatenames = ['measdataGEANT','ref','refconv_GEANT']
    raw240.histstatedesc = ['LiH Data','No Absorber Data', 'GEANT4 conv. No Absorber Data']
    [ref172, refsys172] = raw172.MCSPlot("refG4")
    [ref200, refsys200] = raw200.MCSPlot("refG4")
    [ref240, refsys240] = raw240.MCSPlot("refG4")

    print "\n"
    raw172.normnames = ['ref', 'ref']
    raw200.normnames = ['ref', 'ref']
    raw240.normnames = ['ref', 'ref']
    raw172.histstatenames = ['ref','refMC']
    raw172.histstatedesc = ['No Absorber Data','No Absorber MC']
    raw200.histstatenames = ['ref','refMC']
    raw200.histstatedesc = ['No Absorber Data','No Absorber MC']
    raw240.histstatenames = ['ref','refMC']
    raw240.histstatedesc = ['No Absorber Data','No Absorber MC']
    [MCref172, MCrefsys172] = raw172.MCSPlot("MCref")
    [MCref200, MCrefsys200] = raw200.MCSPlot("MCref")
    [MCref240, MCrefsys240] = raw240.MCSPlot("MCref")

    print "\n"
    raw172.normnames = ['measdata', 'measdata']
    raw200.normnames = ['measdata', 'measdata']
    raw240.normnames = ['measdata', 'measdata']
    raw172.histstatenames = ['measdataGEANT','measdataMC']
    raw172.histstatedesc = ['LiH Data','LiH MC']
    raw200.histstatenames = ['measdataGEANT','measdataMC']
    raw200.histstatedesc = ['LiH Data','LiH MC']
    raw240.histstatenames = ['measdataGEANT','measdataMC']
    raw240.histstatedesc = ['LiH Data','LiH MC']
    [MCdata172, MCdatasys172] = raw172.MCSPlot("MCdata")
    [MCdata200, MCdatasys200] = raw200.MCSPlot("MCdata")
    [MCdata240, MCdatasys240] = raw240.MCSPlot("MCdata")

    rawnoerror172 = plotgen("LihMuon_03172.root")
    rawnoerror172.sysFiles = []

    rawnoerror200 = plotgen("LihMuon_03200.root")
    rawnoerror200.sysFiles = []

    rawnoerror240 = plotgen("LihMuon_03240.root")
    rawnoerror240.sysFiles = []

    print "\n"
    rawnoerror172.normnames = ['measdata', 'measdata', 'conv']
    rawnoerror200.normnames = ['measdata', 'measdata', 'conv']
    rawnoerror240.normnames = ['measdata', 'measdata', 'conv']
    rawnoerror172.histstatenames = ['measdataGEANT','measdataMC','refconv_GEANT']
    rawnoerror172.histstatedesc = ['LiH Data','LiH MC', 'GEANT4 conv. No Absorber Data']
    rawnoerror200.histstatenames = ['measdataGEANT','measdataMC','refconv_GEANT']
    rawnoerror200.histstatedesc = ['LiH Data','LiH MC', 'GEANT4 conv. No Absorber Data']
    rawnoerror240.histstatenames = ['measdataGEANT','measdataMC','refconv_GEANT']
    rawnoerror240.histstatedesc = ['LiH Data','LiH MC', 'GEANT4 conv. No Absorber Data']
    [MCdatanoerror172, MCdatanoerrorsys172] = rawnoerror172.MCSPlot("MCdatanoerror")
    [MCdatanoerror200, MCdatanoerrorsys200] = rawnoerror200.MCSPlot("MCdatanoerror")
    [MCdatanoerror240, MCdatanoerrorsys240] = rawnoerror240.MCSPlot("MCdatanoerror")

    print "\n"
    raw172.normnames = ['measdata', 'conv', 'conv']
    raw200.normnames = ['measdata', 'conv', 'conv']
    raw240.normnames = ['measdata', 'conv', 'conv']
    raw172.histstatenames = ['measdataGEANT','refconv_GEANT','refconv_Moliere']
    raw172.histstatedesc = ['LiH Data', 'GEANT4 conv. No Absorber Data','Moliere conv. No Absorber Data']
    raw200.histstatenames = ['measdataGEANT','refconv_GEANT','refconv_Moliere']
    raw200.histstatedesc = ['LiH Data', 'GEANT4 conv. No Absorber Data','Moliere conv. No Absorber Data']
    raw240.histstatenames = ['measdataGEANT','refconv_GEANT','refconv_Moliere']
    raw240.histstatedesc = ['LiH Data', 'GEANT4 conv. No Absorber Data','Moliere conv. No Absorber Data']
    [con172, consys172] = raw172.MCSPlot("con")
    [con200, consys200] = raw200.MCSPlot("con")
    [con240, consys240] = raw240.MCSPlot("con")

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


    print 'MCData'
    print '\hline'
    print '\hline'
    print MCdata172[0]
    print MCdata172[1]
    print '\hline'
    print MCdata200[0]
    print MCdata200[1]
    print '\hline'
    print MCdata240[0]
    print MCdata240[1]

    print '\hline'
    print '\hline'
    print MCdata172[3]
    print MCdata200[3]
    print MCdata240[3]
    print '\hline'
    print '\hline\n'

    print 'MCData no error'
    print '\hline'
    print '\hline'
    print MCdatanoerror172[0]
    print MCdatanoerror172[1]
    print '\hline'
    print MCdatanoerror200[0]
    print MCdatanoerror200[1]
    print '\hline'
    print MCdatanoerror240[0]
    print MCdatanoerror240[1]

    print '\hline'
    print '\hline'
    print MCdatanoerror172[3]
    print MCdatanoerror200[3]
    print MCdatanoerror240[3]
    print '\hline'
    print '\hline\n'

    print 'Raw sys summary'
    print '\hline\n'
    print '\hline'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    for i in range(len(raw172.sysFiles)):
        print raw172.sysFiles[i][3], sys200[i]
    print '\hline\n'
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print sys200[-1]
    print '\hline\n'
    print '\hline'

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

    raw172.sysFiles.append(
                      ["MCLihMuon_03172.root",
                       "MCLihMuon_03172.root",
                       1,"decon"]
                      )
    raw200.sysFiles.append(
                      ["MCLihMuon_03200.root",
                       "MCLihMuon_03200.root",
                       1,"decon"]
                      )
    raw240.sysFiles.append(
                      ["MCLihMuon_03240.root",
                       "MCLihMuon_03240.root",
                       1,"decon"]
                      )
    raw172.histvarnames = ['thetaX','thetaY']
    raw200.histvarnames = ['thetaX','thetaY']
    raw240.histvarnames = ['thetaX','thetaY']
    raw172.normnames = ['recoGold', 'GEANT', 'Moliere']
    raw200.normnames = ['recoGold', 'GEANT', 'Moliere']
    raw240.normnames = ['recoGold', 'GEANT', 'Moliere']
    raw172.histstatenames = ['recoGold', 'GEANT', 'Moliere']
    raw172.histstatedesc = ['Deconvolved Data', 'GEANT4 Model', 'Moliere Model']
    raw200.histstatenames = ['recoGold', 'GEANT', 'Moliere']
    raw200.histstatedesc = ['Deconvolved Data', 'GEANT4 Model', 'Moliere Model']
    raw240.histstatenames = ['recoGold', 'GEANT', 'Moliere']
    raw240.histstatedesc = ['Deconvolved Data', 'GEANT4 Model', 'Moliere Model']

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
        print ' \&  ',raw172.sysFiles[i][3], Goldsys172[i]
    print '\hline\n'
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print Goldsys172[-1]
    print '\hline\n'
    print '\hline'

    print '\hline\n'
    print '\hline'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    for i in range(len(raw172.sysFiles)):
        print ' \& ',raw172.sysFiles[i][3], Goldsys200[i]
    print '\hline\n'
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print Goldsys200[-1]
    print '\hline\n'
    print '\hline'

    print '\hline\n'
    print '\hline'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    for i in range(len(raw172.sysFiles)):
        print ' \& ',raw172.sysFiles[i][3], Goldsys240[i]
    print '\hline\n'
    print 'Sum'
    print ' \& $\\Delta\\theta_{X}$ & $\\Delta\\theta_{Y}$ & $\\langle\\theta_{Scatt}^{2}\\rangle$ \\\\'
    print Goldsys240[-1]
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





    Goldfile172 = plotgen("MCLihMuon_03172.root")
    Goldfile172.sysFiles = [
                      ["../TOFsys/LiHMu_3200_tof_lim28.4918072213_u28.691807221.root",
                        "../TOFsys/LiHMu_3200_tof_lim28.6918072213_u28.8918072213.root",
                        70./200.,"TOF"],
                        ["../angdef/angdef_0/LiHMu_3172_15.root",
                        "../angdef/angdef_0/LiHMu_3172_100.root",
                        1/1483,"angdef"],
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
                       ["../TOFsys/LiHMu_3200_tof_lim27.6918072213_u27.8918072213.root",
                        "../TOFsys/LiHMu_3200_tof_lim28.0918072213_u28.2918072213.root",
                        70./400.,"TOF"],
                        ["../angdef/angdef_0/LiHMu_3200_15.root",
                        "../angdef/angdef_0/LiHMu_3200_100.root",
                        1/1483,"angdef"],
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
    Goldfile240.sysFiles = [
                       ["../TOFsys/LiHMu_3200_tof_lim27.0918072213_u27.2918072213.root",
                        "../TOFsys/LiHMu_3200_tof_lim27.2918072213_u27.4918072213.root",
                        70./200.,"TOF"],
                        ["../angdef/angdef_0/LiHMu_3240_15.root",
                        "../angdef/angdef_0/LiHMu_3240_100.root",
                        1/1483,"angdef"],
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

    # Goldfile172.useeventnorm = False
    # Goldfile200.useeventnorm = False
    # Goldfile240.useeventnorm = False
    # Goldfile200.normnames = ['recoGold', 'graph']
    # Goldfile200.normnames = ['recoGold', 'graph']
    Goldfile200.normnames = ['graph','recoGold']
    Goldfile200.histvarnames = ['thetaX','thetaY']
    # Goldfile200.histstatenames = ['recoGoldhold3_sym', 'graph']
    Goldfile200.histstatenames = ['recoGold', 'graph']
    Goldfile200.histstatedesc = ['Deconvolved Data', 'Truth']
    [Goldout200, Goldoutsys200] = Goldfile200.MCSPlot("Gold")



