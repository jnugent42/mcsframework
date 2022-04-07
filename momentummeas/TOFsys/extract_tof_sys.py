import libxml2
import sys, os, subprocess, math, ROOT
from ROOT import TH1D, TCanvas, TFile, TGraphErrors, TMarker, TF1, TLegend, TText


def plotFiles(meas172, measfile, meas240, syslist):
    
    hists  = []
    graphs = []

    m172 = ExtractPars(meas172)
    mval = ExtractPars(measfile)
    m240 = ExtractPars(meas240)
    m172G4 = ExtractParsG4(meas172)
    m200G4 = ExtractParsG4(measfile)
    m240G4 = ExtractParsG4(meas240)
    print mval, m172, m240
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
            sysval[16][0] *= 1.107
            sysval[16][0] += 1.05
            sysval[16][1] *= 1.107
        elif sysfile.find("LiHMu_3200") >= 0:
            pcomp = 200
            '''
p0                        =      1.13876   +/-   0.671571    
p1                        =      1.10475   +/-   0.00407707
            '''
            sysval[16][0] *= 1.104
            sysval[16][1] *= 1.104
            sysval[16][0] += 1.139
        elif sysfile.find("LiHMu_3240") >= 0:
            pcomp = 240
            '''
p0                        =     -9.41424   +/-   0.293949    
p1                        =      1.17466   +/-   0.00142602  
            '''
            sysval[16][0] *= 1.17
            sysval[16][1] *= 1.17
            sysval[16][0] -= 9.41
            # sysval[16][1] *= 0.780
        
        if sysval[2][0] > 2000.: # and pcomp < sysval[-1][0] + 10:
            pointlist.append(sysval)
            print sysval
            
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
    func = TF1("func",'[0]*13.6*sqrt(1 + 105.65 * 105.65/x/x)/(x)',150,280)
    # func = TF1("func",'[0]*13.6*105.65*105.65/sqrt(x*x + 105.65 * 105.65)/x/x',150,280)
    func0 = TF1("func0",'[0]*13.6*sqrt(1 + 105.65 * 105.65/x/x)/(x)',150,280)
    ftof = TF1("ftof",'[0]/x + [1]', 26.8, 29.4)
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
            leg = TLegend(0.15,0.65,0.45,0.89)
        else: 
            leg = TLegend(0.5,0.5,0.89,0.89)
        leg.SetLineColor(10)
        leg.AddEntry(graph,"Data","p")
        
        if graph.GetName() == "rawmeantheta2" or graph.GetName() == "meantheta2" \
           or graph.GetName() == "rawsigmaX" or graph.GetName() == "sigmaX" \
           or graph.GetName() == "rawsigmaY" or graph.GetName() == "sigmaY":
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
                dsdp  = (fup - fdn)/8.
                errsys = dpdt * 0.129  
                print k, j, pointlist[k][1][0], pointlist[k][0][0], x, y, fup, fdn, errsys, graph.GetErrorY(k), math.sqrt(pow(float(pointlist[k][i][1]),2)), \
                    math.sqrt(pow(float(pointlist[k][i][1]),2) + errsys*errsys)
                graph.SetPointError(k, math.sqrt(ex*ex + errsys*errsys), math.sqrt(ey*ey + 0.0*errsys*errsys))
                
            graph.Fit('func')
            
            leg.AddEntry(func, "Fit to #frac{13.6 a}{p#beta}", "l")
            fup = graph.GetFunction('func').Clone()
            fup.SetName('fup')
            fup.SetParameter(0, graph.GetFunction('func').GetParameter('p0') + \
                             graph.GetFunction('func').GetParError(0))
            # fup.SetParameter(1, graph.GetFunction('func').GetParameter('p1') + \
            #                 graph.GetFunction('func').GetParError(1))
            fup.SetLineColor(4)
            leg.AddEntry(fup, "Fit plus/minus error", "l")
            fdn = graph.GetFunction('func').Clone()
            fdn.SetName('fdn')
            fdn.SetParameter(0, graph.GetFunction('func').GetParameter('p0') - \
                             graph.GetFunction('func').GetParError(0))
            # fdn.SetParameter(1, graph.GetFunction('func').GetParameter('p1') - \
            #                  graph.GetFunction('func').GetParError(1))
            fdn.SetLineColor(4)
            # leg.AddEntry(fdn, "Fit minus stat. error", "l")
            fup.Draw('lsame')
            fdn.Draw('lsame')
            func0.SetParameter(0, 253*(1 + 0.036*math.log(0.253*0.253)))
            func0.SetLineColor(8)
            func0.Draw("lsame")
            leg.AddEntry(func0,"PDG prediction","l")
            
        if graph.GetName() == "momentum":
            graph.Fit('ftof')
            print "TOF for 172 MeV/c is ",graph.GetFunction('ftof').Eval(172)
            print "TOF for 200 MeV/c is ",graph.GetFunction('ftof').Eval(200)
            print "TOF for 240 MeV/c is ",graph.GetFunction('ftof').Eval(240)
            leg.AddEntry(ftof, "Fit to #frac{a}{p} + b", "l")
            
        graph.GetHistogram().GetXaxis().SetRangeUser(25,30)
        if j!=2 and j!=3 and j!=4 and j!=6 and j!=8 and j!=9 \
           and j!=10 and j!=12 and j!=len(graphs)-1:
            graph.GetHistogram().GetYaxis().SetRangeUser(10,28)
        elif j==4 or j==6 or j==10 or j==12:
            graph.GetHistogram().GetYaxis().SetRangeUser(-5,10)
            # graph.GetHistogram().GetYaxis().SetRangeUser(float(mval[j][0]) - 30*float(mval[j][1]),
            #                                             float(mval[j][0]) + 30*float(mval[j][1]))
        # m = TMarker((float(mval[1][0]) + float(mval[0][0]), float(mval[j][0]), 25)
        elif j==2:
            graph.GetHistogram().GetYaxis().SetRangeUser(0,10000)
            
        if graph.GetName() == "momentum":
            m = TMarker(1.104 * float(mval[j][0]) + 1.139, (float(mval[1][0]) + float(mval[0][0]))/2., 20)
            m_172 = TMarker(1.107 * float(m172[j][0]) + 1.105, (float(m172[1][0]) + float(m172[0][0]))/2., 20)
            m_240 = TMarker(1.175 * float(m240[j][0]) - 9.41, (float(m240[1][0]) + float(m240[0][0]))/2., 20)
            # mG4_172 = TMarker(scale * float(m172G4[j][0]) + offset, (float(m172G4[1][0]) + float(m172G4[0][0]))/2., 20)
            # mG4_200 = TMarker(scale * float(m200G4[j][0]) + offset, (float(m200G4[1][0]) + float(m200G4[0][0]))/2., 20)
            # mG4_240 = TMarker(scale * float(m240G4[j][0]) + offset, (float(m240G4[1][0]) + float(m240G4[0][0]))/2., 20)
        else:
            m = TMarker((1.104 * float(mval[-1][0]) + 1.139), float(mval[j][0]), 20)
            m_172 = TMarker((1.107 * float(m172[-1][0]) + 1.105), float(m172[j][0]), 20)
            m_240 = TMarker((1.175 * float(m240[-1][0]) - 9.41), float(m240[j][0]), 20)
            
            mG4_172 = TMarker((1.107 * float(m172G4[-1][0]) + 1.105), float(m172G4[j][0]), 21)
            mG4_200 = TMarker((1.104 * float(m200G4[-1][0]) + 1.139), float(m200G4[j][0]), 21)
            mG4_240 = TMarker((1.175 * float(m240G4[-1][0]) - 9.41), float(m240G4[j][0]), 21)
        m.SetMarkerColor(2)
        m.SetMarkerSize(2)
        m.Draw("psame")
        m_172.SetMarkerColor(2)
        m_172.SetMarkerSize(2)
        # m_172.Draw("psame")
        m_240.SetMarkerColor(2)
        m_240.SetMarkerSize(2)
        # m_240.Draw("psame")
        leg.AddEntry(m,"Nominal Momenta","p")
        if graph.GetName() != "momentum" and graph.GetName() != "integralX":
            mG4_172.SetMarkerColor(4)
            mG4_172.SetMarkerSize(2)
            # mG4_172.Draw("psame")
            mG4_200.SetMarkerColor(4)
            mG4_200.SetMarkerSize(2)
            mG4_200.Draw("psame")
            mG4_240.SetMarkerColor(4)
            mG4_240.SetMarkerSize(2)
            # mG4_240.Draw("psame")
            leg.AddEntry(mG4_200,"Moliere Predictions","p")
        leg.Draw("same")
        j+=1
        i+=1
        t1 = TText(0.18,0.215,"MICE isis cycle 2015/04")
        t2 = TText(0.18,0.185,"LiH, MAUS v2.5")
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
    thX_md = rf.Get("thetaX_measdataMoliere")
    # except Exception as inst:
    #     print inst
    #     return vals
    # if thX_md.ClassName() != 'TH1D':
    #    return vals
    thX_md.Rebin(2)
    thY_md = rf.Get("thetaY_measdataMoliere")
    thY_md.Rebin(2)
    th2_md = rf.Get("theta2Scatt_measdataMoliere")
    thX_rc = rf.Get("thetaX_recoMoliere")
    thX_rc.Rebin(2)
    thY_rc = rf.Get("thetaY_recoMoliere")
    thY_md.Rebin(2)
    th2_rc = rf.Get("theta2Scatt_recoMoliere")
    pcalc  = rf.Get("thetaScatt_measdata_vpMoliere")
    
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
    thX_md = rf.Get("thetaX_measuredMoliere")
    thX_md.Rebin(2)
    thY_md = rf.Get("thetaY_measuredMoliere")
    thY_md.Rebin(2)
    th2_md = rf.Get("theta2Scatt_measuredMoliere")
    thX_rc = rf.Get("thetaX_Moliere")
    thX_rc.Rebin(2)
    thY_rc = rf.Get("thetaY_Moliere")
    thY_md.Rebin(2)
    th2_rc = rf.Get("theta2Scatt_Moliere")
    # try:
    pcalc  = rf.Get("thetaScatt_measdata_vpMoliere")
    # except:
    #    pcalc  = rf.Get("thetaScatt_measdata_vpMoliere")

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
    plotFiles(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
