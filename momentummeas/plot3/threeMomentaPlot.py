from ROOT import TH1D, TH2D, TCanvas, TLegend, TFile, TText

f0 = TFile("../lih172/LihMuon_03172.root")
f1 = TFile("../lih200/LihMuon_03200.root")
f2 = TFile("../lih240/LihMuon_03240.root")
f3 = TFile("XePion_03240.root")


# fz0 = TFile("../ZeroMuon_3172.root")
# fz1 = TFile("../ZeroMuon_03200_LiHMC.root")
# fz2 = TFile("../ZeroMuon_03240_LiHMC.root")

h0 = f0.Get("thetaX_data")
h1 = f1.Get("thetaX_data")
h2 = f2.Get("thetaX_data")
h3 = f3.Get("thetaX_data")

hz0 = f0.Get("thetaX_ref")
hz1 = f1.Get("thetaX_ref")
hz2 = f2.Get("thetaX_ref")
hz3 = f3.Get("thetaX_ref")

ht0 = f0.Get("thetaX_GEANT")
ht1 = f1.Get("thetaX_GEANT")
ht2 = f2.Get("thetaX_GEANT")
ht3 = f3.Get("thetaX_truthsim")

hc0 = f0.Get("thetaX_Cobb")
hc1 = f1.Get("thetaX_Cobb")
hc2 = f2.Get("thetaX_Cobb")
hc3 = f3.Get("thetaX_Cobb")

h20 = f0.Get("dataUS_XY")
h21 = f1.Get("dataUS_XY")
h22 = f2.Get("dataUS_XY")

h20_px = h20.ProjectionX("h20_px")
h21_px = h21.ProjectionX("h21_px")
h22_px = h22.ProjectionX("h22_px")

m20 = f0.Get("dataUSref_XY")
m21 = f1.Get("dataUSref_XY")
m22 = f2.Get("dataUSref_XY")

m20_px = m20.ProjectionX("m20_px")
m21_px = m21.ProjectionX("m21_px")
m22_px = m22.ProjectionX("m22_px")

t1 = TText(0.18,0.885,"MICE isis cycle 2015/04")
t2 = TText(0.18,0.85,"LiH, MAUS v2.5")
t1.SetNDC(1)
t1.SetTextSize(0.04)
t1.SetTextFont(42)
t2.SetNDC(1)
t2.SetTextSize(0.03)
t2.SetTextFont(42)
                                                                                
c = TCanvas("c","Scattering")
c.SetLogy(0)
c.SetBottomMargin(0.125)
c.SetLeftMargin(0.15)
c.SetTopMargin(0.07)
c.SetRightMargin(0.05)
leg = TLegend(0.75,0.75,0.925,0.925)
leg.SetLineColor(10)
h0.SetLineColor(4)
h0.SetMarkerColor(4)
h0.SetLineStyle(9)
h0.SetMarkerStyle(20)
h0.SetLineWidth(2)

n = h0.GetEntries()
h0.Scale(1./n)
h0.SetStats(0)
h0.Rebin(1)
print h0.GetXaxis().GetBinWidth(20)
leg.AddEntry(h0,"172 MeV/c","p")
h0.SetTitle(";#Delta#theta_{X} (radians);Probability per mrad")
binwidth = h0.GetXaxis().GetBinWidth(50)
print binwidth
h0.Scale(0.001/binwidth)
h0.GetYaxis().SetTitleOffset(1.05)
h0.GetXaxis().SetRangeUser(-0.06,0.06)
h0.GetYaxis().SetRangeUser(5e-6,0.05)
h0.GetXaxis().SetLabelSize(0.05)
h0.GetXaxis().SetTitleSize(0.05)
h0.GetYaxis().SetLabelSize(0.05)
h0.GetYaxis().SetTitleSize(0.05)
h0.GetYaxis().SetTitleOffset(1.5)
h0.Draw("p")

h1.Scale(0.001/binwidth)
h1.SetLineColor(2)
h1.SetLineStyle(5)
h1.SetMarkerColor(2)
h1.SetMarkerStyle(25)
h1.SetLineWidth(2)
m = h1.GetEntries()
h1.Scale(1./m)
h1.Rebin(1)
leg.AddEntry(h1,"200 MeV/c","p")
h1.Draw("psame")


h2.Scale(0.001/binwidth)
h2.SetLineColor(1)
h2.SetLineStyle(1)
h2.SetLineColor(1)
# h2.SetMarkerStyle(22)
h2.SetLineWidth(2)
h2.SetMarkerColor(1)
h2.SetMarkerStyle(26)
q = h2.GetEntries()
h2.Scale(1./q)
h2.Rebin(1)
leg.AddEntry(h2,"240 MeV/c","p")
h2.Draw("psame")
t1.Draw()
t2.Draw()

leg.Draw("same")
c.Print("LiH_Raw_Combined.eps")
c.Print("LiH_Raw_Combined_pq.jpg")

c.Clear()
leg.Clear()
t3 = TText(0.18,0.885,"MICE isis cycle 2015/03")
t3.SetNDC(1)
t3.SetTextSize(0.04)
t3.SetTextFont(42)
t4 = TText(0.18,0.85,"Xenon, Pion Beam, MAUS v2.5")
t4.SetNDC(1)
t4.SetTextSize(0.03)
t4.SetTextFont(42)
h3.SetLineColor(4)
h3.SetMarkerColor(4)
h3.SetLineStyle(9)
h3.SetMarkerStyle(20)
h3.SetLineWidth(2)

n = h3.GetEntries()
h3.Scale(1./n)
h3.SetStats(0)
h3.Rebin(1)
print h3.GetXaxis().GetBinWidth(20)
leg.AddEntry(h3,"240 MeV/c","p")
h3.SetTitle(";#Delta#theta_{X} (radians);Probability per mrad")
binwidth = h3.GetXaxis().GetBinWidth(50)
print binwidth
h3.Scale(0.001/binwidth)
h3.GetYaxis().SetTitleOffset(1.05)
h3.GetXaxis().SetRangeUser(-0.06,0.06)
h3.GetYaxis().SetRangeUser(5e-6,0.05)
h3.GetXaxis().SetLabelSize(0.05)
h3.GetXaxis().SetTitleSize(0.05)
h3.GetYaxis().SetLabelSize(0.05)
h3.GetYaxis().SetTitleSize(0.05)
h3.GetYaxis().SetTitleOffset(1.5)
h3.Draw("p")
t3.Draw()
t4.Draw()
leg.Draw()
c.Print("Xe_Raw.eps")
c.Print("Xe_Raw_pq.jpg")

c.Clear()
leg.Clear()

t5 = TText(0.18,0.85,"Helium, Pion Beam, MAUS v2.5")
t5.SetNDC(1)
t5.SetTextSize(0.03)
t5.SetTextFont(42)
hz3.SetLineColor(4)
hz3.SetMarkerColor(4)
hz3.SetLineStyle(9)
hz3.SetMarkerStyle(20)
hz3.SetLineWidth(2)

n = hz3.GetEntries()
hz3.Scale(1./n)
hz3.SetStats(0)
hz3.Rebin(1)
print hz3.GetXaxis().GetBinWidth(20)
leg.AddEntry(hz3,"240 MeV/c","p")
hz3.SetTitle(";#Delta#theta_{X} (radians);Probability per mrad")
binwidth = hz3.GetXaxis().GetBinWidth(50)
print binwidth
hz3.Scale(0.001/binwidth)
hz3.GetYaxis().SetTitleOffset(1.05)
hz3.GetXaxis().SetRangeUser(-0.06,0.06)
hz3.GetYaxis().SetRangeUser(5e-6,0.05)
hz3.GetXaxis().SetLabelSize(0.05)
hz3.GetXaxis().SetTitleSize(0.05)
hz3.GetYaxis().SetLabelSize(0.05)
hz3.GetYaxis().SetTitleSize(0.05)
hz3.GetYaxis().SetTitleOffset(1.5)
hz3.Draw("p")
t3.Draw()
t5.Draw()
leg.Draw()
c.Print("He_Raw.eps")
c.Print("He_Raw_pq.jpg")

c.Clear()
leg.Clear()

# c.SetLogy()
hz0.SetLineColor(4)
hz0.SetLineStyle(9)
hz0.SetLineWidth(2)
hz0.SetMarkerColor(4)
hz0.SetMarkerStyle(20)
hz0.SetStats(0)
nz = hz0.GetEntries()
hz0.Scale(1./nz)
hz0.Rebin(1)
hz0.GetXaxis().SetRangeUser(-0.06,0.06)
# hz0.SetMaximum(0.2)
# 
print hz0.GetXaxis().GetBinWidth(20)
leg.AddEntry(hz0,"172 MeV/c","p")
hz0.SetTitle(";#Delta#theta_{X} (radians);Probability per mrad")
binwidth = hz0.GetXaxis().GetBinWidth(50)
print binwidth
hz0.Scale(0.001/binwidth)
hz0.GetYaxis().SetRangeUser(5e-6,0.05)
hz0.GetXaxis().SetLabelSize(0.05)
hz0.GetXaxis().SetTitleSize(0.05)
hz0.GetYaxis().SetLabelSize(0.05)
hz0.GetYaxis().SetTitleSize(0.05)
hz0.GetYaxis().SetTitleOffset(1.5)
hz0.Draw("p")


hz1.Scale(0.001/binwidth)
hz1.SetLineColor(2)
hz1.SetLineStyle(5)
hz1.SetLineWidth(2)
hz1.SetMarkerColor(2)
hz1.SetMarkerStyle(25)
mz = hz1.GetEntries()
hz1.Scale(1./mz)
hz1.Rebin(1)
leg.AddEntry(hz1,"200 MeV/c","p")
hz1.Draw("psame")


hz2.Scale(0.001/binwidth)
hz2.SetLineColor(1)
hz2.SetLineStyle(1)
hz2.SetLineWidth(2)
hz2.SetMarkerColor(1)
hz2.SetMarkerStyle(26)
qz = hz2.GetEntries()
hz2.Scale(1./qz)
hz2.Rebin(1)
leg.AddEntry(hz2,"240 MeV/c","p")
hz2.Draw("psame")

t6 = TText(0.18,0.85,"Zero Abs., Muon Beams, MAUS v2.5")
t6.SetNDC(1)
t6.SetTextSize(0.03)
t6.SetTextFont(42)
t1.Draw()
t6.Draw()
leg.Draw("same")
c.Print("Empty_Raw_Combined.eps")
c.Print("Empty_Raw_Combined_pq.jpg")
leg.Clear()

c.SetLogy()
ht0.SetLineColor(4)
ht0.SetLineStyle(9)
ht0.SetLineWidth(2)
ht0.SetStats(0)
nz = ht0.GetEntries()
ht0.Scale(1./nz)
ht0.Rebin(1)
ht0.GetXaxis().SetRangeUser(-0.06,0.06)
ht0.GetYaxis().SetRangeUser(5e-5,0.1)
print ht0.GetXaxis().GetBinWidth(20)
leg.AddEntry(ht0,"172 MeV/c","l")
ht0.SetTitle(";#Delta#theta_{X} (radians);Probability per mrad")
binwidth = ht0.GetXaxis().GetBinWidth(50)
print binwidth
ht0.Scale(0.001/binwidth)
ht0.GetXaxis().SetLabelSize(0.05)
ht0.GetXaxis().SetTitleSize(0.05)
ht0.GetYaxis().SetLabelSize(0.05)
ht0.GetYaxis().SetTitleSize(0.05)
ht0.GetYaxis().SetTitleOffset(1.5)
ht0.Draw()

ht1.SetLineColor(2)
ht1.SetLineStyle(5)
ht1.SetLineWidth(2)
mz = ht1.GetEntries()
ht1.Scale(1./mz)
ht1.Scale(0.001/binwidth)
ht1.Rebin(1)
leg.AddEntry(ht1,"200 MeV/c","l")
ht1.Draw("same")

ht2.SetLineColor(1)
ht2.SetLineStyle(1)
ht2.SetLineWidth(2)
qz = ht2.GetEntries()
ht2.Scale(1./qz)
ht2.Scale(0.001/binwidth)
ht2.Rebin(1)
leg.AddEntry(ht2,"240 MeV/c","l")
ht2.Draw("same")

t7 = TText(0.18,0.885,"MICE GEANT4 Simulation")
t7.SetNDC(1)
t7.SetTextSize(0.04)
t7.SetTextFont(42)
t8 = TText(0.18,0.85,"LiH, Muon Beams, MAUS v2.5")
t8.SetNDC(1)
t8.SetTextSize(0.03)
t8.SetTextFont(42)
t7.Draw()
t8.Draw()
leg.Draw("same")
c.Print("GEANT4_LiH_Combined.eps")
c.Print("GEANT4_LiH_Combined_pq.jpg")
leg.Clear()

ht3.SetLineColor(1)
ht3.SetLineStyle(1)
ht3.SetLineWidth(2)
qz = ht3.GetEntries()
ht3.Scale(1./qz)
ht3.Scale(0.001/binwidth)
ht3.Rebin(1)
ht3.GetXaxis().SetRangeUser(-0.06,0.06)
ht3.GetYaxis().SetRangeUser(5e-5,0.1)
ht3.SetStats(0)
leg.AddEntry(ht3,"240 MeV/c","l")
ht3.Draw()
t7.Draw()
t4.Draw()
leg.Draw("same")
c.Print("GEANT4_Xe.eps")
c.Print("GEANT4_Xe_pq.jpg")
leg.Clear()

hc0.SetLineColor(4)
hc0.SetLineStyle(9)
hc0.SetLineWidth(2)
hc0.SetStats(0)
nz = hc0.GetEntries()
hc0.Scale(1./nz)
hc0.Rebin(1)
hc0.GetXaxis().SetRangeUser(-0.06,0.06)
hc0.GetYaxis().SetRangeUser(5e-6,0.1)
print hc0.GetXaxis().GetBinWidth(20)
leg.AddEntry(hc0,"172 MeV/c","l")
hc0.SetTitle(";#Delta#theta_{X} (radians);Probability per mrad")
binwidth = hc0.GetXaxis().GetBinWidth(50)
print binwidth
hc0.Scale(0.001/binwidth)
hc0.GetXaxis().SetLabelSize(0.05)
hc0.GetXaxis().SetTitleSize(0.05)
hc0.GetYaxis().SetLabelSize(0.05)
hc0.GetYaxis().SetTitleSize(0.05)
hc0.GetYaxis().SetTitleOffset(1.5)
hc0.Draw()

hc1.SetLineColor(2)
hc1.SetLineStyle(5)
hc1.SetLineWidth(2)
mz = hc1.GetEntries()
hc1.Scale(1./mz)
hc1.Scale(0.001/binwidth)
hc1.Rebin(1)
leg.AddEntry(hc1,"200 MeV/c","l")
hc1.Draw("same")

hc2.SetLineColor(1)
hc2.SetLineStyle(1)
hc2.SetLineWidth(2)
qz = hc2.GetEntries()
hc2.Scale(1./qz)
hc2.Scale(0.001/binwidth)
hc2.Rebin(1)
leg.AddEntry(hc2,"240 MeV/c","l")
hc2.Draw("same")


t9 = TText(0.18,0.885,"MICE Carlisle-Cobb Simulation")
t9.SetNDC(1)
t9.SetTextSize(0.04)
t9.SetTextFont(42)
t10 = TText(0.18,0.85,"LiH, Muon Beams, MAUS v2.5")
t10.SetNDC(1)
t10.SetTextSize(0.03)
t10.SetTextFont(42)
t9.Draw()
t10.Draw()
t9.Draw()
t10.Draw()
leg.Draw("same")
c.Print("Cobb-Carlisle_Combined.eps")
c.Print("Cobb-Carlisle_Combined_pq.jpg")
leg.Clear()

hc3.SetLineColor(1)
hc3.SetLineStyle(1)
hc3.SetLineWidth(2)
qz = hc3.GetEntries()
hc3.Scale(1./qz)
hc3.Scale(0.001/binwidth)
hc3.Rebin(1)
hc3.GetXaxis().SetRangeUser(-0.06,0.06)
hc3.GetYaxis().SetRangeUser(5e-5,0.1)
hc3.SetStats(0)
leg.AddEntry(hc3,"240 MeV/c","l")
hc3.Draw()
t9.Draw()
t4.Draw()
leg.Draw("same")
c.Print("Carlisle-Cobb_Xe.eps")
c.Print("Carlisle-Cobb_Xe_pq.jpg")
leg.Clear()


c.Clear()
c.SetLogy(0)
leg.Clear()
leg.AddEntry(h20_px,"Data","l")
leg.AddEntry(m20_px,"MC","l")
h20_px.SetLineColor(1)
m20_px.SetLineColor(2)
m20_px.SetLineStyle(2)
m20_px.GetXaxis().SetLabelSize(0.05)
m20_px.GetXaxis().SetTitleSize(0.05)
m20_px.GetYaxis().SetLabelSize(0.05)
m20_px.GetYaxis().SetTitleSize(0.05)
n20 = h20_px.GetEntries()
h20_px.Scale(1./n20)
l20 = m20_px.GetEntries()
m20_px.Scale(1./l20)
m20_px.Draw()
h20_px.Draw("same")
leg.Draw("same")
# c.Print("Momentum_Dist_172MeV.eps")

c.Clear()
leg.Clear()
leg.AddEntry(h21_px,"Data","l")
leg.AddEntry(m21_px,"MC","l")
h21_px.SetLineColor(1)
m21_px.SetLineColor(2)
m21_px.SetLineStyle(2)
m21_px.GetXaxis().SetLabelSize(0.05)
m21_px.GetXaxis().SetTitleSize(0.05)
m21_px.GetYaxis().SetLabelSize(0.05)
m21_px.GetYaxis().SetTitleSize(0.05)
n21 = h21_px.GetEntries()
h21_px.Scale(1./n21)
l21 = m21_px.GetEntries()
m21_px.Scale(1./l21)
m21_px.Draw()
h21_px.Draw("same")
leg.Draw("same")
# c.Print("Momentum_Dist_200MeV.eps")

c.Clear()
leg.Clear()
leg.AddEntry(h22_px,"Data","l")
leg.AddEntry(m22_px,"MC","l")
h22_px.SetLineColor(1)
m22_px.SetLineColor(2)
m22_px.SetLineStyle(2)
m22_px.GetXaxis().SetLabelSize(0.05)
m22_px.GetXaxis().SetTitleSize(0.05)
m22_px.GetYaxis().SetLabelSize(0.05)
m22_px.GetYaxis().SetTitleSize(0.05)
n22 = h22_px.GetEntries()
h22_px.Scale(1./n22)
l22 = m22_px.GetEntries()
m22_px.Scale(1./l22)
m22_px.Draw()
h22_px.Draw("same")
leg.Draw("same")
# c.Print("Momentum_Dist_240MeV.eps")


c.Clear()
leg.Clear()
tof01_0 = f0.Get("tof10")
# tof01_0.SetLabel(2,"US Track Found")
tof01_1 = f1.Get("tof10")
tof01_2 = f2.Get("tof10")
tof01_0.SetFillColor(2)
tof01_1.SetFillColor(4)
tof01_2.SetFillColor(8)
tof01_0.SetLineColor(2)
tof01_1.SetLineColor(4)
tof01_2.SetLineColor(8)
tof01_0.SetFillStyle(3003)
tof01_1.SetFillStyle(3004)
tof01_2.SetFillStyle(3005)
leg.AddEntry(tof01_0,"172 MeV/c","f")
leg.AddEntry(tof01_1,"200 MeV/c","f")
leg.AddEntry(tof01_2,"240 MeV/c","f")
# tof01_0.SetTitle("LiH Abs., Muon Beams, March 2016, MAUS v2.5")
tof01_0.GetXaxis().SetRangeUser(24,33)
#tof01_0.GetYaxis().SetRangeUser(0,1e5)
tof01_0.GetYaxis().SetTitle("Events per 0.2 ns")
tof01_0.GetXaxis().SetLabelSize(0.05)
tof01_0.GetXaxis().SetTitleSize(0.05)
tof01_0.GetYaxis().SetLabelSize(0.05)
tof01_0.GetYaxis().SetTitleSize(0.05)
tof01_0.SetStats(0)
tof01_0.Draw()
tof01_1.Draw("same")
tof01_2.Draw("same")
leg.Draw()
# t1.Draw()
t2.Draw()
c.Print("TOF10_dist.eps")
c.Print("TOF10_dist_pq.jpg")

leg.Clear()
c.SetLogy()
c.SetBottomMargin(0.15)
accept_0 = f0.Get("cuts_accept")
accept_1 = f1.Get("cuts_accept")
accept_2 = f2.Get("cuts_accept")
accept_0.SetFillColor(2)
accept_1.SetFillColor(4)
accept_2.SetFillColor(8)
accept_0.SetLineColor(2)
accept_1.SetLineColor(4)
accept_2.SetLineColor(8)
accept_0.SetFillColor(2)
accept_1.SetFillColor(4)
accept_2.SetFillColor(8)
accept_0.SetBarWidth(0.25)
accept_1.SetBarWidth(0.25)
accept_2.SetBarWidth(0.25)
accept_0.SetBarOffset(0.0)
accept_1.SetBarOffset(0.3)
accept_2.SetBarOffset(0.6)
accept_0.SetMinimum(1e3)
accept_0.SetMaximum(8e6)
accept_0.GetXaxis().SetRange(1,4)
# accept_0.SetTitle("LiH Abs., Muon Beams, March 2016, MAUS v2.5")
leg.AddEntry(accept_0,"172 MeV/c","f")
leg.AddEntry(accept_1,"200 MeV/c","f")
leg.AddEntry(accept_2,"240 MeV/c","f")
accept_0.GetXaxis().SetLabelSize(0.05)
accept_0.GetXaxis().SetTitleSize(0.05)
accept_0.GetYaxis().SetLabelSize(0.05)
accept_0.GetYaxis().SetTitleSize(0.05)
accept_0.GetXaxis().SetTitleOffset(1.25)
accept_0.SetStats(0)
accept_0.Draw("bar")
accept_1.Draw("barsame")
accept_2.Draw("barsame")
leg.Draw()
# t1.Draw()
t2.Draw()
c.Print("LiH_Accept_Cuts.eps")
c.Print("LiH_Accept_Cuts_pq.jpg")
c.SetLogy(0)

c.Clear()

# c.SetTopMargin(0.9)
# c.SetLogz()

c.SetRightMargin(0.15)
# c.SetBottomMargin(0.10)
XY_10 = f1.Get("dataUS_alltof_XY")
# XY_10.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_10.SetStats(0)
XY_10.SetTitle("")
# XY_10.SetMaximum(600)
XY_10.SetMinimum(50)
XY_10.GetXaxis().SetLabelSize(0.05)
XY_10.GetXaxis().SetTitleSize(0.05)
XY_10.GetYaxis().SetLabelSize(0.05)
XY_10.GetYaxis().SetTitleSize(0.05)
XY_10.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUS_200_alltof_XY.eps")
c.Print("dataUS_200_alltof_XY_pq.jpg")

c.Clear()
XY_11 = f1.Get("dataUS_prerad_XY")
# XY_11.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_11.SetStats(0)
XY_11.SetTitle("")
# XY_11.SetMaximum(600)
XY_11.SetMinimum(2)
XY_11.GetXaxis().SetLabelSize(0.05)
XY_11.GetXaxis().SetTitleSize(0.05)
XY_11.GetYaxis().SetLabelSize(0.05)
XY_11.GetYaxis().SetTitleSize(0.05)
XY_11.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUS_200_prerad_XY.eps")
c.Print("dataUS_200_prerad_XY_pq.jpg")

c.Clear()
XY_12 = f1.Get("dataUS_XY")
# XY_12.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_12.SetStats(0)
XY_12.SetTitle("")
# XY_12.SetMaximum(600)
XY_12.SetMinimum(0)
XY_12.GetXaxis().SetLabelSize(0.05)
XY_12.GetXaxis().SetTitleSize(0.05)
XY_12.GetYaxis().SetLabelSize(0.05)
XY_12.GetYaxis().SetTitleSize(0.05)
XY_12.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUS_200_XY.eps")
c.Print("dataUS_200_XY_pq.jpg")

c.Clear()
dXdzdYdz_10 = f1.Get("dataUS_alltof_dXdzdYdz")
# dXdzdYdz_10.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_10.SetStats(0)
dXdzdYdz_10.SetTitle("")
# dXdzdYdz_10.SetMaximum(800)
dXdzdYdz_10.SetMinimum(20)
dXdzdYdz_10.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_10.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_10.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_10.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_10.GetXaxis().SetNdivisions(508)
dXdzdYdz_10.GetYaxis().SetNdivisions(508)
dXdzdYdz_10.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUS_200_alltof_dXdzdYdz.eps")
c.Print("dataUS_200_alltof_dXdzdYdz_pq.jpg")

c.Clear()
dXdzdYdz_11 = f1.Get("dataUS_prerad_dXdzdYdz")
# dXdzdYdz_11.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5)"
dXdzdYdz_11.SetStats(0)
dXdzdYdz_11.SetTitle("")
# dXdzdYdz_11.SetMaximum(800)
dXdzdYdz_11.SetMinimum(3)
dXdzdYdz_11.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_11.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_11.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_11.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_11.GetXaxis().SetNdivisions(508)
dXdzdYdz_11.GetYaxis().SetNdivisions(508)
dXdzdYdz_11.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUS_200_prerad_dXdzdYdz.eps")
c.Print("dataUS_200_prerad_dXdzdYdz_pq.jpg")

c.Clear()
dXdzdYdz_12 = f1.Get("dataUS_dXdzdYdz")
# dXdzdYdz_12.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_12.SetStats(0)
dXdzdYdz_12.SetTitle("")
# dXdzdYdz_12.SetMaximum(800)
dXdzdYdz_12.SetMinimum(0)
dXdzdYdz_12.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_12.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_12.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_12.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_12.GetXaxis().SetNdivisions(508)
dXdzdYdz_12.GetYaxis().SetNdivisions(508)
dXdzdYdz_12.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUS_200_dXdzdYdz.eps")
c.Print("dataUS_200_dXdzdYdz_pq.jpg")
# ================================================================


c.SetRightMargin(0.15)
# c.SetBottomMargin(0.10)
XY_10 = f1.Get("dataDS_alltof_XY")
# XY_10.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_10.SetStats(0)
XY_10.SetTitle("")
# XY_10.SetMaximum(600)
XY_10.SetMinimum(50)
XY_10.GetXaxis().SetLabelSize(0.05)
XY_10.GetXaxis().SetTitleSize(0.05)
XY_10.GetYaxis().SetLabelSize(0.05)
XY_10.GetYaxis().SetTitleSize(0.05)
XY_10.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDS_200_alltof_XY.eps")
c.Print("dataDS_200_alltof_XY_pq.jpg")

c.Clear()
XY_11 = f1.Get("dataDS_prerad_XY")
# XY_11.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_11.SetStats(0)
XY_11.SetTitle("")
# XY_11.SetMaximum(600)
XY_11.SetMinimum(2)
XY_11.GetXaxis().SetLabelSize(0.05)
XY_11.GetXaxis().SetTitleSize(0.05)
XY_11.GetYaxis().SetLabelSize(0.05)
XY_11.GetYaxis().SetTitleSize(0.05)
XY_11.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDS_200_prerad_XY.eps")
c.Print("dataDS_200_prerad_XY_pq.jpg")

c.Clear()
XY_12 = f1.Get("dataDS_XY")
# XY_12.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_12.SetStats(0)
XY_12.SetTitle("")
# XY_12.SetMaximum(600)
XY_12.SetMinimum(0)
XY_12.GetXaxis().SetLabelSize(0.05)
XY_12.GetXaxis().SetTitleSize(0.05)
XY_12.GetYaxis().SetLabelSize(0.05)
XY_12.GetYaxis().SetTitleSize(0.05)
XY_12.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDS_200_XY.eps")
c.Print("dataDS_200_XY_pq.jpg")

c.Clear()
dXdzdYdz_10 = f1.Get("dataDS_alltof_dXdzdYdz")
# dXdzdYdz_10.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_10.SetStats(0)
dXdzdYdz_10.SetTitle("")
# dXdzdYdz_10.SetMaximum(800)
dXdzdYdz_10.SetMinimum(20)
dXdzdYdz_10.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_10.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_10.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_10.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_10.GetXaxis().SetNdivisions(508)
dXdzdYdz_10.GetYaxis().SetNdivisions(508)
dXdzdYdz_10.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDS_200_alltof_dXdzdYdz.eps")
c.Print("dataDS_200_alltof_dXdzdYdz_pq.jpg")

c.Clear()
dXdzdYdz_11 = f1.Get("dataDS_prerad_dXdzdYdz")
# dXdzdYdz_11.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5)"
dXdzdYdz_11.SetStats(0)
dXdzdYdz_11.SetTitle("")
# dXdzdYdz_11.SetMaximum(800)
dXdzdYdz_11.SetMinimum(3)
dXdzdYdz_11.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_11.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_11.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_11.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_11.GetXaxis().SetNdivisions(508)
dXdzdYdz_11.GetYaxis().SetNdivisions(508)
dXdzdYdz_11.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDS_200_prerad_dXdzdYdz.eps")
c.Print("dataDS_200_prerad_dXdzdYdz_pq.jpg")

c.Clear()
dXdzdYdz_12 = f1.Get("dataDS_dXdzdYdz")
# dXdzdYdz_12.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_12.SetStats(0)
dXdzdYdz_12.SetTitle("")
# dXdzdYdz_12.SetMaximum(800)
dXdzdYdz_12.SetMinimum(0)
dXdzdYdz_12.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_12.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_12.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_12.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_12.GetXaxis().SetNdivisions(508)
dXdzdYdz_12.GetYaxis().SetNdivisions(508)
dXdzdYdz_12.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDS_200_dXdzdYdz.eps")
c.Print("dataDS_200_dXdzdYdz_pq.jpg")

# ====================================================



t2 = TText(0.18,0.85,"Zero Abs., February 2016, MAUS v2.5")
t2.SetNDC(1)
t2.SetTextSize(0.03)
t2.SetTextFont(42)

c.Clear()
XY_10 = f1.Get("dataUSref_alltof_XY")
# XY_10.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_10.SetStats(0)
XY_10.SetTitle("")
# XY_10.SetMaximum(400)
XY_10.SetMinimum(10)
XY_10.GetXaxis().SetLabelSize(0.05)
XY_10.GetXaxis().SetTitleSize(0.05)
XY_10.GetYaxis().SetLabelSize(0.05)
XY_10.GetYaxis().SetTitleSize(0.05)
XY_10.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUSref_200_alltof_XY.eps")
c.Print("dataUSref_200_alltof_XY_pq.jpg")

c.Clear()
XY_11 = f1.Get("dataUSref_prerad_XY")
#XY_11.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_11.SetStats(0)
XY_11.SetTitle("")
# XY_11.SetMaximum(400)
XY_11.SetMinimum(2)
XY_11.GetXaxis().SetLabelSize(0.05)
XY_11.GetXaxis().SetTitleSize(0.05)
XY_11.GetYaxis().SetLabelSize(0.05)
XY_11.GetYaxis().SetTitleSize(0.05)
XY_11.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUSref_200_prerad_XY.eps")
c.Print("dataUSref_200_prerad_XY_pq.jpg")

c.Clear()
XY_12 = f1.Get("dataUSref_XY")
# XY_12.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_12.SetStats(0)
XY_12.SetTitle("")
# XY_12.SetMaximum(400)
XY_12.SetMinimum(1)
XY_12.GetXaxis().SetLabelSize(0.05)
XY_12.GetXaxis().SetTitleSize(0.05)
XY_12.GetYaxis().SetLabelSize(0.05)
XY_12.GetYaxis().SetTitleSize(0.05)
XY_12.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUSref_200_XY.eps")
c.Print("dataUSref_200_XY_pq.jpg")

c.Clear()
dXdzdYdz_10 = f1.Get("dataUSref_alltof_dXdzdYdz")
# dXdzdYdz_10.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_10.SetStats(0)
dXdzdYdz_10.SetTitle("")
# dXdzdYdz_10.SetMaximum(300)
dXdzdYdz_10.SetMinimum(10)
dXdzdYdz_10.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_10.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_10.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_10.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_10.GetXaxis().SetNdivisions(508)
dXdzdYdz_10.GetYaxis().SetNdivisions(508)
dXdzdYdz_10.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUSref_200_alltof_dXdzdYdz.eps")
c.Print("dataUSref_200_alltof_dXdzdYdz_pq.jpg")

c.Clear()
dXdzdYdz_11 = f1.Get("dataUSref_prerad_dXdzdYdz")
# dXdzdYdz_11.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_11.SetStats(0)
dXdzdYdz_11.SetTitle("")
# dXdzdYdz_11.SetMaximum(300)
dXdzdYdz_11.SetMinimum(2)
dXdzdYdz_11.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_11.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_11.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_11.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_11.GetXaxis().SetNdivisions(508)
dXdzdYdz_11.GetYaxis().SetNdivisions(508)
dXdzdYdz_11.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUSref_200_prerad_dXdzdYdz.eps")
c.Print("dataUSref_200_prerad_dXdzdYdz_pq.jpg")

c.Clear()
dXdzdYdz_12 = f1.Get("dataUSref_dXdzdYdz")
# dXdzdYdz_12.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_12.SetStats(0)
dXdzdYdz_12.SetTitle("")
# dXdzdYdz_12.SetMaximum(300)
dXdzdYdz_12.SetMinimum(0)
dXdzdYdz_12.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_12.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_12.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_12.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_12.GetXaxis().SetNdivisions(508)
dXdzdYdz_12.GetYaxis().SetNdivisions(508)
dXdzdYdz_12.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataUSref_200_dXdzdYdz.eps")
c.Print("dataUSref_200_dXdzdYdz_pq.jpg")

# =====================================================


c.Clear()
XY_10 = f1.Get("dataDSref_alltof_XY")
# XY_10.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_10.SetStats(0)
XY_10.SetTitle("")
# XY_10.SetMaximum(400)
XY_10.SetMinimum(10)
XY_10.GetXaxis().SetLabelSize(0.05)
XY_10.GetXaxis().SetTitleSize(0.05)
XY_10.GetYaxis().SetLabelSize(0.05)
XY_10.GetYaxis().SetTitleSize(0.05)
XY_10.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDSref_200_alltof_XY.eps")
c.Print("dataDSref_200_alltof_XY_pq.jpg")

c.Clear()
XY_11 = f1.Get("dataDSref_prerad_XY")
#XY_11.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_11.SetStats(0)
XY_11.SetTitle("")
# XY_11.SetMaximum(400)
XY_11.SetMinimum(2)
XY_11.GetXaxis().SetLabelSize(0.05)
XY_11.GetXaxis().SetTitleSize(0.05)
XY_11.GetYaxis().SetLabelSize(0.05)
XY_11.GetYaxis().SetTitleSize(0.05)
XY_11.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDSref_200_prerad_XY.eps")
c.Print("dataDSref_200_prerad_XY_pq.jpg")

c.Clear()
XY_12 = f1.Get("dataDSref_XY")
# XY_12.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
XY_12.SetStats(0)
XY_12.SetTitle("")
# XY_12.SetMaximum(400)
XY_12.SetMinimum(1)
XY_12.GetXaxis().SetLabelSize(0.05)
XY_12.GetXaxis().SetTitleSize(0.05)
XY_12.GetYaxis().SetLabelSize(0.05)
XY_12.GetYaxis().SetTitleSize(0.05)
XY_12.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDSref_200_XY.eps")
c.Print("dataDSref_200_XY_pq.jpg")

c.Clear()
dXdzdYdz_10 = f1.Get("dataDSref_alltof_dXdzdYdz")
# dXdzdYdz_10.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_10.SetStats(0)
dXdzdYdz_10.SetTitle("")
# dXdzdYdz_10.SetMaximum(300)
dXdzdYdz_10.SetMinimum(10)
dXdzdYdz_10.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_10.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_10.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_10.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_10.GetXaxis().SetNdivisions(508)
dXdzdYdz_10.GetYaxis().SetNdivisions(508)
dXdzdYdz_10.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDSref_200_alltof_dXdzdYdz.eps")
c.Print("dataDSref_200_alltof_dXdzdYdz_pq.jpg")

c.Clear()
dXdzdYdz_11 = f1.Get("dataDSref_prerad_dXdzdYdz")
# dXdzdYdz_11.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_11.SetStats(0)
dXdzdYdz_11.SetTitle("")
# dXdzdYdz_11.SetMaximum(300)
dXdzdYdz_11.SetMinimum(2)
dXdzdYdz_11.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_11.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_11.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_11.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_11.GetXaxis().SetNdivisions(508)
dXdzdYdz_11.GetYaxis().SetNdivisions(508)
dXdzdYdz_11.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDSref_200_prerad_dXdzdYdz.eps")
c.Print("dataDSref_200_prerad_dXdzdYdz_pq.jpg")

c.Clear()
dXdzdYdz_12 = f1.Get("dataDSref_dXdzdYdz")
# dXdzdYdz_12.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
dXdzdYdz_12.SetStats(0)
dXdzdYdz_12.SetTitle("")
# dXdzdYdz_12.SetMaximum(300)
dXdzdYdz_12.SetMinimum(0)
dXdzdYdz_12.GetXaxis().SetLabelSize(0.05)
dXdzdYdz_12.GetXaxis().SetTitleSize(0.05)
dXdzdYdz_12.GetYaxis().SetLabelSize(0.05)
dXdzdYdz_12.GetYaxis().SetTitleSize(0.05)
dXdzdYdz_12.GetXaxis().SetNdivisions(508)
dXdzdYdz_12.GetYaxis().SetNdivisions(508)
dXdzdYdz_12.Draw("colz")
# t1.Draw()
t2.Draw()
c.Print("dataDSref_200_dXdzdYdz.eps")
c.Print("dataDSref_200_dXdzdYdz_pq.jpg")

morehist = ["projposUSDSdiff_GEANT",
            "thetaXUS_thetaXDSGEANT"]

for hist in morehist:
    c.Clear()
    h = f1.Get(hist)
    print hist
    # h.SetTitle("TkU Station 1, LiH 200 MeV/c, March 2016, MAUS v2.5")
    h.SetStats(0)
    h.SetTitle("")
    # h.SetMaximum(300)
    h.SetMinimum(0)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetXaxis().SetNdivisions(508)
    h.GetYaxis().SetNdivisions(508)
    h.Draw("colz")
    # t1.Draw()
    t2.Draw()
    c.Print(hist + "_200.eps")
    c.Print(hist + "_200_pq.jpg")
    


# if __name__ == '__main__':

    
