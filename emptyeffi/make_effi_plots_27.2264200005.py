#!/usr/bin/env python
# This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
#
# MAUS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MAUS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MAUS.  If not, see <http://www.gnu.org/licenses/>.
#
# pylint: disable = W0311, E1101, W0613, W0621, C0103, C0111, W0702, W0611
# pylint: disable = R0914, R0912, R0915, W0603, W0612, C0302

# Import MAUS Framework (Required!)
import MAUS

# Generic Python imports
import sys
import os
import argparse
import math
from math import sqrt
import array

# Third Party library import statements
import json
import event_loader
import analysis
from analysis import tools
from analysis import covariances
from analysis import hit_types
import ROOT

def calculate_efficiency() :
    t1 = ROOT.TText(0.28,0.185,"MICE preliminary [simulation]")
    t2 = ROOT.TText(0.28,0.15,"ISIS Cycle 2015/04")
    t1.SetNDC(1)
    t1.SetTextSize(0.04)
    t1.SetTextFont(42)
    t2.SetNDC(1)
    t2.SetTextSize(0.03)
    t2.SetTextFont(42)

    pEffx = ROOT.TEfficiency()
    EEE = "EEE"
    inf = ROOT.TFile("plot_"+str(27.2264200005)+"_"+str(EEE)+".root","READ")
    dir = inf.Get("downstream")
    ds_recon_theta_x = dir.Get("recon_theta_x")
    ds_mc_theta_x = dir.Get("MC_theta_x;1")
    pEffx = ROOT.TEfficiency(ds_recon_theta_x,ds_mc_theta_x)
    inf.Close()
    # f = ROOT.TFile("tracker_resolution_plots_"+str(27.2264200005)+"_"+str(EEE)+".root","RECREATE")
    f = ROOT.TFile("tracker_resolution_plots_"+str(27.2264200005)+".root","RECREATE")
    c5 = ROOT.TCanvas()
    pEffx_graph = pEffx.CreateGraph()
    pEffx_graph.SetName("Effx_graph")
    pEffx_graph.GetXaxis().SetTitle("#theta_{x} (mrad)")
    pEffx_graph.GetYaxis().SetTitle("Efficiency")
    pEffx_graph.Draw("ap")
    f1 = ROOT.TF1("f1","pol2",-0.040,0.040)
    t1.Draw("same")
    t2.Draw("same")
    t1.Paint()
    t2.Paint()
    c5.SaveAs("pEff_x.pdf")
    pEffx_graph.Write()
    f.Close()

    pEffy = ROOT.TEfficiency()
    inf = ROOT.TFile("plot_"+str(27.2264200005)+"_"+str(EEE)+".root","READ")
    dir = inf.Get("downstream")
    ds_recon_theta_y = dir.Get("recon_theta_y")
    ds_mc_theta_y = dir.Get("MC_theta_y")
    pEffy = ROOT.TEfficiency(ds_recon_theta_y,ds_mc_theta_y)
    inf.Close()
    # f = ROOT.TFile("tracker_resolution_plots_"+str(27.2264200005)+"_"+str(EEE)+".root","UPDATE")
    f = ROOT.TFile("tracker_resolution_plots_"+str(27.2264200005)+".root","UPDATE")
    c7 = ROOT.TCanvas()
    pEffy_graph = pEffy.CreateGraph()
    pEffy_graph.SetName("Effy_graph")
    pEffy_graph.GetXaxis().SetTitle("#theta_{y} (mrad)")
    pEffy_graph.GetYaxis().SetTitle("Efficiency")
    pEffy_graph.Draw("ap")
    t1.Draw("same")
    t2.Draw("same")
    t1.Paint()
    t2.Paint()
    c7.SaveAs("pEff_y.pdf")
    pEffy_graph.Write()
    f.Close()

    pEffscatt = ROOT.TEfficiency()
    inf = ROOT.TFile("plot_"+str(27.2264200005)+"_"+str(EEE)+".root","READ")
    dir = inf.Get("downstream")
    ds_recon_theta_scatt = dir.Get("recon_theta_scatt")
    ds_mc_theta_scatt = dir.Get("MC_theta_scatt")
    pEffscatt = ROOT.TEfficiency(ds_recon_theta_scatt,ds_mc_theta_scatt)
    inf.Close()
    # f = ROOT.TFile("tracker_resolution_plots_"+str(27.2264200005)+"_"+str(EEE)+".root","UPDATE")
    f = ROOT.TFile("tracker_resolution_plots_"+str(27.2264200005)+".root","UPDATE")
    c17 = ROOT.TCanvas()
    pEffscatt.Draw()
    c17.SaveAs("pEff_scatt.pdf")
    pEffscatt_graph = pEffscatt.CreateGraph()
    pEffscatt_graph.SetName("Effscatt_graph")
    pEffscatt_graph.Write()
    f.Close()

    pEff2scatt = ROOT.TEfficiency()
    inf = ROOT.TFile("plot_"+str(27.2264200005)+"_"+str(EEE)+".root","READ")
    dir = inf.Get("downstream")
    ds_recon_theta_2scatt = dir.Get("recon_theta_2scatt")
    ds_mc_theta_2scatt = dir.Get("MC_theta_2scatt")
    pEff2scatt = ROOT.TEfficiency(ds_recon_theta_2scatt,ds_mc_theta_2scatt)
    inf.Close()
    # f = ROOT.TFile("tracker_resolution_plots_"+str(27.2264200005)+"_"+str(EEE)+".root","UPDATE")
    f = ROOT.TFile("tracker_resolution_plots_"+str(27.2264200005)+".root","UPDATE")
    c17 = ROOT.TCanvas()
    pEff2scatt.Draw()
    c17.SaveAs("pEff_2scatt.pdf")
    pEff2scatt_graph = pEff2scatt.CreateGraph()
    pEff2scatt_graph.SetName("Eff2scatt_graph")
    pEff2scatt_graph.Write()
    f.Close()

    c1 = ROOT.TCanvas()
    c1.SaveAs('recon_theta.pdf')
    c1.Clear()
    c4= ROOT.TCanvas()
    c2 = ROOT.TCanvas()
    c2.SaveAs('MC_theta.pdf')

if __name__ == "__main__" :
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  calculate_efficiency()

  print
  print "Complete."
  print

