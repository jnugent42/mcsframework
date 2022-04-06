#!/usr/bin/env python

import libxml2
import os, subprocess
from ROOT import TRandom3, TMath
import CompCobbData400, CompCobbData


def updatesysval(infile, outfile, name, value):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/sys"):
        if node.prop("name").find(name) >= 0:
            node.setProp("value", str(value))
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def updatecutval(infile, outfile, name, value):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/cuts"):
        if node.prop("name").find(name) >= 0:
            node.setProp("value", str(value))
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()


def updatefilename(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("outfile") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def submitUnfolding(xmlfile):
    execdir = "/afs/phas.gla.ac.uk/user/j/jnugent/workarea"
    cmd = [os.path.join(execdir,"Unfolding"), xmlfile] 
    proc = subprocess.Popen(cmd)
    return proc

def run_anal():
    '''
TOF for 172 MeV/c is  29.2147513286
TOF for 200 MeV/c is  28.4365596446
TOF for 240 MeV/c is  27.6398395872
    '''
    llim_172 = 29.2147513286 - 0.1
    ulim_172 = 29.2147513286 + 0.1
    
    llim_200 = 28.4365596446 - 0.1
    ulim_200 = 28.4365596446 + 0.1
    
    llim_240 = 27.6398395872 - 0.1
    ulim_240 = 27.6398395872 + 0.1
    
    llim_p240 = 27.8327860284 - 0.3
    ulim_p240 = 27.8327860284 + 0.3
    
    print llim_172, ulim_172
    print llim_200, ulim_200
    print llim_240, ulim_240
    
    updatefilename("LiHMu_3172MC_2.xml", "LiHMu_3172MC_0.xml", "LiHMuon_03172MC_0.root") 
    updatecutval("LiHMu_3172MC_0.xml", "LiHMu_3172MC_0.xml", "TOF_ll", llim_172)
    updatecutval("LiHMu_3172MC_0.xml", "LiHMu_3172MC_0.xml", "TOF_ul", ulim_172)
    
    updatefilename("LiHMu_3200MC_2.xml", "LiHMu_3200MC_0.xml", "LiHMuon_03200MC_0.root") 
    updatecutval("LiHMu_3200MC_0.xml", "LiHMu_3200MC_0.xml", "TOF_ll", llim_200)
    updatecutval("LiHMu_3200MC_0.xml", "LiHMu_3200MC_0.xml", "TOF_ul", ulim_200)
    
    updatefilename("LiHMu_3240MC_2.xml", "LiHMu_3240MC_0.xml", "LiHMuon_03240MC_0.root") 
    updatecutval("LiHMu_3240MC_0.xml", "LiHMu_3240MC_0.xml", "TOF_ll", llim_240)
    updatecutval("LiHMu_3240MC_0.xml", "LiHMu_3240MC_0.xml", "TOF_ul", ulim_240)
    
    updatefilename("XePion_3240MC.xml", "XePion_3240MC_0.xml", "XePion_3240_0.root") 
    updatecutval("XePion_3240MC_0.xml", "XePion_3240MC_0.xml", "TOF_ll", llim_p240)
    updatecutval("XePion_3240MC_0.xml", "XePion_3240MC_0.xml", "TOF_ul", ulim_p240)

    proc1 = submitUnfolding("LiHMu_3172MC_0.xml")
    proc2 = submitUnfolding("LiHMu_3200MC_0.xml")
    proc3 = submitUnfolding("LiHMu_3240MC_0.xml")
    proc4 = submitUnfolding("XePion_3240MC_0.xml")

    proc1.wait() 
    CompCobbData.dataModelComps("LiHMuon_03172MC_0.root", "../Cobb_results/LiH-6p5-172-B-hists.dat", "../Cobb_results/LiH-6p5-172-B-hists.root")
    proc2.wait()
    CompCobbData.dataModelComps("LiHMuon_03200MC_0.root", "../Cobb_results/LiH-6p5-200-B-hists.dat", "../Cobb_results/LiH-6p5-200-B-hists.root")
    proc3.wait()
    CompCobbData400.dataModelComps("LiHMuon_03240MC_0.root", "../Cobb_results/LiH-6p5-240-hists.dat", "../Cobb_results/LiH-6p5-240-hists.root")
    
    updatefilename("LiHMu_3172_2.xml", "LiHMu_3172_0.xml", "LiHMuon_03172_0.root") 
    updatecutval("LiHMu_3172_0.xml", "LiHMu_3172_0.xml", "TOF_ll", llim_172)
    updatecutval("LiHMu_3172_0.xml", "LiHMu_3172_0.xml", "TOF_ul", ulim_172)
    updatesysval("LiHMu_3172_0.xml", "LiHMu_3172_0.xml", "niter", 10)
    
    updatefilename("LiHMu_3200_2.xml", "LiHMu_3200_0.xml", "LiHMuon_03200_0.root") 
    updatecutval("LiHMu_3200_0.xml", "LiHMu_3200_0.xml", "TOF_ll", llim_200)
    updatecutval("LiHMu_3200_0.xml", "LiHMu_3200_0.xml", "TOF_ul", ulim_200)
    updatesysval("LiHMu_3200_0.xml", "LiHMu_3200_0.xml", "niter", 10)

    updatefilename("LiHMu_3240_2.xml", "LiHMu_3240_0.xml", "LiHMuon_03240_0.root") 
    updatecutval("LiHMu_3240_0.xml", "LiHMu_3240_0.xml", "TOF_ll", llim_240)
    updatecutval("LiHMu_3240_0.xml", "LiHMu_3240_0.xml", "TOF_ul", ulim_240)
    updatesysval("LiHMu_3240_0.xml", "LiHMu_3240_0.xml", "niter", 10)

    updatefilename("XePion_3240_0.xml", "XePion_3240_0.xml", "XePion_3240_0.root") 
    updatecutval("XePion_3240_0.xml", "XePion_3240_0.xml", "TOF_ll", llim_p240)
    updatecutval("XePion_3240_0.xml", "XePion_3240_0.xml", "TOF_ul", ulim_p240)
    updatesysval("XePion_3240_0.xml", "XePion_3240_0.xml", "niter", 10)
    
    proc1 = submitUnfolding("LiHMu_3172_0.xml")
    proc2 = submitUnfolding("LiHMu_3200_0.xml")
    proc3 = submitUnfolding("LiHMu_3240_0.xml")
    proc4 = submitUnfolding("XePion_3240_0.xml")
    
    proc1.wait()
    proc2.wait()
    proc3.wait()
    proc4.wait()
    # submitUnfolding("LiHMu_3172_pcut.xml")
    # submitUnfolding("LiHMu_3200_pcut.xml")
    # submitUnfolding("LiHMu_3240_pcut.xml")
    
if __name__ == '__main__':
    print "running analysis"
    run_anal()
