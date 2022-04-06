import libxml2
import os, subprocess
import numpy as np
from ROOT import TRandom3, TMath


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
    execdir = "/data/neutrino03/jnugent/Unfolding/scanoffset"
    cmd = [os.path.join(execdir,"MCSUnfolding"), xmlfile]
    proc = subprocess.Popen(cmd)
    proc.wait()

def print_batch_submission(configfile, settings):

    testscript = '''#!/bin/bash
cd %(working_dir)s
. %(maus_root_dir)s/local_env.sh
%(execcmd)s %(xmlfile)s
'''% settings

    outfilename = '''%(shellscript)s'''% settings
    outfile0 = open(outfilename, 'w+')
    outfile0.write(testscript)
    outfile0.close()

    batch_data = '''
universe       = vanilla
executable     = %(shellscript)s
output         = %(working_dir)s/%(name)s.out
error          = %(working_dir)s/%(name)s.err
log            = test.log
requirements   = OpSysAndVer == "CentOS7"

queue
'''% settings

    outfile = open(configfile, 'w+')
    outfile.write(batch_data)
    outfile.close()

def submit_to_batch(xmlfile):
    execdir = "/data/neutrino03/jnugent/Unfolding/scanoffset"
    execcmd = os.path.join(execdir,"../MCSUnfolding")
    maus_root_dir = "/data/neutrino03/jnugent/Unfolding/scanoffset"
    working_dir = os.getcwd()
    name = xmlfile[:-4]
    settings = {"xmlfile":xmlfile, "name":name, "working_dir":working_dir, "maus_root_dir":maus_root_dir, "execcmd":execcmd, "shellscript":xmlfile[:-4]+".sh"}
    batch_file = os.path.join(working_dir,  name + ".job")
    print_batch_submission(batch_file, settings)
    #cmd = ['qsub', batch_file]
    cmd = ['condor_submit', batch_file]
    q = subprocess.Popen(cmd)
    q.wait()


llim_172 = 1.
ulim_172 = 28.9

llim_200 = 1.
ulim_200 = 27.9

llim_240 = 1.
ulim_240 = 27.3

for i in np.arange(-10,10,1):
    '''
    updatefilename("offset/LiHMu_3172_0.xml",
                   "offset/LiHMu_3172_"+str(i)+".xml",
                   "offset/LiHMu_3172_"+str(i)+".root")
    updatesysval("offset/LiHMu_3172_"+str(i)+".xml",
                 "offset/LiHMu_3172_"+str(i)+".xml",
                 "rot_ang", i)
    submit_to_batch("offset/LiHMu_3172_"+str(i)+".xml")
    '''
    updatefilename("offset/LiHMu_3200_0.xml",
                   "offset/LiHMu_3200_"+str(i)+".xml",
                   "offset/LiHMu_3200_"+str(i)+".root")
    updatecutval("offset/LiHMu_3200_"+str(i)+".xml",
                 "offset/LiHMu_3200_"+str(i)+".xml",
                 "offset", i)
    submit_to_batch("offset/LiHMu_3200_"+str(i)+".xml")
    '''
    updatefilename("offset/LiHMu_3240_0.xml",
                   "offset/LiHMu_3240_"+str(i)+".xml",
                   "offset/LiHMu_3240_"+str(i)+".root")
    updatesysval("offset/LiHMu_3240_"+str(i)+".xml",
                 "offset/LiHMu_3240_"+str(i)+".xml",
                 "rot_ang", i)
    submit_to_batch("offset/LiHMu_3240_"+str(i)+".xml")
    '''
    '''
    updatefilename("offset/XePion_3240.xml",
                   "offset/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                   "offset/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".root")
    updatesysval("offset/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "offset/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "resX", llim_240 + i*0.0109)
    updatesysval("offset/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "offset/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml",
                 "resY", llim_240 + i*0.0109)
    submit_to_batch("offset/XePion_3240_res_lim"+str(llim_240 + i*0.0109)+".xml")

'''
