import libxml2
import os, subprocess
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

def updatetrkreffiname(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("trkreffiname") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def updateemptytrkreffiname(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("trkreffiemptyname") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def submitUnfolding(xmlfile):
    execdir = "/data/neutrino03/jnugent/Unfolding"
    cmd = [os.path.join(execdir,"MCSUnfolding"), xmlfile]
    proc = subprocess.Popen(cmd)
    return proc

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
log            = %(working_dir)s/%(name)s.log
requirements   = OpSysAndVer == "CentOS7"
request_memory = 6 GB

queue
'''% settings

    outfile = open(configfile, 'w+')
    outfile.write(batch_data)
    outfile.close()

def submit_to_batch(xmlfile):
    execdir = "/data/neutrino03/jnugent/Unfolding/USDS_mom_cor"
    execcmd = os.path.join(execdir,"../MCSUnfolding")
    maus_root_dir = "/data/neutrino03/jnugent/Unfolding/USDS_mom_cor"
    working_dir = os.getcwd()
    name = xmlfile[:-4]
    settings = {"shellscript":xmlfile[:-4]+".sh", "xmlfile":xmlfile, "name":name, "working_dir":working_dir, "maus_root_dir":maus_root_dir, "execcmd":execcmd}
    batch_file = os.path.join(working_dir,  name + ".job")
    print_batch_submission(batch_file, settings)
    cmd = ['condor_submit', batch_file]
    q = subprocess.Popen(cmd)
    q.wait()

'''
llim_172 = 28.850290925 - 0.1
ulim_172 = 28.850290925 + 0.1

llim_200 = 28.1818843981 - 0.1
ulim_200 = 28.1818843981 + 0.1

llim_240 = 27.4975634301 - 0.1
ulim_240 = 27.4975634301 + 0.1

llpi_240 = 27.2
ulpi_240 = 27.8
'''

# llim_200 = 31.0472422527 - 0.1
# ulim_200 = 31.0472422527 + 0.1
llim_200=27.9918072213 - 0.1
ulim_200=27.9918072213 + 0.1
refllim_200 = 28.1196160799 - 0.1
refulim_200 = 28.1196160799 + 0.1

for i in range(-5,7):
    print llim_200, ulim_200
    updatefilename("LiHMu_3200_5.xml",
                   "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                   "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".root")
    updatetrkreffiname( \
	            "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
	            "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
	            "/data/neutrino03/jnugent/Unfolding/tracker_efficiency/tracker_resolution_plots_"+str(llim_200 + i*0.2)+".root")
    updateemptytrkreffiname( \
	            "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
	            "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
	            "/data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_"+str(refllim_200 + i*0.2)+".root")
    updatecutval("LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "TOF_ll", llim_200 + i*0.2)
    updatecutval("LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "TOF_ul", ulim_200 + i*0.2)
    updatecutval("LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "TOF_ll_ref", refllim_200 + i*0.2)
    updatecutval("LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "TOF_ul_ref", refulim_200 + i*0.2)
    submit_to_batch("LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml")

l=0
refcut = [28.7943657389,27.2264200005]
for i in [28.6023106126,27.1643870826]:
    updatefilename("LiHMu_3200_5.xml",
                   "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
                   "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".root")
    updatetrkreffiname( \
	            "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
	            "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
	            "/data/neutrino03/jnugent/Unfolding/tracker_efficiency/tracker_resolution_plots_"+str(i)+".root")
    updateemptytrkreffiname( \
	            "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
	            "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
	            "/data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_"+str(refllim_200 + i)+".root")
    updatecutval("LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
                 "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
                 "TOF_ll", i)
    updatecutval("LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
                 "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
                 "TOF_ul", i+0.2)
    updatecutval("LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
                 "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
                 "TOF_ll_ref", refcut[l])
    updatecutval("LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
                 "LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml",
                 "TOF_ul_ref", refcut[l] + 0.2)
    submit_to_batch("LiHMu_3200_tof_lim"+str(i)+"_u"+str(i+0.2)+".xml")
