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

def updatefilename(infile, outfile, newname, name):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find(name) >= 0:
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
request_memory = 3 GB

queue
'''% settings

    outfile = open(configfile, 'w+')
    outfile.write(batch_data)
    outfile.close()

def submit_to_batch(xmlfile):
    execdir = "/data/neutrino03/jnugent/Unfolding/tracker_effi_sel"
    execcmd = os.path.join(execdir,"../MCSUnfolding")
    maus_root_dir = "/data/neutrino03/jnugent/Unfolding/tracker_effi_sel"
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
llim_200 = 28.0537703996 - 0.1
ulim_200 = 28.0537703996 + 0.1
refllim_200 = 28.1920848764 - 0.1
refulim_200 = 28.1920848764 + 0.1

for i in range(-5,8):
    '''
    print llim_172, ulim_172
    updatefilename("LiHMu_3172_5.xml",
                   "LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                   "LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".root")
    updatetrkreffiname( \
	            "LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
	            "LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
	            "/home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_"+str(llim_172 + i*0.2)+".root")
    updatecutval("LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                 "LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                 "TOF_ll", llim_172 + i*0.2)
    updatecutval("LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                 "LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                 "TOF_ul", ulim_172 + i*0.2)
    submit_to_batch("LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml")
    # proc0 = submitUnfolding("LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml")
     '''
    print llim_200, ulim_200
    updatefilename("LiHMu_3200_5.xml",
                   "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                   "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".root",
                   "outfile")
    updatetrkreffiname( \
	            "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
	            "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
	            "/home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_"+str(llim_200 + i*0.2)+".root")
    updateemptytrkreffiname( \
	            "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
	            "LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
	            "/home/ppe/j/jnugent/workarea/emptyeffi/tracker_resolution_plots_"+str(refllim_200 + i*0.2)+".root")
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
    for j in range(0,13):
        updatefilename("lihmu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                    "lihmu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                    "/data/neutrino03/jnugent/Unfolding/tracker_effi/output*/",
                    "datafile")
        submit_to_batch("LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml")
    # proc1 = submitUnfolding("LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml")
    '''
    print llim_240, ulim_240
    updatefilename("LiHMu_3240_5.xml",
                   "LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                   "LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".root")
    updatetrkreffiname( \
	            "LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
	            "LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
	            "/home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_"+str(llim_240 + i*0.2)+".root")
    updatecutval("LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                 "LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                 "TOF_ll", llim_240 + i*0.2)
    updatecutval("LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                 "LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                 "TOF_ul", ulim_240 + i*0.2)
    submit_to_batch("LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml")
    # proc2 = submitUnfolding("LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml")
    '''
    '''
    proc2.wait()
    proc1.wait()
    proc0.wait()


    print llpi_240, ulpi_240
    updatefilename("XePion_3240.xml",
                   "XePion_3240_tof_lim"+str(llpi_240 + i*0.2)+"_u"+str(ulpi_240 + i*0.2)+".xml",
                   "XePion_3240_tof_lim"+str(llpi_240 + i*0.2)+"_u"+str(ulpi_240 + i*0.2)+".root")
    updatecutval("XePion_3240_tof_lim"+str(llpi_240 + i*0.2)+"_u"+str(ulpi_240 + i*0.2)+".xml",
                 "XePion_3240_tof_lim"+str(llpi_240 + i*0.2)+"_u"+str(ulpi_240 + i*0.2)+".xml",
                 "TOF_ll", llpi_240 + i*0.2)
    updatecutval("XePion_3240_tof_lim"+str(llpi_240 + i*0.2)+"_u"+str(ulpi_240 + i*0.2)+".xml",
                 "XePion_3240_tof_lim"+str(llpi_240 + i*0.2)+"_u"+str(ulpi_240 + i*0.2)+".xml",
                 "TOF_ul", ulpi_240 + i*0.2)
    submit_to_batch("XePion_3240_tof_lim"+str(llpi_240 + i*0.2)+"_u"+str(ulpi_240 + i*0.2)+".xml")
    '''
