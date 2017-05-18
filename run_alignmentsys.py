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
    execdir = "/data/neutrino02/rbayes/MICE/MCS_selection/Unfolding/v2.5.0"
    cmd = [os.path.join(execdir,"MCSUnfolding"), xmlfile] 
    proc = subprocess.Popen(cmd)
    return proc

def print_batch_submission(configfile, settings):
    batch_data = '''
#PBS -N %(name)s
#PBS -q medium6
#PBS -l walltime=5:59:00,mem=3000Mb
#PBS -e %(working_dir)s/logs/%(name)s.err
#PBS -o %(working_dir)s/logs/%(name)s.out

cd %(working_dir)s
. %(maus_root_dir)s/local_env.sh
%(execcmd)s %(xmlfile)s
'''% settings
    
    outfile = open(configfile, 'w+')
    outfile.write(batch_data)
    outfile.close()

def submit_to_batch(xmlfile):
    execdir = "/data/neutrino02/rbayes/MICE/MCS_selection/Unfolding/v2.5.0"
    execcmd = os.path.join(execdir,"MCSUnfolding")
    maus_root_dir = "/data/neutrino02/rbayes/MICE/MCS_selection/Unfolding/v2.5.0"
    working_dir = os.getcwd()
    name = xmlfile[:-4]
    settings = {"xmlfile":xmlfile, "name":name, "working_dir":working_dir, "maus_root_dir":maus_root_dir, "execcmd":execcmd}
    batch_file = os.path.join(working_dir,  name + ".job")
    print_batch_submission(batch_file, settings)
    cmd = ['qsub', batch_file]
    q = subprocess.Popen(cmd)
    q.wait()


llim_172 = 29.272394075 - 0.1
ulim_172 = 29.272394075 + 0.1

llim_200 = 28.4129102539 - 0.1
ulim_200 = 28.4129102539 + 0.1

llim_240 = 27.5329625324 - 0.1
ulim_240 = 27.5329625324 + 0.1
'''
updatefilename("LiHMu_3172MC_1.xml", "LiHMu_3172MC_4.xml", "LiHMuon_03172MC_4.root") 
updatecutval("LiHMu_3172MC_4.xml", "LiHMu_3172MC_4.xml", "TOF_ll", llim_172)
updatecutval("LiHMu_3172MC_4.xml", "LiHMu_3172MC_4.xml", "TOF_ul", ulim_172)

updatefilename("LiHMu_3200MC_1.xml", "LiHMu_3200MC_4.xml", "LiHMuon_03200MC_4.root") 
updatecutval("LiHMu_3200MC_4.xml", "LiHMu_3200MC_4.xml", "TOF_ll", llim_200)
updatecutval("LiHMu_3200MC_4.xml", "LiHMu_3200MC_4.xml", "TOF_ul", ulim_200)

updatefilename("LiHMu_3240MC_1.xml", "LiHMu_3240MC_4.xml", "LiHMuon_03240MC_4.root") 
updatecutval("LiHMu_3240MC_4.xml", "LiHMu_3240MC_4.xml", "TOF_ll", llim_240)
updatecutval("LiHMu_3240MC_4.xml", "LiHMu_3240MC_4.xml", "TOF_ul", ulim_240)

#proc1 = submitUnfolding("LiHMu_3172MC_4.xml")
#proc2 = submitUnfolding("LiHMu_3200MC_4.xml")
#proc3 = submitUnfolding("LiHMu_3240MC_4.xml")

#proc1.wait() 
#CompCobbData.dataModelComps("LiHMuon_03172MC_4.root", "../../Cobb_results/LiH-6p5-172-B-hists.dat", "../../Cobb_results/LiH-6p5-172-B-hists.root")
#proc2.wait()
#CompCobbData.dataModelComps("LiHMuon_03200MC_4.root", "../../Cobb_results/LiH-6p5-200-B-hists.dat", "../../Cobb_results/LiH-6p5-200-B-hists.root")
#proc3.wait()
#CompCobbData400.dataModelComps("LiHMuon_03240MC_4.root", "../../Cobb_results/LiH-6p5-240-hists.dat", "../../Cobb_results/LiH-6p5-240-hists.root")

updatefilename("LiHMu_3172_2.xml", "LiHMu_3172_9.xml", "LiHMuon_03172_9.root") 
updatecutval("LiHMu_3172_9.xml", "LiHMu_3172_9.xml", "TOF_ll", llim_172)
updatecutval("LiHMu_3172_9.xml", "LiHMu_3172_9.xml", "TOF_ul", ulim_172)
updatesysval("LiHMu_3172_9.xml", "LiHMu_3172_9.xml", "niter", 10)

updatefilename("LiHMu_3200_2.xml", "LiHMu_3200_9.xml", "LiHMuon_03200_9.root") 
updatecutval("LiHMu_3200_9.xml", "LiHMu_3200_9.xml", "TOF_ll", llim_200)
updatecutval("LiHMu_3200_9.xml", "LiHMu_3200_9.xml", "TOF_ul", ulim_200)
updatesysval("LiHMu_3200_9.xml", "LiHMu_3200_9.xml", "niter", 10)

updatefilename("LiHMu_3240_2.xml", "LiHMu_3240_9.xml", "LiHMuon_03240_9.root") 
updatecutval("LiHMu_3240_9.xml", "LiHMu_3240_9.xml", "TOF_ll", llim_240)
updatecutval("LiHMu_3240_9.xml", "LiHMu_3240_9.xml", "TOF_ul", ulim_240)
updatesysval("LiHMu_3240_9.xml", "LiHMu_3240_9.xml", "niter", 10)

submit_to_batch("LiHMu_3172_9.xml")
submit_to_batch("LiHMu_3200_9.xml")
submit_to_batch("LiHMu_3240_9.xml")
submit_to_batch("XePion_3240.xml")

proc1.wait()
proc2.wait()
proc3.wait()
''' 
# submitUnfolding("LiHMu_3172_pcut.xml")
# submitUnfolding("LiHMu_3200_pcut.xml")
# submitUnfolding("LiHMu_3240_pcut.xml")


rndm = TRandom3()

# for i in range(100):
    # updatecutval("LiHMu_3172.xml", "LiHMu_3172_"+str(i)+".xml", "TOF_ll", llim_172)
    # updatecutval("LiHMu_3172_"+str(i)+".xml", "LiHMu_3172_"+str(i)+".xml", "TOF_ul", ulim_172)
    # updatesysval("LiHMu_3172_"+str(i)+".xml", "LiHMu_3172_"+str(i)+".xml", "alXUS",0.356 + rndm.Gaus(0, 0.24))
    # updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "alYUS", -1.13 + rndm.Gaus(0, 0.24))
    # updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "alXDS", 0.126 + rndm.Gaus(0, 0.24))
    # updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "alYDS", -0.118 + rndm.Gaus(0, 0.24))
    # updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thYUS", 0.0031 +  rndm.Gaus(0, 0.004))
    # updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thXUS", 0.005 + rndm.Gaus(0, 0.004))
    # updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thYDS", 0.003 + rndm.Gaus(0, 0.004))
    # updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thXDS", 0.0151 + rndm.Gaus(0, 0.004))
    # updatefilename("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".root")
    # submit_to_batch("LiHMu_3172_"+str(i)+".xml")


for i in range(500):
    updatecutval("LiHMu_3200_6.xml", "LiHMu_3200_"+str(i)+".xml", "TOF_ll", llim_200)
    # updatecutval("LiHMu_3200_"+str(i)+".xml", "LiHMu_3200_"+str(i)+".xml", "TOF_ul", ulim_200)
    # updatesysval("LiHMu_3200_"+str(i)+".xml", "LiHMu_3200_"+str(i)+".xml", "alXUS", 0.356 + rndm.Gaus(0, 0.24))
    # updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "alYUS",  -1.13 + rndm.Gaus(0, 0.24))
    # updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "alXDS",  0.126 + rndm.Gaus(0, 0.24))
    # updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "alYDS",  -0.118 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thYUS",  -0.0031 +  0.01*rndm.Rndm() - 0.005)
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thXUS",  -0.005 + 0.01*rndm.Rndm() - 0.005)
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thYDS",  -0.003 + 0.01*rndm.Rndm()- 
0.005)
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thXDS",  -0.0151 + 0.03*rndm.Rndm() - 0.015)
    updatefilename("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".root")
    submit_to_batch("LiHMu_3200_"+str(i)+".xml")

'''
for i in range(100):
    updatecutval("LiHMu_3240.xml", "LiHMu_3240_"+str(i)+".xml", "TOF_ll", llim_240)
    updatecutval("LiHMu_3240_"+str(i)+".xml", "LiHMu_3240_"+str(i)+".xml", "TOF_ul", ulim_240)
    updatesysval("LiHMu_3240_"+str(i)+".xml", "LiHMu_3240_"+str(i)+".xml", "alXUS", 0.356 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "alYUS",  -1.13 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "alXDS",  0.126 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "alYDS",  -0.118 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thYUS", 0.0031 +  rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thXUS", 0.005 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thYDS", 0.003 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thXDS", 0.0151 + rndm.Gaus(0, 0.004))
    updatefilename("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".root")
    submit_to_batch("LiHMu_3240_"+str(i)+".xml")


for i in range(100):
    updatecutval("XePion_3240.xml", "XePion_3240_"+str(i)+".xml", "TOF_ll", llim_240)
    updatecutval("XePion_3240_"+str(i)+".xml", "XePion_3240_"+str(i)+".xml", "TOF_ul", ulim_240)
    updatesysval("XePion_3240_"+str(i)+".xml", "XePion_3240_"+str(i)+".xml", "alXUS", 0.15 + rndm.Gaus(0, 0.24))
    updatesysval("XePion_3240_"+str(i)+".xml","XePion_3240_"+str(i)+".xml", "alYUS",  0.066 + rndm.Gaus(0, 0.24))
    updatesysval("XePion_3240_"+str(i)+".xml","XePion_3240_"+str(i)+".xml", "alXDS",  0.099 + rndm.Gaus(0, 0.24))
    updatesysval("XePion_3240_"+str(i)+".xml","XePion_3240_"+str(i)+".xml", "alYDS",  -0.103 + rndm.Gaus(0, 0.24))
    updatesysval("XePion_3240_"+str(i)+".xml","XePion_3240_"+str(i)+".xml", "thYUS",  -0.0019 +  rndm.Gaus(0, 0.004))
    updatesysval("XePion_3240_"+str(i)+".xml","XePion_3240_"+str(i)+".xml", "thXUS",  -0.005 + rndm.Gaus(0, 0.004))
    updatesysval("XePion_3240_"+str(i)+".xml","XePion_3240_"+str(i)+".xml", "thYDS",  -0.000 + rndm.Gaus(0, 0.004))
    updatesysval("XePion_3240_"+str(i)+".xml","XePion_3240_"+str(i)+".xml", "thXDS",  -0.00391 + rndm.Gaus(0, 0.004))
    updatefilename("XePion_3240_"+str(i)+".xml","XePion_3240_"+str(i)+".xml","XePion_3240_"+str(i)+".root")
    submit_to_batch("XePion_3240_"+str(i)+".xml")
'''

#  LocalWords:  ulim
