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

def updatedatafilename(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("datafile") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def updatetrainfilename(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("trainfile") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def updatemodelfilename(infile, outfile, newname):
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("spec/file"):
        if node.prop("id").find("modelfile") >= 0:
            node.setProp("name", newname)
    f = open(outfile,"w")
    doc.saveTo(f)
    f.close()
    doc.freeDoc()

def submitUnfolding(xmlfile):
    execdir = "/afs/phas.gla.ac.uk/user/j/jnugent/workarea/Unfolding/fiducial"
    cmd = [os.path.join(execdir,"../MCSUnfolding"), xmlfile] 
    proc = subprocess.Popen(cmd)
    proc.wait()

def print_batch_submission(configfile, settings):
#PBS -N %(name)s
#PBS -q long6
#PBS -l walltime=23:59:00,mem=3000Mb
#PBS -e %(working_dir)s/logs/%(name)s.err
#PBS -o %(working_dir)s/logs/%(name)s.out
    batch_data = '''
cd %(working_dir)s
. %(maus_root_dir)s/local_env.sh
%(execcmd)s %(xmlfile)s
'''% settings
    
    outfile = open(configfile, 'w+')
    outfile.write(batch_data)
    outfile.close()

def submit_to_batch(xmlfile):
    execdir = "/afs/phas.gla.ac.uk/user/j/jnugent/workarea/Unfolding/fiducial"
    execcmd = os.path.join(execdir,"../MCSUnfolding")
    maus_root_dir = "/afs/phas.gla.ac.uk/user/j/jnugent/workarea/Unfolding/fiducial"
    working_dir = os.getcwd()
    name = xmlfile[:-4]
    settings = {"xmlfile":xmlfile, "name":name, "working_dir":working_dir, "maus_root_dir":maus_root_dir, "execcmd":execcmd}
    batch_file = os.path.join(working_dir,  name + ".job")
    print_batch_submission(batch_file, settings)
    #cmd = ['qsub', batch_file]
    cmd = [batch_file]
    q = subprocess.Popen(cmd, shell=True)
    q.wait()

for i in range(8):
    for j in range(-20,-13):
        '''
        updatefilename("LiHMu_3172_3.xml", "LiHMu_3172_R_"+str(i)+"_G_"+str(j)+".xml", "LiHMu_3172_R_"+str(i)+"_G_"+str(j)+".root")
        updatecutval("LiHMu_3172_R_"+str(i)+"_G_"+str(j)+".xml", "LiHMu_3172_R_"+str(i)+"_G_"+str(j)+".xml", "fid_rad", 100 + i*10)
        updatecutval("LiHMu_3172_R_"+str(i)+"_G_"+str(j)+".xml", "LiHMu_3172_R_"+str(i)+"_G_"+str(j)+".xml", "fid_grad", (100 + i*10)*0.0003 + 0.002*j)
        submit_to_batch("LiHMu_3172_R_"+str(i)+"_G_"+str(j)+".xml")
        '''
        updatefilename("LiHMu_3200_0.xml", "LiHMu_3200_R_"+str(i)+"_G_"+str(j)+".xml", "LiHMu_3200_R_"+str(i)+"_G_"+str(j)+".root")
        updatecutval("LiHMu_3200_R_"+str(i)+"_G_"+str(j)+".xml", "LiHMu_3200_R_"+str(i)+"_G_"+str(j)+".xml", "fid_rad", 100 + i*10)
        updatecutval("LiHMu_3200_R_"+str(i)+"_G_"+str(j)+".xml", "LiHMu_3200_R_"+str(i)+"_G_"+str(j)+".xml", "fid_grad", (100 + i*10)*0.0003 + 0.002*j)
        submit_to_batch("LiHMu_3200_R_"+str(i)+"_G_"+str(j)+".xml")
        '''
        updatefilename("LiHMu_3240_3.xml", "LiHMu_3240_R_"+str(i)+"_G_"+str(j)+".xml", "LiHMu_3240_R_"+str(i)+"_G_"+str(j)+".root")
        updatecutval("LiHMu_3240_R_"+str(i)+"_G_"+str(j)+".xml", "LiHMu_3240_R_"+str(i)+"_G_"+str(j)+".xml", "fid_rad", 100 + i*10)
        updatecutval("LiHMu_3240_R_"+str(i)+"_G_"+str(j)+".xml", "LiHMu_3240_R_"+str(i)+"_G_"+str(j)+".xml", "fid_grad", 0.002*j)
        submit_to_batch("LiHMu_3240_R_"+str(i)+"_G_"+str(j)+".xml")
        
        updatefilename("XePion_3240.xml", "XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "XePion_3240_R_"+str(i)+"_G_"+str(j)+".root")
        updatecutval("XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "fid_rad", 100 + i*10)
        updatecutval("XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml", "fid_grad", 0.005*j)
        submit_to_batch("XePion_3240_R_"+str(i)+"_G_"+str(j)+".xml")
        '''
# submitUnfolding("LiHMu_3172.xml")
# submitUnfolding("LiHMu_3200.xml")
# submitUnfolding("LiHMu_3240.xml")

'''
rndm = TRandom3()
for i in range(30):
    updatesysval("LiHMu_3172.xml", "LiHMu_3172_"+str(i)+".xml", "alXUS", -0.356 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "alYUS",  1.13 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "alXDS", -0.126 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "alYDS",  0.118 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thXUS", -0.0031 +  rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thYUS", -0.005 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thXDS", -0.003 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml", "thYDS", -0.0151 + rndm.Gaus(0, 0.004))
    updatefilename("LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".xml","LiHMu_3172_"+str(i)+".root")
    submitUnfolding("LiHMu_3172_"+str(i)+".xml")


for i in range(30):
    updatesysval("LiHMu_3200.xml", "LiHMu_3200_"+str(i)+".xml", "alXUS", -0.356 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "alYUS",  1.13 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "alXDS", -0.126 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "alYDS",  0.118 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thXUS", -0.0031 +  rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thYUS", -0.005 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thXDS", -0.003 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml", "thYDS", -0.0151 + rndm.Gaus(0, 0.004))
    updatefilename("LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".xml","LiHMu_3200_"+str(i)+".root")
    submitUnfolding("LiHMu_3200_"+str(i)+".xml")


for i in range(30):
    updatesysval("LiHMu_3240.xml", "LiHMu_3240_"+str(i)+".xml", "alXUS", -0.356 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "alYUS",  1.13 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "alXDS", -0.126 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "alYDS",  0.118 + rndm.Gaus(0, 0.24))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thXUS", -0.0031 +  rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thYUS", -0.005 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thXDS", -0.003 + rndm.Gaus(0, 0.004))
    updatesysval("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml", "thYDS", -0.0151 + rndm.Gaus(0, 0.004))
    updatefilename("LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".xml","LiHMu_3240_"+str(i)+".root")
    submitUnfolding("LiHMu_3240_"+str(i)+".xml")
'''
