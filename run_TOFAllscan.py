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
    print xmlfile
    execdir = "/data/neutrino02/rbayes/MICE/MCS_selection/Unfolding"
    cmd = [os.path.join(execdir,"MCSUnfolding"), xmlfile] 
    proc = subprocess.Popen(cmd)
    # proc.wait()
    return proc
    
def print_batch_submission(configfile, settings):
    batch_data = '''
#PBS -N %(name)s
#PBS -q long6
#PBS -l walltime=23:59:00,mem=3000Mb
#PBS -e %(working_dir)s/logs/%(name)s.err
#PBS -o %(working_dir)s/logs/%(name)s.out

cd %(working_dir)s
. %(maus_root_dir)s/local_env.sh
# %(execcmd)s %(xmlfileMC)s
# python %(modelcmd)s %(mcfile)s %(ccsource) %(modelfile)
%(execcmd)s %(xmlfile)s
'''% settings
    
    outfile = open(configfile, 'w+')
    outfile.write(batch_data)
    outfile.close()

def print_MCbatch_submission(configfile, settings):
    batch_data = '''
#PBS -N %(name)s
#PBS -q long6
#PBS -l walltime=23:59:00,mem=3000Mb
#PBS -e %(working_dir)s/logs/%(name)s.err
#PBS -o %(working_dir)s/logs/%(name)s.out

cd %(working_dir)s
. %(maus_root_dir)s/local_env.sh
%(execcmd)s %(xmlfileMC)s
python %(modelcmd)s %(mcfile)s %(ccsource) %(modelfile)
%(execcmd)s %(xmlfile)s
'''% settings
    
    outfile = open(configfile, 'w+')
    outfile.write(batch_data)
    outfile.close()

def submit_to_batch(xmlfile, runMC=False):
    execdir = "/data/neutrino02/rbayes/MICE/MCS_selection/Unfolding"
    execcmd = os.path.join(execdir,"MCSUnfolding")
    maus_root_dir = "/data/neutrino02/rbayes/MICE/MCS_selection/Unfolding/v2.6.1"
    working_dir = os.getcwd()
    name = xmlfile[:-4]
    if xmlfile.find("3240") >= 0:
        modelcmd = "CompCobbData400.py"
        ccsource = os.path.join(execdir, "Cobb_results/LiH-6p5-240-hists.dat")
    elif xmlfile.find("3200") >= 0:
        modelcmd = "CompCobbData.py"
        ccsource = os.path.join(execdir, "Cobb_results/LiH-6p5-200-B-hists.dat")
    elif xmlfile.find("3172") >= 0:
        modelcmd = "CompCobbData.py"
        ccsource = os.path.join(execdir, "Cobb_results/LiH-6p5-172-B-hists.dat")
    else:
        print "Unrecognized case"
        exit(1)
    settings = {"xmlfile":xmlfile, "xmlfileMC":name+"MC.xml", "mcfile":name+"MC.root", \
                    "name":name, "working_dir":working_dir, "maus_root_dir":maus_root_dir, \
                    "execcmd":execcmd, "modelcmd":modelcmd, "modelfile":name+"_mod.root", \
                    "ccsource":ccsource}
    batch_file = os.path.join(working_dir,  name + ".job")
    if runMC: print_MCbatch_submission(batch_file, settings)
    else: print_batch_submission(batch_file, settings)
    cmd = ['qsub', batch_file]
    q = subprocess.Popen(cmd)
    q.wait()


llim_172 = 29.2147513286 - 0.1
ulim_172 = 29.2147513286 + 0.1

llim_200 = 28.4365596446 - 0.1
ulim_200 = 28.4365596446 + 0.1

llim_240 = 27.6398395872 - 0.1
ulim_240 = 27.6398395872 + 0.1

llim_p240 = 27.8327860284 - 0.3
ulim_p240 = 27.8327860284 + 0.3


for i in range(-10,10):

    print llim_172, ulim_172
    updatefilename("LiHMu_3172_2.xml",
                   "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                   "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".root")
    updatemodelfilename( \
        "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
        "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
        "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"_mod.root")
    updatedatafilename(\
        "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
        "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml", 
        "LiHMuon_AllData/reduced_tree_datat1.root")
    updatetrainfilename(\
        "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
        "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
        "ZeroMuon_AllData/reduced_tree_datat1.root")


    updatecutval("TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                 "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                 "TOF_ll", llim_172 + i*0.2)
    updatecutval("TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                 "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml",
                 "TOF_ul", ulim_172 + i*0.2)

    updatefilename("LiHMu_3172MC_2.xml",
                   "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.xml",
                   "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.root")
    updatedatafilename(\
        "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.xml",
        "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.xml", 
        "LiHMuon_AllData/reduced_tree_datat1.root")

    updatecutval("TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.xml",
                 "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.xml",
                 "TOF_ll", llim_172 + i*0.2)
    updatecutval("TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.xml",
                 "TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.xml",
                 "TOF_ul", ulim_172 + i*0.2)
    submit_to_batch("TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml")
        
    print llim_200, ulim_200
    updatefilename("LiHMu_3200_2.xml",
                   "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                   "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".root")
    updatedatafilename(\
        "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
        "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml", 
        "LiHMuon_AllData/reduced_tree_datat1.root")
    updatetrainfilename(\
        "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
        "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
        "ZeroMuon_AllData/reduced_tree_datat1.root")
    updatemodelfilename( \
                   "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                   "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                   "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"_mod.root")
    updatecutval("TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "TOF_ll", llim_200 + i*0.2)
    updatecutval("TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml",
                 "TOF_ul", ulim_200 + i*0.2)

    updatefilename("LiHMu_3200MC_2.xml",
                   "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.xml",
                   "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.root")
    updatedatafilename(\
        "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.xml",
        "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.xml", 
        "LiHMuon_AllData/reduced_tree_datat1.root")
    updatecutval("TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.xml",
                 "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.xml",
                 "TOF_ll", llim_200 + i*0.2)
    updatecutval("TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.xml",
                 "TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.xml",
                 "TOF_ul", ulim_200 + i*0.2)
    submit_to_batch("TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml")
    
    print llim_240 + i*0.2, ulim_240 + i*0.2
    updatefilename("LiHMu_3240_2.xml",
                   "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                   "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".root")
    updatedatafilename(\
        "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
        "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml", 
        "LiHMuon_AllData/reduced_tree_datat1.root")
    updatetrainfilename(\
        "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
        "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
        "ZeroMuon_AllData/reduced_tree_datat1.root")
    updatemodelfilename( \
                   "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                   "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                   "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"_mod.root")

    updatecutval("TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                 "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                 "TOF_ll", llim_240 + i*0.2)
    updatecutval("TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                 "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml",
                 "TOF_ul", ulim_240 + i*0.2)

    updatefilename("LiHMu_3240MC_2.xml",
                   "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.xml",
                   "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.root")
    updatedatafilename(\
        "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.xml",
        "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.xml", 
        "LiHMuon_AllData/reduced_tree_datat1.root")
    updatecutval("TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.xml",
                 "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.xml",
                 "TOF_ll", llim_240 + i*0.2)
    updatecutval("TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.xml",
                 "TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.xml",
                 "TOF_ul", ulim_240 + i*0.2)
    submit_to_batch("TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml")
    
    
    # proc11 = submitUnfolding("TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml")
    # proc12 = submitUnfolding("TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml")
    # proc13 = submitUnfolding("TOFAllscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml")

    # proc11.wait()
    # proc12.wait()
    # proc13.wait()
    '''
    
    proc11 = submitUnfolding("TOFAllscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.xml")
    proc12 = submitUnfolding("TOFAllscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.xml")
    proc13 = submitUnfolding("TOFscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.xml")
    try:
        proc11.wait()
        CompCobbData.dataModelComps("TOFscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.root",
                                    "Cobb_results/LiH-6p5-172-B-hists.dat", 
                                    "TOFscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"_mod.root")
        proc21 = submitUnfolding("TOFscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml")
        proc21.wait()
        # submit_to_batch("TOFscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml")
    except (AttributeError, ZeroDivisionError):
        print "Cannot filter model information please resubmit TOFscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+"MC.root"
        # submit_to_batch("TOFscan/LiHMu_3172_tof_lim"+str(llim_172 + i*0.2)+"_u"+str(ulim_172 + i*0.2)+".xml", True)


    
    try:
        proc12.wait()
        CompCobbData.dataModelComps("TOFscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.root",
                                    "Cobb_results/LiH-6p5-200-B-hists.dat", 
                                    "TOFscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"_mod.root")
        proc22 = submitUnfolding("TOFscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml")
        # 
        proc22.wait()
    except (AttributeError, ZeroDivisionError):
        print "Cannot filter model information please resubmit TOFscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+"MC.root"
        
        # submit_to_batch("TOFscan/LiHMu_3200_tof_lim"+str(llim_200 + i*0.2)+"_u"+str(ulim_200 + i*0.2)+".xml", True)

    
        
    try:
        proc13.wait()
        CompCobbData400.dataModelComps("TOFscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.root",
                                       "Cobb_results/LiH-6p5-240-hists.dat", 
                                       "TOFscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"_mod.root")
        proc23 = submitUnfolding("TOFscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml")
        # submit_to_batch("TOFscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+".xml")
        proc23.wait()
    except (AttributeError, ZeroDivisionError):
        print "Cannot filter model information please resubmit TOFscan/LiHMu_3240_tof_lim"+str(llim_240 + i*0.2)+"_u"+str(ulim_240 + i*0.2)+"MC.root"
    
        '''
