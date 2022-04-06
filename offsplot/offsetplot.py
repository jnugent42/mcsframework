import libxml2
import os, subprocess
from ROOT import TRandom3, TMath
import matplotlib.pyplot as plt
from ROOT import TROOT

TROOT.gROOT.SetBatch(1)

def getval(infile, name):
    x = [[],[]]
    doc = libxml2.parseFile(infile)
    for node in doc.xpathEval("offset/"+name):
        x[0].append(node.prop("id"))
        x[1].append(node.prop("value"))
    doc.freeDoc()
    return x

# x = [[],[]]

xmlname1 = "/data/neutrino03/jnugent/Unfolding/mom_offset.xml"
xmlname2 = "/data/neutrino03/jnugent/Unfolding/USDS_mom_offset.xml"

x = getval(xmlname1,"LiH")
xref = getval(xmlname1,"ref")
y = getval(xmlname2,"LiH")
yref = getval(xmlname2,"ref")

fig, axs = plt.subplots(2,2,figsize=(15,15))
fig.suptitle('Offsets')
axs[0,0].plot(x[0],x[1],'ro')
axs[0,0].set_ylabel('offset (MeV/c)')
axs[0,0].set_title('LiH DS')
axs[0,0].set_xlabel('TOF centre window (ns)')

axs[0,1].plot(xref[0],xref[1],'ro')
axs[0,1].set_title('ref DS')
axs[0,1].set_ylabel('offset (MeV/c)')
axs[0,1].set_xlabel('TOF centre window (ns)')

axs[1,0].plot(y[0],y[1],'ro')
axs[1,0].set_title('LiH US')
axs[1,0].set_ylabel('offset (MeV/c)')
axs[1,0].set_xlabel('TOF centre window (ns)')

axs[1,1].plot(yref[0],yref[1],'ro')
axs[1,1].set_title('ref US')
axs[1,1].set_ylabel('offset (MeV/c)')
axs[1,1].set_xlabel('TOF centre window (ns)')

fig.savefig("offsets.pdf")

