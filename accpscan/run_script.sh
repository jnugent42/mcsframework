#!/bin/bash 

#python extract_tof_sys.py ../lih172/LihMuon_3172.xml ../lih200/LihMuon_3200.xml ../lih240/LihMuon_3240.xml *.xml
python extract_accpt.py ../lih172/LihMuon_3172.xml ../lih200/LihMuon_3200.xml ../lih240/LihMuon_3240.xml rot_ang_0/LiHMu_3200_*xml
#python extract_tof_sys.py ../lih172/LihMuon_3172.xml ../lih200/LihMuon_3200.xml ../lih240/LihMuon_3240.xml LiHMu_3172_tof*xml
