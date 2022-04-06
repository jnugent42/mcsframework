#!/bin/bash 

#python extract_tof_sys.py ../lih172/LihMuon_3172.xml ../lih200/LihMuon_3200.xml ../lih240/LihMuon_3240.xml *.xml
#python extract_rot_ang.py ../lih172/LihMuon_3172.xml ../lih200/LihMuon_3200.xml ../lih240/LihMuon_3240.xml rot_ang_0/LiHMu_3200_*xml
python extract_rot_ang.py ../MClih172/0LihMuon_3172.xml ../MClih200/0LihMuon_3200.xml ../MClih240/0LihMuon_3240.xml rot_ang_0/LiHMu_3200_*xml
#python extract_tof_sys.py ../lih172/LihMuon_3172.xml ../lih200/LihMuon_3200.xml ../lih240/LihMuon_3240.xml LiHMu_3172_tof*xml
