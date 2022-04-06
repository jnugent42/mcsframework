#!/bin/bash

#python extract_rot_ang.py ../lih172/LihMuon_3172.xml ../lih200/LihMuon_3200.xml ../lih240/LihMuon_3240.xml rot_ang_0/LiHMu_3200_*xml
mv rot_ang_0/LiHMu_3200_0.xml .
python extract_mod_rot_ang.py rot_ang_0/LiHMu_3200_*xml
mv LiHMu_3200_0.xml rot_ang_0
