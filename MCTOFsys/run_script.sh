#!/bin/bash

mv LiHMu_3200_5.xml xml
cp ../lih172/LihMuon_03172.root .
cp ../lih200/LihMuon_03200.root .
cp ../lih240/LihMuon_03240.root .
python extract_tof_sys.py ../lih172/LihMuon_3172.xml ../lih200/LihMuon_3200.xml ../lih240/LihMuon_3240.xml *.xml
cp xml/LiHMu_3200_5.xml .
