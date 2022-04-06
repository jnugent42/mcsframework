#!/bin/bash

cp backup/mom_offset.xml .
str2=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
cd MCTOFsys
python run_TOFsys.py
cd ..
str1=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
until [  "$str1" == "$str2" ];
do
  sleep 10m
  echo "sleeping"
  str1=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
  if [  "$str1" == "$str2" ]; then
    break
  fi
done
xmllint -format mom_offset.xml | cat > hold
mv hold mom_offset.xml

cp backup/USDS_mom_offset.xml .
cd USDS_mom_cor
python run_TOFsys.py
cd ..
str1=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
until [  "$str1" == "$str2" ];
do
  sleep 10m
  echo "sleeping"
  str1=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
  if [  "$str1" == "$str2" ]; then
    break
  fi
done
xmllint -format USDS_mom_offset.xml | cat > hold
mv hold USDS_mom_offset.xml

str1=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
cd scanrotang
python run_rot_ang.py
cd ..
str1=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
until [  "$str1" == "$str2" ];
do
  sleep 10m
  echo "sleeping"
  str1=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
  if [  "$str1" == "$str2" ]; then
    break
  fi
done

cd scanrotang
mv rot_ang_0/LiHMu_3200_0.xml .
./run_script.sh
mv LiHMu_3200_0.xml rot_ang_0
cd ..

cd fiducial
python run_MCSbyxml.py
cd ..

cd alignment
python run_alignmentsys.py
cd ..

cd TOFsys
python run_TOFsys.py
cd ..

cd angdef
python run_angdefsys.py
cd ..

cd batch
condor_submit batchscript_hack
cd ..

cd scanradial
python run_accp_radial.py
cd ..

cd batch
./submitsplitjobbashscript.sh
cd ..

: '
cd MCTOFsys
python run_TOFsys.py
cd ..
'

str1=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
until [  "$str1" == "$str2" ];
do
  sleep 10m
  echo "sleeping"
  str1=$(condor_q | grep "Total for jnugent: 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
  if [  "$str1" == "$str2" ]; then
    break
  fi
done

cd MCSsys
cp ../lih172/LihMuon_03172.root .
cp ../lih200/LihMuon_03200.root .
cp ../lih240/LihMuon_03240.root .

cp ../MClih172/LihMuon_03172.root MCLihMuon_03172.root
cp ../MClih200/LihMuon_03200.root MCLihMuon_03200.root
cp ../MClih240/LihMuon_03240.root MCLihMuon_03240.root

python extract_MCSsys.py
cd ..

cd MCDatacom
python MCDatacom.py
cd ..

cd plot3
python threeMomentaPlot.py
cd ..

cd scanradial
mv accp_rad_0/LiHMu_3200_0.xml .
./run_script.sh
mv LiHMu_3200_0.xml accp_rad_0
cd ..

cd TOFsys
mv LiHMu_3200_5.xml xml
cp ../lih172/LihMuon_03172.root .
cp ../lih200/LihMuon_03200.root .
cp ../lih240/LihMuon_03240.root .
./run_script.sh
cp xml/LiHMu_3200_5.xml .
cd ..

cd offsplot
python offsetplot.py
cd ..
