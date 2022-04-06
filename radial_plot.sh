#!/bin/bash

cd scanradial
mv accp_rad_0/LiHMu_3200_0.xml .
./run_script.sh
mv LiHMu_3200_0.xml accp_rad_0
cd ..
