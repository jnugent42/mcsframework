#!/bin/bash

find *.xml -not -name LiHMu_3172.xml -not -name LiHMu_3200_6.xml -not -name LiHMu_3240.xml -exec rm {} \;
find *.sh -not -name cleaup.sh -not -name local_env.sh -exec rm {} \;
rm *.job
rm *.out
rm *.err
rm *.log
rm *.root
