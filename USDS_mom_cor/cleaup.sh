#!/bin/bash

find *.xml -not -name LiHMu_3172_0.xml -not -name LiHMu_3200_5.xml -not -name LiHMu_3240_0.xml -exec rm {} \;
find *.sh -not -name cleaup.sh -not -name local_env.sh -not -name run_script.sh -exec rm {} \;
rm *.job
rm *.out
rm *.err
rm *.log
rm *.root
