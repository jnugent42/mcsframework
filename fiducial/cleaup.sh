#!/bin/bash

find *.xml -not -name LiHMu_3172_0.xml -not -name LiHMu_3200_0.xml -not -name LiHMu_3240_0.xml -exec rm {} \;
rm *.job
rm *.out
rm *.err
rm *.log
rm *.root
