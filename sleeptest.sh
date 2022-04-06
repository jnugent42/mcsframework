#!/bin/bash

str1=$(condor_q | grep "0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
str2="0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended"
until [  "$str1" == "$str2" ];
do
  sleep 1s
  echo "sleeping"
  if [  "$str1" == "$str2" ]; then
    echo "end"
    break
  fi
done
