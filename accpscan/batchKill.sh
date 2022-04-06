#!/bin/bash
#first ssh into ppepbs
for i in {43754..43833}
do
	condor_rm ${i}
done;
