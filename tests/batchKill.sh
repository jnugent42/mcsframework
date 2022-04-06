#!/bin/bash
#first ssh into ppepbs
for i in {66031..66043}
do
	condor_rm ${i}
	#condor_release ${i}
done;
