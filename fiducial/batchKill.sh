#!/bin/bash
#first ssh into ppepbs
for i in {141365..141453}
do
	condor_rm ${i}
	#condor_release ${i}
done;
