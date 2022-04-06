#!/bin/bash
#first ssh into ppepbs
for i in {774342..774359}
do
	condor_rm ${i}
	# condor_release ${i}
done;
