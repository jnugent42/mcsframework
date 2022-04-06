#!/bin/bash
#first ssh into ppepbs
for i in {432508..432540}
do
	condor_rm ${i}
	#condor_hold ${i}
done;
