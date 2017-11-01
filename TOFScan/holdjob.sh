#!/bin/bash
#first ssh into ppepbs
for i in {13499..13799}
do
	#condor_hold ${i}
	condor_release ${i}
done;
