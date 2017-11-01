#!/bin/bash
#first ssh into ppepbs
for i in {2741..2760}
do
	condor_rm ${i}
done;
