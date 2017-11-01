#!/bin/bash
#first ssh into ppepbs
for i in {12180..12237}
do
	condor_rm ${i}
done;
