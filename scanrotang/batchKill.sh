#!/bin/bash
#first ssh into ppepbs
for i in {58923..59039}
do
	condor_rm ${i}
done;
