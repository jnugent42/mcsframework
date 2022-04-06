#!/bin/bash

for i in {96812..96861}
do
	#condor_rm -forcex ${i}
	condor_rm ${i}
done;
