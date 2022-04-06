#!/bin/bash

for i in {97175..97193}
do
	condor_rm ${i}
done;
