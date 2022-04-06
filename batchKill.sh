#!/bin/bash
condor_q | grep -oP '(?<=jnugent ID: )[0-9]+' | xargs condor_rm
#for i in {0..3000}
#do
#	condor_rm ${i}
#	#condor_release ${i}
#done;
