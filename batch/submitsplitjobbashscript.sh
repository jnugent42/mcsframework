#!/bin/bash

pwd
for xml in nopion 0;
do
	sed -e "s/YY/$xml/g" mcsub_script_YY.sh >& mcsub_script_$xml.sh
       sed -e "s/YY/$xml/g" mcbatchscript_YY >& mcbatchscript_$xml
        for bs in 172 200 240;
        #for bs in 200;
	do
          echo ${xml}_${bs}
	   sed -e "s/XX/${bs}/g" mcsub_script_${xml}.sh >& mcsub_script_${xml}_${bs}.sh
           sed -e "s/XX/${bs}/g" mcbatchscript_${xml} >& mcbatchscript_${xml}_${bs}
           condor_submit mcbatchscript_${xml}_${bs}
        done;
done;

for bs in 172 200 240 rot200;
do
    sed -e "s/XX/$bs/g" sub_script_XX.sh >& sub_script_$bs.sh
    sed -e "s/XX/$bs/g" batchscript_XX >& batchscript_$bs
    condor_submit batchscript_$bs
done;

