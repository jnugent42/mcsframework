#!/bin/bash

pwd

#for bs in actgold generateMCS goldeffbayes99 goldeffnegresp goldeffresamp1 ptruth truthset;
#for bs in bayes10it bayes100it bayes1000it bayes10000it goldeffresamp1 bayesres20 bayesres200 bayesres2000 generateMCS generateMCS100 generateMCS1000 generateMCS10000 gold10 gold100 gold1000 gold10000 truthset;
#for bs in gold10befoacc gold100befoacc gold1000befoacc gold10000befoacc gold10noacc gold100noacc gold1000noacc gold10000noacc gold10 gold100 gold1000 gold10000 truthset;
#for bs in newmodel oldmodel;
for bs in SL6;
do
   sed -e "s/XX/$bs/g" sub_script_XX.sh >& sub_script_$bs.sh
   sed -e "s/XX/$bs/g" batchscript_XX >& batchscript_$bs
   condor_submit batchscript_$bs
done;

