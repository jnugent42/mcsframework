#!/bin/bash

for ifname in *.xml; do
	#sed -e "s/ZeroMuon\_AllData\/reduced\_tree\_datat1\.root/\/data\/neutrino02\/jnugent\/data\_reduced\/allempty\/reduced\_tree\_Abs\_03172\.root/g" $ifname >& holder
	#sed -e "s/\/afs\/phas\.gla\.ac\.uk\/user\//\/home\/ppe\//" $ifname >& holder
	sed -e "s/<file id=\"model2\" name=\"Cobb\"\/>/<file id=\"model3\" name=\"Moliere\"\/>\n<file id=\"model2\" name=\"Cobb\"\/>/" $ifname >& holder
	mv holder $ifname
	sed -e "s/<cuts name=\"mode\" value=\"2\"\/>/<cuts name=\"mode\" value=\"-1\"\/>/" $ifname >& holder
	mv holder $ifname
done
