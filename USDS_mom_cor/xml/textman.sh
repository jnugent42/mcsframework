#!/bin/bash

for ifname in *.xml; do
	#sed -e "s/ZeroMuon\_AllData\/reduced\_tree\_datat1\.root/\/data\/neutrino02\/jnugent\/data\_reduced\/allempty\/reduced\_tree\_Abs\_03172\.root/g" $ifname >& holder
	sed -e "s/TOFAllScan/TOFAllscan/" $ifname >& holder
	mv holder $ifname
done
