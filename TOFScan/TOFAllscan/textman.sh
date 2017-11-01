#!/bin/bash

for ifname in *.xml; do
	#sed -e "s/ZeroMuon\_AllData\/reduced\_tree\_datat1\.root/\/data\/neutrino02\/jnugent\/data\_reduced\/allempty\/reduced\_tree\_Abs\_03172\.root/g" $ifname >& holder
	sed -e "s/LiHMuon\_AllData\/reduced\_tree\_datat1\.root/\/data\/neutrino06\/jnugent\/v2\.9\.1\/allLiH\/reduced\_tree\_\.root/g" $ifname >& holder
	mv holder $ifname
	sed -e "s/ZeroMuon\_AllData\/reduced\_tree\_datat1\.root/\/data\/neutrino06\/jnugent\/v2\.9\.1\/allZeroAbs\/reduced\_tree\_\.root/g" $ifname >& holder
	mv holder $ifname
done
