#!/bin/bash
cd /data/neutrino03/jnugent/Unfolding/TOFScan
. /data/neutrino03/jnugent/Unfolding/TOFScan/local_env.sh
/data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3200_tof_lim29.5365596446_u29.7365596446.xml
# python CompCobbData.py LiHMu_3200_tof_lim29.5365596446_u29.7365596446MC.root LiHMu_3200_tof_lim29.5365596446_u29.7365596446_mod.root /data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3200_tof_lim29.5365596446_u29.7365596446.xml
