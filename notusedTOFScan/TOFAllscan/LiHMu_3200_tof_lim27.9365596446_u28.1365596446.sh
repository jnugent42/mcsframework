#!/bin/bash
cd /data/neutrino03/jnugent/Unfolding/TOFScan
. /data/neutrino03/jnugent/Unfolding/TOFScan/local_env.sh
/data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3200_tof_lim27.9365596446_u28.1365596446.xml
# python CompCobbData.py LiHMu_3200_tof_lim27.9365596446_u28.1365596446MC.root LiHMu_3200_tof_lim27.9365596446_u28.1365596446_mod.root /data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3200_tof_lim27.9365596446_u28.1365596446.xml
