#!/bin/bash
cd /data/neutrino03/jnugent/Unfolding/TOFScan
. /data/neutrino03/jnugent/Unfolding/TOFScan/local_env.sh
/data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3240_tof_lim28.3398395872_u28.5398395872.xml
# python CompCobbData400.py LiHMu_3240_tof_lim28.3398395872_u28.5398395872MC.root LiHMu_3240_tof_lim28.3398395872_u28.5398395872_mod.root /data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3240_tof_lim28.3398395872_u28.5398395872.xml
