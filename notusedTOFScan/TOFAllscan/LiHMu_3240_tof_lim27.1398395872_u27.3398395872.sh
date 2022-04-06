#!/bin/bash
cd /data/neutrino03/jnugent/Unfolding/TOFScan
. /data/neutrino03/jnugent/Unfolding/TOFScan/local_env.sh
/data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3240_tof_lim27.1398395872_u27.3398395872.xml
# python CompCobbData400.py LiHMu_3240_tof_lim27.1398395872_u27.3398395872MC.root LiHMu_3240_tof_lim27.1398395872_u27.3398395872_mod.root /data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3240_tof_lim27.1398395872_u27.3398395872.xml
