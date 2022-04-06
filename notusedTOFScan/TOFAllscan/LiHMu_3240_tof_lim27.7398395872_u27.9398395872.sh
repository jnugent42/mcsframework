#!/bin/bash
cd /data/neutrino03/jnugent/Unfolding/TOFScan
. /data/neutrino03/jnugent/Unfolding/TOFScan/local_env.sh
/data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3240_tof_lim27.7398395872_u27.9398395872.xml
# python CompCobbData400.py LiHMu_3240_tof_lim27.7398395872_u27.9398395872MC.root LiHMu_3240_tof_lim27.7398395872_u27.9398395872_mod.root /data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3240_tof_lim27.7398395872_u27.9398395872.xml
