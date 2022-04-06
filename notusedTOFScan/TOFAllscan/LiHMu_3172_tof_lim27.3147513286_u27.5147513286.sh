#!/bin/bash
cd /data/neutrino03/jnugent/Unfolding/TOFScan
. /data/neutrino03/jnugent/Unfolding/TOFScan/local_env.sh
/data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3172_tof_lim27.3147513286_u27.5147513286.xml
# python CompCobbData.py LiHMu_3172_tof_lim27.3147513286_u27.5147513286MC.root LiHMu_3172_tof_lim27.3147513286_u27.5147513286_mod.root /data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3172_tof_lim27.3147513286_u27.5147513286.xml
