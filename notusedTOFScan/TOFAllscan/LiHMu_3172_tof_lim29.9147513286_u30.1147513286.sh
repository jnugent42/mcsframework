#!/bin/bash
cd /data/neutrino03/jnugent/Unfolding/TOFScan
. /data/neutrino03/jnugent/Unfolding/TOFScan/local_env.sh
/data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3172_tof_lim29.9147513286_u30.1147513286.xml
# python CompCobbData.py LiHMu_3172_tof_lim29.9147513286_u30.1147513286MC.root LiHMu_3172_tof_lim29.9147513286_u30.1147513286_mod.root /data/neutrino03/jnugent/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3172_tof_lim29.9147513286_u30.1147513286.xml
