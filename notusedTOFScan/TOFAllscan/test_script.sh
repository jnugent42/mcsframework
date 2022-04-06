#!/bin/bash
cd /home/ppe/j/jnugent/workarea/Unfolding/TOFScan
. /home/ppe/j/jnugent/workarea/Unfolding/TOFScan/local_env.sh
/home/ppe/j/jnugent/workarea/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3172_tof_lim28.9147513286_u29.1147513286MC.xml
# python CompCobbData.py TOFAllscan/LiHMu_3172_tof_lim28.9147513286_u29.1147513286MC.root TOFAllscan/LiHMu_3172_tof_lim28.9147513286_u29.1147513286_mod.root /home/ppe/j/jnugent/workarea/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3172_tof_lim28.9147513286_u29.1147513286.xml
