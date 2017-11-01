#!/bin/bash
cd /home/ppe/j/jnugent/workarea/Unfolding/TOFScan
. /home/ppe/j/jnugent/workarea/Unfolding/TOFScan/local_env.sh
/home/ppe/j/jnugent/workarea/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3172_tof_lim29.1147513286_u29.3147513286.xml
# python CompCobbData.py LiHMu_3172_tof_lim29.1147513286_u29.3147513286MC.root LiHMu_3172_tof_lim29.1147513286_u29.3147513286_mod.root /home/ppe/j/jnugent/workarea/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3172_tof_lim29.1147513286_u29.3147513286.xml
