#!/bin/bash
cd /home/ppe/j/jnugent/workarea/Unfolding/TOFScan
. /home/ppe/j/jnugent/workarea/Unfolding/TOFScan/local_env.sh
/home/ppe/j/jnugent/workarea/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3200_tof_lim27.7365596446_u27.9365596446.xml
# python CompCobbData.py LiHMu_3200_tof_lim27.7365596446_u27.9365596446MC.root LiHMu_3200_tof_lim27.7365596446_u27.9365596446_mod.root /home/ppe/j/jnugent/workarea/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3200_tof_lim27.7365596446_u27.9365596446.xml
