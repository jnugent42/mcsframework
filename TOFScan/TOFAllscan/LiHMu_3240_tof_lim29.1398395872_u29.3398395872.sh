#!/bin/bash
cd /home/ppe/j/jnugent/workarea/Unfolding/TOFScan
. /home/ppe/j/jnugent/workarea/Unfolding/TOFScan/local_env.sh
/home/ppe/j/jnugent/workarea/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3240_tof_lim29.1398395872_u29.3398395872.xml
# python CompCobbData400.py LiHMu_3240_tof_lim29.1398395872_u29.3398395872MC.root LiHMu_3240_tof_lim29.1398395872_u29.3398395872_mod.root /home/ppe/j/jnugent/workarea/Unfolding/TOFScan/../MCSUnfolding TOFAllscan/LiHMu_3240_tof_lim29.1398395872_u29.3398395872.xml
