#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_28.6023106126.py --not_require_all_planes --track_algorithm 0 --output_filename plots_28.6023106126_6200_uu_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{6200..uu}_sim.root

