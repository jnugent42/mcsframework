#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_28.6023106126.py --not_require_all_planes --track_algorithm 0 --output_filename plots_28.6023106126_8400_8619_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{8400..8619}_sim.root

