#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_29.6918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_29.6918072213_4880_5099_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{4880..5099}_sim.root

