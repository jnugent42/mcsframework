#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_28.0918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_28.0918072213_15440_15659_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{15440..15659}_sim.root

