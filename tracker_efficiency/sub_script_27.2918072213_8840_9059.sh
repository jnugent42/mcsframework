#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_27.2918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.2918072213_8840_9059_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{8840..9059}_sim.root

