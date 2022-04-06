#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_29.8918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_29.8918072213_5320_5539_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{5320..5539}_sim.root

