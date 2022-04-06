#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_28.0918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_28.0918072213_8620_8839_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{8620..8839}_sim.root

