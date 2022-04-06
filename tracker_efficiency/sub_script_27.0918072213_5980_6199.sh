#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_27.0918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.0918072213_5980_6199_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{5980..6199}_sim.root

