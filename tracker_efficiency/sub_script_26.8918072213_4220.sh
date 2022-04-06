#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_26.8918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_26.8918072213_4220_uu_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{4220..uu}_sim.root

