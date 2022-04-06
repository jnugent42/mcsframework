#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_27.8918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.8918072213_12140_12359_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{12140..12359}_sim.root

