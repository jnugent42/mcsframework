#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_27.6918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.6918072213_13020_13239_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{13020..13239}_sim.root

