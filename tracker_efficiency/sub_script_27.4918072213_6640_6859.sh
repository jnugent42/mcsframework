#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_27.4918072213.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.4918072213_6640_6859_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{6640..6859}_sim.root

