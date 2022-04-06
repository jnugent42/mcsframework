#!/bin/bash

source /data/neutrino/jnugent/MAUS-v3.3.2/env.sh
python /home/ppe/j/jnugent/workarea/tracker_efficiency/tracker_resolution_plots_27.1643870826.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.1643870826_13240_uu_EEE /data/neutrino/jnugent/LiHMC/*LiH*/{13240..uu}_sim.root

