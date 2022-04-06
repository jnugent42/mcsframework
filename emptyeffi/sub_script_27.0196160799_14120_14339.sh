#!/bin/bash

source /data/neutrino06/jnugent/MAUS-v3.1.2/env.sh
python /data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_27.0196160799.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.0196160799_14120_14339_EEE /data/neutrino01/jnugent/LiHMC/*ZeroAbs*/{14120..14339}_sim.root

