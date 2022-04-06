#!/bin/bash

source /data/neutrino06/jnugent/MAUS-v3.1.2/env.sh
python /data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_29.6196160799.py --not_require_all_planes --track_algorithm 0 --output_filename plots_29.6196160799_13900_14119_EEE /data/neutrino01/jnugent/LiHMC/*ZeroAbs*/{13900..14119}_sim.root

