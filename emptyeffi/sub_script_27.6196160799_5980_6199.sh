#!/bin/bash

source /data/neutrino06/jnugent/MAUS-v3.1.2/env.sh
python /data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_27.6196160799.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.6196160799_5980_6199_EEE /data/neutrino01/jnugent/LiHMC/*ZeroAbs*/{5980..6199}_sim.root

