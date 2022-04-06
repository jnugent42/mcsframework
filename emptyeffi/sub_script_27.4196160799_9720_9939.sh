#!/bin/bash

source /data/neutrino06/jnugent/MAUS-v3.1.2/env.sh
python /data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_27.4196160799.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.4196160799_9720_9939_EEE /data/neutrino01/jnugent/LiHMC/*ZeroAbs*/{9720..9939}_sim.root

