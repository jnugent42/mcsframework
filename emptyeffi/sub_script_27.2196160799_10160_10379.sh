#!/bin/bash

source /data/neutrino06/jnugent/MAUS-v3.1.2/env.sh
python /data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_27.2196160799.py --not_require_all_planes --track_algorithm 0 --output_filename plots_27.2196160799_10160_10379_EEE /data/neutrino01/jnugent/LiHMC/*ZeroAbs*/{10160..10379}_sim.root

