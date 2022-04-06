#!/bin/bash

source /data/neutrino06/jnugent/MAUS-v3.1.2/env.sh
python /data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_XX.py --not_require_all_planes --track_algorithm 0 --output_filename plots_XX_LL_uu_EEE /data/neutrino01/jnugent/LiHMC/*ZeroAbs*/{LL..uu}_sim.root

