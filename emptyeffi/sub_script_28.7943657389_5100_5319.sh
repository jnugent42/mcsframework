#!/bin/bash

source /data/neutrino06/jnugent/MAUS-v3.1.2/env.sh
python /data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_28.7943657389.py --not_require_all_planes --track_algorithm 0 --output_filename plots_28.7943657389_5100_5319_EEE /data/neutrino01/jnugent/LiHMC/*ZeroAbs*/{5100..5319}_sim.root

