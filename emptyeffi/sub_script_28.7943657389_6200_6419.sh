#!/bin/bash

source /data/neutrino06/jnugent/MAUS-v3.1.2/env.sh
python /data/neutrino03/jnugent/Unfolding/emptyeffi/tracker_resolution_plots_28.7943657389.py --not_require_all_planes --track_algorithm 0 --output_filename plots_28.7943657389_6200_6419_EEE /data/neutrino01/jnugent/LiHMC/*ZeroAbs*/{6200..6419}_sim.root

