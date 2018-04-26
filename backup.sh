#!/bin/bash

find /data/neutrino03/jnugent/Unfolding -maxdepth 1 -not -name "Unfolding" -not -name "materials" -not -name "TOFScan" -not -name "TOFsys" -not -name "fiducial" -not -name "alignment" -not -name "momentummeas" -exec cp -vr {} . \;
