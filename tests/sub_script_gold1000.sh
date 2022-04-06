#!/bin/bash

cd /data/neutrino03/jnugent/Unfolding/tests/gold1000
. ../../local_env.sh
pwd
ls
#mkdir /tmp/johndata
date
#cp /data/neutrino02/jnugent/coectedoutput/* /tmp/johndata
date
/data/neutrino03/jnugent/Unfolding/tests/gold1000/MCSUnfolding 0LihMuon_3200.xml
date
cp /tmp/johndata/LihMuon_03200.root /data/neutrino03/jnugent/Unfolding/tests/gold1000/
date
