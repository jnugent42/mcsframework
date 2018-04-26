#!/bin/bash

cd MClih172
../MCSUnfolding LihMuon_3172.xml
../MCSUnfolding nopionLihMuon_3172.xml
../MCSUnfolding acceptLihMuon_3172.xml
cd ..

cd MClih200
../MCSUnfolding LihMuon_3200.xml
../MCSUnfolding nopionLihMuon_3200.xml
../MCSUnfolding acceptLihMuon_3200.xml
cd ..

cd MClih240
../MCSUnfolding LihMuon_3240.xml
../MCSUnfolding nopionLihMuon_3240.xml
../MCSUnfolding acceptLihMuon_3240.xml
cd ..

cd lih172
../MCSUnfolding LihMuon_3172.xml
cd ..

cd lih200
../MCSUnfolding LihMuon_3200.xml
cd ..

cd lih240
../MCSUnfolding LihMuon_3240.xml
cd ..

