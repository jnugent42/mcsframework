#!/bin/bash

llim_172=$(bc <<< "29.2147513286-0.1")
ulim_172=$(bc <<< "29.2147513286+0.1")

llim_200=$(bc <<< "28.1196160799-0.1")
ulim_200=$(bc <<< "28.1196160799+0.1")
# llim_200=$(bc <<< "31.0472422527 - 0.1")
# ulim_200=$(bc <<< "31.0472422527 + 0.1")

llim_240=$(bc <<< "27.6398395872-0.1")
ulim_240=$(bc <<< "27.6398395872+0.1")

for i in {-5..10} ;
do

	file=$(bc <<< "$llim_200+$i*0.2")
	fileu=$(bc <<< "$ulim_200+$i*0.2")
	sed -e "s/TOF_ll/$file/g" make_effi_plots_XX.py >& make_effi_plots_${file}.py
    python make_effi_plots_${file}.py

done;

file=$(bc <<< "28.8943657389-0.1")
sed -e "s/TOF_ll/$file/g" make_effi_plots_XX.py >& make_effi_plots_${file}.py
python make_effi_plots_${file}.py

file=$(bc <<< "27.3264200005-0.1")
sed -e "s/TOF_ll/$file/g" make_effi_plots_XX.py >& make_effi_plots_${file}.py
python make_effi_plots_${file}.py

# counter=0
# for k in $(seq 0.055 0.01 0.195) ;
# do

#     counter="$(echo $counter+1 | bc)"
# 	sed -e "s/TOF_ll/$llim_200/g" make_effi_plots_XX.py >& make_effi_plots_${llim_200}.py
# 	sed -e "s/EEE/$counter/g" make_effi_plots_${llim_200}.py >& make_effi_plots_${llim_200}_${counter}.py
#     python make_effi_plots_${llim_200}_${counter}.py

# done;
