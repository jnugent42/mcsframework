#!/bin/bash

llim_172=$(bc <<< "29.2147513286-0.1")
ulim_172=$(bc <<< "29.2147513286+0.1")

llim_200=$(bc <<< "27.9918072213-0.1")
ulim_200=$(bc <<< "27.9918072213+0.1")
# llim_200=$(bc <<< "31.0472422527 - 0.1")
# ulim_200=$(bc <<< "31.0472422527 + 0.1")

llim_240=$(bc <<< "27.6398395872-0.1")
ulim_240=$(bc <<< "27.6398395872+0.1")

# counter=0
# for k in $(seq 0.055 0.01 0.195) ;
# do

#     counter="$(echo $counter+1 | bc)"
#     ls plots_${llim_200}_*_${counter}.root | xargs hadd plot_${llim_200}_${counter}.root

# done;

for i in {-5..10} ;
do

	file=$(bc <<< "$llim_200+$i*0.2")
	fileu=$(bc <<< "$ulim_200+$i*0.2")
    ls plots_${file}_*_EEE.root | xargs hadd plot_${file}_EEE.root

done;

file=$(bc <<< "28.7023106126-0.1")
fileu=$(bc <<< "28.7023106126+0.1")
ls plots_${file}_* | xargs hadd plot_${file}_EEE.root

file=$(bc <<< "27.2643870826-0.1")
fileu=$(bc <<< "27.2643870826+0.1")
ls plots_${file}_* | xargs hadd plot_${file}_EEE.root
