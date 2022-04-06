#!/bin/bash

llim_172=$(bc <<< "29.2147513286-0.1")
ulim_172=$(bc <<< "29.2147513286+0.1")

# llim_200=$(bc <<< "28.0537703996-0.1")
# ulim_200=$(bc <<< "28.0537703996+0.1")
llim_200=$(bc <<< "27.9918072213 - 0.1")
ulim_200=$(bc <<< "27.9918072213 + 0.1")

llim_240=$(bc <<< "27.6398395872-0.1")
ulim_240=$(bc <<< "27.6398395872+0.1")

for i in {-5..10} ;
do

	# file=$(bc <<< "$llim_172+$i*0.2")
	# fileu=$(bc <<< "$ulim_172+$i*0.2")
	# sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
	# sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
	# sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
        # sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
        # condor_submit batchscript_$file

	file=$(bc <<< "$llim_200+$i*0.2")
	fileu=$(bc <<< "$ulim_200+$i*0.2")
	sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
	sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
	sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
    sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
	for i in {4000..16000..220}
	do
		sed -e "s/LL/$i/g" sub_script_${file}.sh >& sub_script_${file}_${i}.sh
                sed -e "s/LL/$i/g" batchscript_${file} >& batchscript_${file}_${i}
		j=$(bc <<< $i+219)
		sed -e "s/uu/$j/g" sub_script_${file}_${i}.sh >& sub_script_${file}_${i}_${j}.sh
                sed -e "s/uu/$j/g" batchscript_${file}_${i} >& batchscript_${file}_${i}_${j}
                condor_submit batchscript_${file}_${i}_${j}
	done

        # condor_submit batchscript_$file

	# file=$(bc <<< "$llim_240+$i*0.2")
	# fileu=$(bc <<< "$ulim_240+$i*0.2")
	# sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
	# sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
	# sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
        # sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
        # condor_submit batchscript_$file

done;

: '
file=$(bc <<< "28.0537703996-0.1")
fileu=$(bc <<< "28.0537703996+0.1")
counter=0
for k in $(seq 0.055 0.01 0.195) ;
do

    counter="$(echo $counter+1 | bc)"
	sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_${fileu}.py
	sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_${fileu}.py >& tracker_resolution_plots_${file}.py
	sed -e "s/EEE/$k/g" tracker_resolution_plots_${file}.py >& tracker_resolution_plots_${file}_${counter}.py
	sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_${file}.sh
    sed -e "s/XX/$file/g" batchscript_XX >& batchscript_${file}
	sed -e "s/EEE/$counter/g" sub_script_${file}.sh >& sub_script_${file}_${counter}.sh
    sed -e "s/EEE/$counter/g" batchscript_${file} >& batchscript_${file}_${counter}
	for i in {4000..16000..110}
	do
		sed -e "s/LL/$i/g" sub_script_${file}_${counter}.sh >& sub_script_${file}_${counter}_${i}.sh
                sed -e "s/LL/$i/g" batchscript_${file}_${counter} >& batchscript_${file}_${counter}_${i}
		j=$(bc <<< $i+109)
		sed -e "s/uu/$j/g" sub_script_${file}_${counter}_${i}.sh >& sub_script_${file}_${counter}_${i}_${j}.sh
                sed -e "s/uu/$j/g" batchscript_${file}_${counter}_${i} >& batchscript_${file}_${counter}_${i}_${j}
                condor_submit batchscript_${file}_${counter}_${i}_${j}
	done

done;
'
: '
llim_172=$(bc <<< "28.850290925-0.1")
ulim_172=$(bc <<< "28.850290925+0.1")

llim_200=$(bc <<< "28.0818843981-0.1")
ulim_200=$(bc <<< "28.0818843981+0.1")

llim_240=$(bc <<< "27.3975634301-0.1")
ulim_240=$(bc <<< "27.3975634301+0.1")

for i in {-10..10} ;
do
	file=$(bc <<< "$llim_172+$i*0.2")
	fileu=$(bc <<< "$ulim_172+$i*0.2")
	sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
	sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
	sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
        sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
        condor_submit batchscript_$file

	file=$(bc <<< "$llim_200+$i*0.2")
	fileu=$(bc <<< "$ulim_200+$i*0.2")
	sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
	sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
	sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
        sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
        condor_submit batchscript_$file

	file=$(bc <<< "$llim_240+$i*0.2")
	fileu=$(bc <<< "$ulim_240+$i*0.2")
	sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
	sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
	sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
        sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
        condor_submit batchscript_$file
done;

file=29.1538507523
fileu=29.3538507523
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
condor_submit batchscript_$file

file=28.449106608
fileu=28.649106608
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
condor_submit batchscript_$file

file=27.4456328413
fileu=27.6456328413
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
condor_submit batchscript_$file

file=$(bc <<< "29.272394075-0.1")
fileu=$(bc <<< "29.272394075+0.1")
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
condor_submit batchscript_$file

file=$(bc <<< "28.4129102539-0.1")
fileu=$(bc <<< "28.4129102539+0.1")
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
condor_submit batchscript_$file

file=$(bc <<< "27.5329625324-0.1")
fileu=$(bc <<< "27.5329625324+0.1")
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
condor_submit batchscript_$file

file=$(bc <<< "27.4975634301-0.1")
fileu=$(bc <<< "27.4975634301+0.1")
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
condor_submit batchscript_$file

file=$(bc <<< "28.850290925-0.1")
fileu=$(bc <<< "28.850290925+0.1")
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
condor_submit batchscript_$file

'
# These two
# file=$(bc <<< "28.7371139354-0.1")
# fileu=$(bc <<< "28.7371139354+0.1")
file=$(bc <<< " 28.7023106126-0.1")
fileu=$(bc <<< " 28.7023106126+0.1")
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
for i in {4000..16000..220}
do
	sed -e "s/LL/$i/g" sub_script_${file}.sh >& sub_script_${file}_${i}.sh
            sed -e "s/LL/$i/g" batchscript_${file} >& batchscript_${file}_${i}
	j=$(bc <<< $i+219)
	sed -e "s/uu/$j/g" sub_script_${file}_${i}.sh >& sub_script_${file}_${i}_${j}.sh
            sed -e "s/uu/$j/g" batchscript_${file}_${i} >& batchscript_${file}_${i}_${j}
    echo batchscript_${file}_${i}_${j}
    condor_submit batchscript_${file}_${i}_${j}
done

# file=$(bc <<< "27.3541567795-0.1")
# fileu=$(bc <<< "27.3541567795+0.1")
file=$(bc <<< "27.2643870826-0.1")
fileu=$(bc <<< "27.2643870826+0.1")
sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
for i in {4000..16000..220}
do
	sed -e "s/LL/$i/g" sub_script_${file}.sh >& sub_script_${file}_${i}.sh
            sed -e "s/LL/$i/g" batchscript_${file} >& batchscript_${file}_${i}
	j=$(bc <<< $i+219)
	sed -e "s/uu/$j/g" sub_script_${file}_${i}.sh >& sub_script_${file}_${i}_${j}.sh
            sed -e "s/uu/$j/g" batchscript_${file}_${i} >& batchscript_${file}_${i}_${j}
    condor_submit batchscript_${file}_${i}_${j}
done

# file=$(bc <<< "28.0537703996-0.1")
# fileu=$(bc <<< "28.0537703996+0.1")
# sed -e "s/TOF_ul = 30/TOF_ul = $fileu/g" tracker_resolution_plots_XX.py >& tracker_resolution_plots_$fileu.py
# sed -e "s/TOF_ll = 27/TOF_ll = $file/g" tracker_resolution_plots_$fileu.py >& tracker_resolution_plots_$file.py
# sed -e "s/XX/$file/g" sub_script_XX.sh >& sub_script_$file.sh
# sed -e "s/XX/$file/g" batchscript_XX >& batchscript_$file
# for i in {4000..16000..220}
# do
# 	sed -e "s/LL/$i/g" sub_script_${file}.sh >& sub_script_${file}_${i}.sh
#             sed -e "s/LL/$i/g" batchscript_${file} >& batchscript_${file}_${i}
# 	j=$(bc <<< $i+219)
# 	sed -e "s/uu/$j/g" sub_script_${file}_${i}.sh >& sub_script_${file}_${i}_${j}.sh
#             sed -e "s/uu/$j/g" batchscript_${file}_${i} >& batchscript_${file}_${i}_${j}
#     condor_submit batchscript_${file}_${i}_${j}
# done

