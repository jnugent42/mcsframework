# cd tracker_effi
# python run_trkr_eff.py
# cd ..
# str1=$(condor_q | grep "0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
# str2="0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended"
# until [  "$str1" == "$str2" ];
# do
#   sleep 1h
#   echo "sleeping"
#   if [  "$str1" == "$str2" ]; then
#     break
#   fi
# done

# cd tracker_effi_sel
# python run_trkr_eff.py
# cd ..
# str1=$(condor_q | grep "0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended")
# str2="0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended"
# until [  "$str1" == "$str2" ];
# do
#   sleep 1h
#   echo "sleeping"
#   if [  "$str1" == "$str2" ]; then
#     break
#   fi
#
 done
