#!/bin/bash
for (( i=1; i<=89; i++ ))
do
 cp aust* eq$i/
 cp 5g.sta eq$i/
 cp na.in eq$i/
 cp shna.cmd eq$i/
 cd eq$i/
 shake_na
 cd ../
done
#python3 get_relocated_compare.py
#python3 get_relocated.py
#scp *.txt shubham@es09156:/Users/Shubham/Dropbox/NA_temp/relocated_eqs/
# /Users/matthias/shub/NAPR/Progs
