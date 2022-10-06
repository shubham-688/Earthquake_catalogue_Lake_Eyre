#!/bin/csh

# rm catalog/may/*.txt
perl runREAL.pl
cat catalog/jun/phase_sel_*.txt > catalog/phase_all_ps_57_jun.txt

cat catalog/jun/catalog_sel_*.txt | awk '{if ($16<= 200 && $14>1) print $0}' > catalog/catalog_all_ps_57_sg200_jun.txt
# python3 get_phases.py

#.125_sta_0.3_s0_tdx4.5.txt $S = "6/0/9/1/2.5/2/1.6/2.5"
#-R4.5/10/0.1/1000/90/360 -G5.5/10/0.05/10 -S5/0/6/1/2/2/2.25/2 -V6.35/3.7 phase_all_.1_sta_0.3_ps_6_MU_LC_rsel_2
