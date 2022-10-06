#!/bin/csh
#collates picks from each day into a single monthly txt file.
perl runREAL.pl
cat catalog/jan/phase_sel_*.txt > catalog/phase_all_ps_79_jan.txt
cat catalog/jan/catalog_sel_*.txt | awk '{if ($16<= 200 && $14>1) print $0}' > catalog/catalog_all_ps_79_sg200_jan.txt
