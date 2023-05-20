#!/bin/sh
#####
# This script makes gifs out of jpg stored in folders inside folders.
# It used 'convert'
# First run this to convert pdf to jpg

for day in {331..331}; do
  cd "$day"
  for index in {1..25}; do
    cd trig_"$index"
    for numbers in {1..26}; do
      cpdir=item_"$numbers".pdf
      outdir=item_"$numbers".jpg
      #converts .pdf to .jpg
      convert -density 100 -quiet $cpdir $outdir
      echo " I did $numbers"
      done
    # takes all jgp files and sort them in ascending order of index
    path=`ls *jpg | awk -F"_" '{print $1,$2,$3,$4}' | awk -F"." '{print $1, $2}' | sort -g -k 2,3 | awk '{print $1"_"$2"."$3}'`
    # creates a gif file out of jpg's
    convert -delay 150 -quiet -quality 100 -loop 3 $path "$day"_"$index".gif
    cp "$day"_"$index".gif ../gifs/
    rm -f *.jpg
    # mkdir lala
    cd ../
  done
  cd ~/Research/Lake_eyre_data/sAus_eq/picker
done

#### copies the gifs from individual trig folder to day gifs folder
# for day in {335..340}; do
#   cd "$day"
#   mkdir gifs
#   for index in {1..20}; do
#     cd trig_"$index"
#     cp "$day"_"$index".gif ../gifs/
#     # takes all jgp files and sort them in ascending order of index
#     # mkdir lala
#     cd ../
#   done
#   cd ~/Research/Lake_eyre_data/sAus_eq/picker
# done

# BUBU's original..
# for numbers in {1..25}; do
#
# cpdir=item_"$numbers".pdf
# outdir=item_"$numbers".jpg
#
# convert -density 300 $cpdir $outdir
# mv $outdir jpgs/
#
# echo " I did $numbers"
# done

#Second, run this command inside the jpg/ folder and obtain the list of jpg.
#echo `ls *jpg | awk -F"_" '{print $1,$2,$3,$4}' | awk -F"." '{print $1, $2}' | sort -g -k 2,3 | awk '{print $1"_"$2"."$3}'`
#COPY THE LIST AND PASTE IT AFTER THIS COMMAND:
# convert -delay 150 -quality 150

#BASICALLY
#convert -delay 150 -quality 150 item_1.jpg item_2.jpg item_3.jpg item_4.jpg item_5.jpg item_6.jpg my_gpoint.gif
