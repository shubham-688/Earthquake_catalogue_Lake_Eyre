#!/bin/sh

# First run this to convert pdf to jpg
mkdir jpgs

for numbers in {1..25}; do

cpdir=item_"$numbers".pdf
outdir=item_"$numbers".jpg

convert -density 300 $cpdir $outdir
mv $outdir jpgs/

echo " I did $numbers"
done

#Second, run this command inside the jpg/ folder and obtain the list of jpg.
#echo `ls *jpg | awk -F"" '{print $1,$2,$3,$4}' | awk -F"." '{print $1, $2}' | sort -g -k 2,3 | awk '{print $1""$2"_"$3"."$4}'
#COPY THE LIST AND PASTE IT AFTER THIS COMMAND:
# convert -delay 150 -quality 150

#BASICALLY
#convert -delay 150 -quality 150 item_1.jpg item_2.jpg item_3.jpg item_4.jpg item_5.jpg item_6.jpg my_gpoint.gif
