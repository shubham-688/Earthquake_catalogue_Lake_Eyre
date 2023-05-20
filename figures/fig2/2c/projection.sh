#!/bin/csh
#ogr2ogr -f "GMT" aus5wgd_l.gmt  aus5wgd_l.shp -t_srs WGS84
set PS=map_ml.ps
#set Rmap=-Rg
#set Rmap=-R-135/40/-90/90
set Rmap=-R133/140/-30.75/-25.5 # Lake_eyre+ Flinders
#set J=-JM8/22
set J=-JM12/11 #sa
#set J=-Jx0.023i/0.029i

set topogrd = /Users/shubham/Research/etopo/etopo_sa_15s.grd
#set topogrd = ~/Research/etopo/etopo1.grd
# set topoxyz = ~/Research/etopo/etopo1.xyz.

# set oleron = ~/Documents/ScientificColourMaps5/oleron/oleron.cpt
set tpushuf = ~/Documents/cpt-city/tp/tpushuf
set afrikakarte = ~/Documents/cpt-city/wkp/lilleskut/afrikakarte

#set eq_real_all_18 = ~/Research/Lake_eyre_data/sAus_eq/picker/Real/REAL/catalog/final_cat_good/cat_sg200_2018_19.txt
#set eq_real_all_19 = ~/Research/Lake_eyre_data/sAus_eq/picker/Real/REAL/catalog/final_cat_good/cat_sg200_2019.txt
#set eq_real_all_20 = ~/Research/Lake_eyre_data/sAus_eq/picker/Real/REAL/catalog/final_cat_good/cat_sg200_2020.txt

set eq_all_LE = ~/Research/earthquakes/ga_database/LE_allmag_alltime_FR.txt
set eq_NA = ~/Dropbox/NA_temp/relocated_eqs/round_1/relocated_eq_na_aust.txt
set eq_NA_ml = ~/Dropbox/NA_temp/mag_esti/round_1/relocated_eq_na_aust_mlv.txt

set eq_NA_2 = ~/Dropbox/NA_temp/relocated_eqs/relocated_eq_na.txt
set eq_NA_ml_2 = ~/Dropbox/NA_temp/mag_esti/relocated_eq_na_mlv_2.txt

set mine_NA = ~/Dropbox/NA_temp/relocated_eqs/round_1/relocated_mine_na_aust.txt
set mine_NA_ml = ~/Dropbox/NA_temp/mag_esti/round_1/relocated_mine_na_aust_mlv.txt

set mine_NA_2 = ~/Dropbox/NA_temp/relocated_eqs/relocated_mine_na_aust.txt
set mine_NA_ml_2 = ~/Dropbox/NA_temp/mag_esti/relocated_mine_na_aust_mlv_2.txt

set eq_na_ml_F = relocated_FINAL_mlsorted.txt
set mine_na_ml_F = relocated_mine_na_ml_FINAL.txt

#sort -k4 -n -r relocated_eq_na_mlv_FINAL.txt > relocated_FINAL_mlsorted.txt

gmt6 gmtset MAP_FRAME_TYPE plain MAP_FRAME_WIDTH 3p FONT_ANNOT_PRIMARY 9.5p MAP_TICK_LENGTH_PRIMARY 5p MAP_FRAME_PEN 2.2p


gmt6 psbasemap -BneWS -Bxa1f.5 -Bya1f.5 $J $Rmap -K >! $PS
#gmt6 grdimage $topogrd $J $Rmap -CETOPO1.cpt -K -O -P  >> $PS
#gmt6 pscoast  $Rmap $J -B -Na/.005p -Ia -A1000 -P -Sazure -Glightgrey -Di -O -W.01p -K >> $PS
#pscoast -Rg -J -B15g15 -Dc -A10000 -Glightgrey -P -O -W.01p -K >> map.ps
#gmt6 xyz2grd $Rmap $J -I.017 $topoxyz -Gtopo.grd -V
# gmt6 makecpt -Fgray -Cetopo1 -V > etopocolor.cpt # etopo1, -A

#gmt6 makecpt -C$oleron -T-700/700/50 > oleron1.cpt
#gmt6 makecpt -C$tpushuf -T0/900/50 > tpushuf.cpt
# gmt6 makecpt -C$afrikakarte -T-1320/720/10 > afrikakarte.cpt


#gmt6 makecpt -Cgray.cpt -T8/18/1 -Z -I > graay.cpt
# gmt6 makecpt -Cgray.cpt -T5/13/1 -Z -I > graay.cpt

# gmt6 makecpt -Cdark_red.cpt -T5/12/1 -Z -I > dark_reds.cpt
# gmt6 makecpt -Csandy_scar.cpt -T5/13/1 -Z  > sandy_scarr.cpt
# gmt6 makecpt -Credss.cpt -T5/13  > redss.cpt

#psxy stationSa.txt -Sa.52 -h0 -W0.3+cf -Cabc.cpt $J $Rmap -O -V -K >> $PS
# gmt6 grdimage $topogrd $J $Rmap -Bx -By -FRgray -Cafrikakarte.cpt -I+nt.85 -K -O >> $PS # original color
gmt6 grdimage $topogrd $J $Rmap -Bx -By -Cfes.cpt -I+nt.85 -K -O >> $PS # gray scale

gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -K -Di -O -W.1p >> $PS #-A10+l -Ia
#gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -Sazure -GFloralWhite  -K -Di -O -W.1p >> $PS #-GSeashell FloralWhite Ivory
#geological features
# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/CrustalBoundaries.txt $Rmap $J -W1.05p,gray -O -K -V >> $PS ### crustal boundary
# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/Gawler_Final.gmt $Rmap $J -W0.85p,darkblue,- -O -K >> $PS ### crustal boundary
# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/Tectonic_Zones.gmt $Rmap $J -W0.55p,yellow, -O -K >> $PS ### crustal boundary
# gmt6 psxy ~/Research/earthquakes/faults_ga_neales.txt $Rmap $J -: -W01.15p,yellow -O -K >> $PS ### neo_tectonics features
# gmt6 psxy ~/Research/earthquakes/faults_ga_karari.txt $Rmap $J -: -W01.15p,yellow -O -K >> $PS ### neo_tectonics features

#Stations
gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_BB.txt -: -W.05 -St.22 -Gdarkgreen@20 $J $Rmap -O -K >> $PS
gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_SP.txt -: -W.05 -St.22 -Gdarkgreen@20 $J $Rmap -O -K >> $PS
gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/per_st.txt -: -Ss.25 -Gblack $J $Rmap -O -K >> $PS
awk '{print $1,$2+.18,$3}' ~/Research/Lake_eyre_data/station/old_mix_stuff/per_st.txt | gmt6 pstext $Rmap $J -: -F+f3.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS
# awk '{print $1-.1,$2,$3}' ~/Research/Lake_eyre_data/station/stations_all_SP.txt | gmt6 pstext $Rmap $J -: -F+f2p,Helvetica-Bold -Gwhite -O -P -K >> $PS
# awk '{print $1-.1,$2,$3}' ~/Research/Lake_eyre_data/saation/stations_all_BB.txt | gmt6 pstext $Rmap $J -: -F+f2p,Helvetica-Bold -Gwhite -O -P -K >> $PS
####
########
gmt6 psxy $Rmap $J -W.95p,navy@5,- -V -O -K << EOF >> $PS
138.216202 -30.447201
134.027 -26
EOF

### Earthquake

# awk '{print $1,$2}' $eq_all_LE | gmt6 psxy $Rmap $J -: -Sc.05 -Gblack -B -O -P -K >> $PS #GA catalog all LE eq white circles
# for solid star with white boundary -Sa.32 -W.15,white -Gdarkred
#
#cat
# awk '{print $1,$2,$3,$4*.14}' $eq_NA_ml | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Cdark_reds.cpt -t20 -B -O -P -K >> $PS # eq
# awk '{print $1,$2,$4*.14}' $mine_NA_ml | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Gwhite@80 -B -O -P -K >> $PS # mine

# awk '{print $1,$2,$3,$4*.14}' $eq_NA_ml_2 | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Cdark_reds.cpt -t20 -B -O -P -K >> $PS # eq
# awk '{print $1,$2,$4*.14}' $mine_NA_ml_2 | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Gwhite@80 -B -O -P -K >> $PS # mine

awk '{print $1,$2,$3,$4*.14}' relocated_cross_sectn.txt | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Cdark_reds.cpt -t20 -B -O -P -K >> $PS # eq

# awk '{print $2,$1}' proj.txt | gmt6 psxy $Rmap $J -Sc.2 -W.15,black -Gblack -t90 -B -O -P -V -K >> $PS # eq

### GA_18_20_eq.txt
#awk '{print $1,$2,$3*.1}' GA_18_20_eq.txt | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Gwhite@10 -B -O -P -V -K >> $PS # mine
####
#######projection
# gmt6 project relocated_cross_sectn.txt -: -C135.688/-27.7944 -A320.83 -Fxypqrs > proj.txt
####
# draws ellipses
#echo 136.35 -28.38 320 5.5 2.2 | gmt6 psxy $Rmap $J -Se -W1,darkred,- -B -O -P -K >> $PS

#FocalMechanism #darkolivegreen for BW
# darkslategray for color topo

# gmt6 psmeca $Rmap $J -Sa.8c/6/u -W.5,gray18 -Gdarkolivegreen -Fr -V -M -Fp -C -O -K << END >> $PS
# # lon lat depth str dip slip mag plon plat # 21; 27;31; 40;
# 135.789 -26.198 9.567 339.7 80.9 32.2 1.969 136.2 -26.5 2 Ml
# 137.076 -25.664 6.350 153.9 74.4 -168.7 3.6788 138.1 -26.1 3.7 Ml
# 135.887 -28.227 9.674 303.6 75.7 15.2 1.212 135 -28.3 1.2 Ml
# 135.704 -27.882 10.934 339.5 77.8 37.3 3.22 136.3 -27.4 3.2 Ml
# 135.842 -28.101 9.931 325.4 63.3 15.5 1.6 134.4 -27.6 1.6 Ml
# 136.662 -27.568 9.934 339.7 79.8 27.6 1.4 137.3 -27.568 1.4 Ml
# END
#
echo 137.3 -28.4 Lake | gmt6 pstext $Rmap $J -F+f7p,Courier -O -P -K  >> $PS
echo 137.3 -28.55 Eyre | gmt6 pstext $Rmap $J -F+f7p,Courier -O -P -K  >> $PS

#gmt6 psxy $Rmap $J mine_sa.txt -Sc.3 -Gwhite -B -O -P -K >> $PS
# gmt6 psxy $Rmap $J mine_sa.txt -Sx.28 -W2.5,black -B -O -P -K >> $PS
#gmt6 psxy $Rmap $J mine_sa.txt -Sx.3 -W1.5,SlateGray -B -O -P -K >> $PS

# echo 133.65 -25.3 C-C\' | gmt6 pstext $Rmap $J -F+f7.5p,Helvetica,navy+a-60 -Gwhite@40 -O -P -K  >> $PS
####
gmt6 psbasemap -B -V $J $Rmap -O -P -K >> $PS
gmt6 gmtset FONT_ANNOT_PRIMARY 8p MAP_FRAME_PEN .8p FONT_LABEL 7.5p
# gmt6 psscale -Dx11.2c/1.4c+w2.6c/.24c+e -O -K -G-100/700 -Cafrikakarte.cpt -Bx200 -By+l"m"  >> $PS #-G-100/1000

gmt6 psscale -Dx11.2c/2.2c+w3c/.24c+e -O -K -Cdark_reds.cpt -Bx3 -Bx+l"Depth (km)" >> $PS
#gmt6 pslegend -Dx1.2c/4c+w2.4c/1.6c+o-1c/-.5c -F+gwhite+p.1 -O -K $J $Rmap << EOF >> $PS
gmt6 pslegend -Dx1.2c/2c+w1.8c/1.5c+o-1c/-.5c -F+gwhite+p.1 -O $J $Rmap << EOF >> $PS
S 0.2c t 0.3c darkgreen - 0.2i 5G
#S 0.2c t 0.3c SteelBlue - 0.2i SP
S 0.2c s 0.22c Black - 0.2i ANSN
S 0.2c c 0.15c darkred@40 - 0.2i Mlv 1
S 0.2c c 0.36c darkred@40 - 0.2i Mlv 3
#S 0.2c a 0.35c darkgreen - 0.2i Eq #2
# S 0.2c c 0.1c black .25p,black 0.2i GA eq
#S 0.2c x 0.2c - .95p,firebrick 0.2i Eq (R_19 < 200)
#S 0.2c x 0.2c - .95p,SlateGray  0.2i Mines
EOF

####LEGEND
# gmt6 pslegend -DJLT+w2.15c/.27c+o-2.15c/-.27c -F+gwhite@10+p.05 -O -K $J $Rmap << EOF >> $PS
# #S 0.2c s 0.3c DarkOrange - 0.2i 1Ds-1Dp
# EOF

# echo Oct'20 - May'22 | gmt6 pstext $J $Rmap -F+cTL+f7.5p -O -P -V >> $PS
# echo GA catalog 1900 - | gmt6 pstext $J $Rmap -F+cTL+f9p,darkpurple -O -P -V >> $PS


gmt6 ps2raster -A -Tj -E720 -P -Z -Vq $PS
open map_ml.jpg
# cp map.pdf ~/Dropbox/map_na_aust_ml1.pdf
###
