#!/bin/csh
#ogr2ogr -f "GMT" aus5wgd_l.gmt  aus5wgd_l.shp -t_srs WGS84
set PS=map.ps
#set Rmap=-Rg
#set Rmap=-R-135/40/-90/90
set Rmap=-R133/140/-30.75/-25.5 # Lake_eyre+ Flinders
#set J=-JM8/22
set J=-JM12/11 #sa
#set J=-Jx0.023i/0.029i

set topogrd = /Users/shubham/Research/etopo/etopo_sa_15s.grd
#set topogrd = ~/Research/etopo/etopo1.grd
# set topoxyz = ~/Research/etopo/etopo1.xyz.

set oleron = ~/Documents/ScientificColourMaps5/oleron/oleron.cpt
set tpushuf = ~/Documents/cpt-city/tp/tpushuf
set afrikakarte = ~/Documents/cpt-city/wkp/lilleskut/afrikakarte

#set eq_real_all_18 = ~/Research/Lake_eyre_data/sAus_eq/picker/Real/REAL/catalog/final_cat_good/cat_sg200_2018_19.txt
#set eq_real_all_19 = ~/Research/Lake_eyre_data/sAus_eq/picker/Real/REAL/catalog/final_cat_good/cat_sg200_2019.txt
#set eq_real_all_20 = ~/Research/Lake_eyre_data/sAus_eq/picker/Real/REAL/catalog/final_cat_good/cat_sg200_2020.txt

set eq_all_LE = ~/Research/earthquakes/ga_database/LE_allmag_alltime_FR.txt
set eq_NA = ~/Dropbox/NA_temp/relocated_eqs/relocated_eq_na_aust.txt
set eq_NA_ml = ~/Dropbox/NA_temp/mag_esti/relocated_eq_na_aust_mlv.txt
set mine_NA = ~/Dropbox/NA_temp/relocated_eqs/relocated_mine_na_aust.txt
set mine_NA_ml = ~/Dropbox/NA_temp/mag_esti/relocated_mine_na_aust_mlv.txt

set eq_na_ml_F = relocated_FINAL_mlsorted.txt
set mine_na_ml_F = relocated_mine_na_ml_FINAL.txt


gmt6 gmtset MAP_FRAME_TYPE plain MAP_FRAME_WIDTH 4p FONT_ANNOT_PRIMARY 9.5p MAP_TICK_LENGTH_PRIMARY 4p MAP_FRAME_PEN 2.2p


gmt6 psbasemap -BneWS -Bxa1f.25 -Bya1f.25 $J $Rmap -K >! $PS

gmt6 makecpt -Cdark_red.cpt -T5/12/1 -Z -I > dark_reds.cpt


# gmt6 grdimage $topogrd $J $Rmap -Bx -By -FRgray -Cafrikakarte.cpt -I+nt.85 -K -O >> $PS # original color
gmt6 grdimage $topogrd $J $Rmap -Bx -By -Cfes.cpt -I+nt.85 -K -O >> $PS # gray scale

gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -K -Di -O -W.1p >> $PS #-A10+l -Ia


# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/CrustalBoundaries.txt $Rmap $J -W1.05p,gray -O -K -V >> $PS ### crustal boundary
# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/Gawler_Final.gmt $Rmap $J -W0.85p,darkblue,- -O -K >> $PS ### crustal boundary

# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/Tectonic_Zones.gmt $Rmap $J -W0.55p,yellow, -O -K >> $PS ### crustal boundary

# gmt6 psxy ~/Research/earthquakes/faults_ga_neales.txt $Rmap $J -: -W01.15p,yellow -O -K >> $PS ### neo_tectonics features
# gmt6 psxy ~/Research/earthquakes/faults_ga_margaret.txt $Rmap $J -: -W01.15p,yellow -O -K >> $PS ### neo_tectonics features
# gmt6 psxy ~/Research/earthquakes/faults_ga_cober.txt $Rmap $J -: -W01.15p,yellow -O -K >> $PS ### neo_tectonics features
# gmt6 psxy ~/Research/earthquakes/faults_ga_karari.txt $Rmap $J -: -W01.15p,yellow -O -K >> $PS ### neo_tectonics features
#Stations
gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_BB.txt -: -W.05 -St.22 -GDarkSlateBlue@20 $J $Rmap -O -K >> $PS
gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_SP.txt -: -W.05 -St.22 -GDarkSlateBlue@20 $J $Rmap -O -K >> $PS
gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/per_st.txt -: -Ss.25 -Gblack $J $Rmap -O -K >> $PS
awk '{print $1,$2+.18,$3}' ~/Research/Lake_eyre_data/station/old_mix_stuff/per_st.txt | gmt6 pstext $Rmap $J -: -F+f3.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS

#######

awk '{print $1,$2}' $eq_all_LE | gmt6 psxy $Rmap $J -: -Sc.035 -Gblack -B -O -P -K >> $PS #GA catalog all LE eq white circles

awk '{print $1,$2,$3,$4*.12}' $eq_na_ml_F | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Cdark_reds.cpt -t10 -B -O -P -K >> $PS # eq
# awk '{print $1,$2,$4*.12}' $mine_na_ml_F | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Gwhite@80 -B -O -P -K >> $PS # mine
### GA_18_20_eq.txt
#awk '{print $1,$2,$3*.1}' GA_18_20_eq.txt | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Gwhite@10 -B -O -P -V -K >> $PS # mine
####
#FocalMechanism #darkolivegreen for BW
# darkslategray for color topo
###round_1 fms
gmt6 psmeca $Rmap $J -Sa.6c/5/u -W.5,gray18 -Gdarkolivegreen -Fr -V -M -Fp -C -O -K << END >> $PS
# lon lat depth str dip slip mag plon plat # 21;31; 40;45;46
135.789 -26.198 9.567 342.6 79.4 30.2 1.969 136.2 -26.5 2 Ml
135.887 -28.227 9.674 303.6 75.7 15.2 1.212 135 -28.3 1.2 Ml
135.704 -27.882 10.934 339.5 77.8 37.3 3.22 136.3 -27.4 3.2 Ml
135.842 -28.101 9.931 323.3 65.6 8.2 1.6 134.4 -27.6 1.6 Ml
136.662 -27.568 9.934 334.9 82.8 32.5 1.4 136.8 -27 1.4 Ml
END
gmt6 psmeca $Rmap $J -Sa.6c/5/u -W.5,gray18 -Gdarkkhaki -Fr -V -M -Fp -C -O -K << end >> $PS
# 27 (1);
137.076 -25.664 6.350 153.0 76.4 -170.4 3.6788 138 -26.23 3.7 Ml
137.076 -25.664 6.350 106.0 55.1 107.6 3.6788 138.1 -25.8
end
#round_2 fms
gmt6 psmeca $Rmap $J -Sa.6c/5/u -W.5,gray18 -Gdarkolivegreen -Fr -V -M -Fp -C -O -K << eND >> $PS
#10, 18, 38, 68, 77, 78
135.703 -28.100 8.873 312.6 58 2 1.844 135.2 -29 1.8 Ml
136.310 -28.695 9.626 123.9 78.4 4.6 1.998 135.8 -29.2 2 Ml
135.941 -28.226 9.917 115.5 86.9 -31.15 1.828 137 -27.8 1.8 Ml
135.639 -28.855 9.115 301.9 53 4 1.725 134.6 -30 1.7 Ml
136.193 -28.499 9.897 304.3 62.4 26.3 2.675 136 -30 2.7 Ml
136.209 -28.502 10.320 316.6 75.1 35.8 2.113 136.5 -29.5 2.1 Ml
eND
gmt6 psmeca $Rmap $J -Sa.6c/5/u -W.5,gray18 -Gdarkkhaki -Fr -V -M -Fp -C -O -K << enD >> $PS
#41, 51, 59,
136.442 -28.243 5.166 275.5 77.6 158.7 1.973 137.7 -28.1 2 Ml
136.468 -28.489 9.749 171.3 72.7 1.9 1.779 138 -28.9 1.8 Ml
136.051 -28.414 10.504 185.3 77.4 -33.6 1.873 134 -28.6 1.9 Ml
enD

echo 137.3 -28.4 Lake | gmt6 pstext $Rmap $J -F+f7p,Courier -O -P -K  >> $PS
echo 137.3 -28.55 Eyre | gmt6 pstext $Rmap $J -F+f7p,Courier -O -P -K  >> $PS

#gmt6 psxy $Rmap $J mine_sa.txt -Sc.3 -Gwhite -B -O -P -K >> $PS
# gmt6 psxy $Rmap $J mine_sa.txt -Sx.28 -W2.5,black -B -O -P -K >> $PS
#gmt6 psxy $Rmap $J mine_sa.txt -Sx.3 -W1.5,SlateGray -B -O -P -K >> $PS

gmt6 psbasemap -B -V $J $Rmap -O -P -K >> $PS

# gmt6 psscale -Dx10c/.6c+w1.5c/.18c+e -O -Cgraay.cpt -Bx5 -By+l"Depth (km)" -K -P >> $PS
#################
#pstext station_Sa.txt $Rmap $J -P -F+f6p,Helvetica,=0.6p,darkgreen -O -B -K >> map.ps

 #-A10+l -Ia
gmt6 gmtset FONT_ANNOT_PRIMARY 8p MAP_FRAME_PEN .8p FONT_LABEL 7.5p
# gmt6 psscale -Dx11.2c/1.4c+w2.6c/.24c+e -O -K -G-100/700 -Cafrikakarte.cpt -Bx200 -By+l"m"  >> $PS #-G-100/1000

gmt6 psscale -Dx11.2c/2.2c+w4c/.24c+e -O -K -Cdark_reds.cpt -Bx3 -Bx+l"Depth (km)" >> $PS
#gmt6 pslegend -Dx1.2c/4c+w2.4c/1.6c+o-1c/-.5c -F+gwhite+p.1 -O -K $J $Rmap << EOF >> $PS
gmt6 pslegend -Dx1.2c/2c+w1.8c/1.8c+o-1c/-.5c -F+gwhite+p.1 -O $J $Rmap << EOF >> $PS
S 0.2c t 0.3c DarkSlateBlue - 0.2i 5G
#S 0.2c t 0.3c SteelBlue - 0.2i SP #S 0.2c t 0.3c darkgreen - 0.2i 5G
S 0.2c s 0.22c Black - 0.2i ANSN
S 0.2c c 0.15c darkred@40 - 0.2i Mlv 1
S 0.2c c 0.36c darkred@40 - 0.2i Mlv 3
S 0.2c c 0.1c black .25p,black 0.2i GA eq
#S 0.2c a 0.35c darkgreen - 0.2i Eq #2
# S 0.2c c 0.1c black .25p,black 0.2i GA eq
#S 0.2c x 0.2c - .95p,firebrick 0.2i Eq (R_19 < 200)
#S 0.2c x 0.2c - .95p,SlateGray  0.2i Mines
EOF

####LEGEND
# gmt6 pslegend -DJLT+w2.15c/.27c+o-2.15c/-.27c -F+gwhite@10+p.05 -O $J $Rmap << EOF >> $PS
# # #S 0.2c s 0.3c DarkOrange - 0.2i 1Ds-1Dp
# EOF

#pstext $J $Rmap -P -O -K <<End >> $PS
#-40 -7 18 0 31 1 (b)
#End

gmt6 ps2raster -A -Tj -E720 -P -Z -Vq $PS
open map.jpg
# cp map.pdf ~/Dropbox/map_na_aust_ml1.pdf
###
