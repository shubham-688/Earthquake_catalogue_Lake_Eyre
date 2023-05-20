#!/bin/csh

set PS=map_ss.ps
#set Rmap=-Rg
#set Rmap=-R-135/40/-90/90
set Rmap=-R134.75/138/-30/-26 # Lake_eyre+ Flinders
#set J=-JM8/22
set J=-JM12/12 #sa
#set J=-Jx0.023i/0.029i

set topogrd = /Users/shubham/Research/etopo/etopo1.grd
#set topogrd = ~/Research/etopo/etopo1.grd
# set topoxyz = ~/Research/etopo/etopo1.xyz.

#set oleron = ~/Documents/ScientificColourMaps5/oleron/oleron.cpt
set tpushuf = ~/Documents/cpt-city/tp/tpushuf
set afrikakarte = ~/Documents/cpt-city/wkp/lilleskut/afrikakarte


set springs = ~/Research/geological_data_straya/hydrological_data/Vector_data/Hydrography/springs.gmt
set rivers = ~/Research/geological_data_straya/hydrological_data/Vector_data/Hydrography/watercourselines.gmt
set GAB = ~/Research/geological_data_straya/artesian_basin/81672/GAB_Hydrological_Boundary.gmt
set faults = ~/Research/geological_data_straya/shapefiles_surafce_geology_australia/faults.gmt
set tect_inl = ~/Research/Lake_eyre_data/geological_ft/Tectonic_Zones_Inlier.gmt

set eq_na_ml_F = relocated_FINAL_mlsorted.txt

# set eq_all_LE = ~/Research/earthquakes/ga_database/LE_allmag_alltime_FR.txt
set lake_eyre = ~/Research/geological_data_straya/hydrological_data/Vector_data/lake_eyre.gmt

gmt6 gmtset MAP_FRAME_TYPE plain MAP_FRAME_WIDTH 2p FONT_ANNOT_PRIMARY 11.5p MAP_TICK_LENGTH_PRIMARY 3p MAP_FRAME_PEN 1.2p
# gmt6 gmtset FONT_ANNOT_SECONDARY 6.5p

gmt6 psbasemap -BneWS -Bxa1f.25 -Bya1f.25 -V $J $Rmap -K >! $PS

#gmt6 makecpt -C$oleron -T-700/700/50 > oleron1.cpt
# gmt6 makecpt -C$tpushuf -T0/900/50 > tpushuf.cpt

gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -K -Di -O -W.1p >> $PS #-A10+l -Ia

#gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -Sazure -GFloralWhite  -K -Di -O -W.1p >> $PS #-GSeashell FloralWhite Ivory

# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/CrustalBoundaries.txt $Rmap $J -W0.35p,black,- -O -K >> $PS ### crustal boundary
# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/Gawler_Final.gmt $Rmap $J -W0.85p,darkblue,- -O -K >> $PS ### crustal boundary
gmt6 psxy $GAB $Rmap $J -W1.2p,navy,- -V -O -K  >> $PS # great artesian basin

gmt6 psxy $Rmap $J -W17p,darkseagreen@60 -V -O -K << EOF >> $PS
136.2 -28.65
135.25 -27.25
EOF

gmt6 psxy $Rmap $J -W17p,darkseagreen@60 -V -O -K << EOF >> $PS
136.2 -28.65
138.4 -30.5
EOF
gmt6 psxy $tect_inl $Rmap $J -W0.75p,dimgrey -Ggrey@20 -O -K >> $PS ### crustal boundary

gmt6 psxy  $Rmap $J -W0.008p,SteelBlue@10 $rivers -O -K  >> $PS # RIVERS

gmt6 psxy $faults $Rmap $J -W1.5p,black -O -K >> $PS ### faults
gmt6 psxy $springs -W.005,royalblue -Sa.25 -Groyalblue@25 $J $Rmap -O -K >> $PS # springs GAB


# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/Tectonic_Zones.gmt $Rmap $J -W0.55p,yellow, -O -K >> $PS ### crustal boundary

# gmt6 grdcontour moho.grd $J $Rmap -C2.5 -A2.5+f6p+ukm -S15 -Wa.5 -Wc.25 -T -O -P -K >> $PS # -S smooth factor

# awk '{print $2,$3}' ~/Research/Lake_eyre_data/station/marla.txt | gmt6 psxy -: -Si.11 -GKhaki $J $Rmap -O -K >> $PS ### Marla

# awk '{print $1,$2}' $eq_all_LE | gmt6 psxy $Rmap $J -: -Sc.05 -Gblack@40 -B -O -P -V -K >> $PS #GA catalog all LE eq white circles
# for solid star with white boundary -Sa.32 -W.15,white -Gdarkred
gmt6 psxy $lake_eyre $Rmap $J -W0.02p,dimgray -O -K -GsteelBlue@60 >> $PS # lake eyre boundary

#Stations
# gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_BB.txt -: -W.05 -St.22 -GDarkSlateBlue@20 $J $Rmap -O -K >> $PS
# gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_SP.txt -: -W.05 -St.22 -GDarkSlateBlue@20 $J $Rmap -O -K >> $PS
# gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/per_st.txt -: -Ss.25 -Gblack $J $Rmap -O -K >> $PS
# awk '{print $1,$2+.18,$3}' ~/Research/Lake_eyre_data/station/old_mix_stuff/per_st.txt | gmt6 pstext $Rmap $J -: -F+f3.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS

####
#awk '{if ($8+$9>107) print $8,$9}' $eq_real_all | gmt6 psxy $Rmap $J -: -Sa.5 -W.17,white -Gdarkred@30 -B -O -P -V -K >> $PS

awk '{print $1,$2,$3,$4*.2}' $eq_na_ml_F | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Cdark_reds.cpt -t10 -B -O -P -K >> $PS # eq

# awk '{print $1,$2,$3,$4*.2}' relocated_eq_na_aust_mlv.txt | gmt6 psxy $Rmap $J -: -Sc-Sa.5 -W.17,white -Cdark_reds.cpt -t5 -B -O -P -K >> $PS # eq

echo 136.2 -29.85 GAB | gmt6 pstext $Rmap $J -F+f16p,navy -O -P -K  >> $PS
echo 137.3 -28.4 Lake | gmt6 pstext $Rmap $J -F+f11p,Courier -O -P -K  >> $PS
echo 137.3 -28.55 Eyre | gmt6 pstext $Rmap $J -F+f11p,Courier -O -P -K  >> $PS

#faults
echo 135.88 -27.85 Kingston Fault | gmt6 pstext $Rmap $J -F+f10p,dimgray+a-40 -Gwhite@20 -V -O -P -K  >> $PS
echo 135.85 -28.65 Levi Fault | gmt6 pstext $Rmap $J -F+f10p,dimgray+a-60 -Gwhite@20 -V -O -P -K  >> $PS


# gmt6 psxy $Rmap $J mine_sa.txt -Sc.22 -Gwhite -B -O -P -K >> $PS
# # gmt6 psxy $Rmap $J mine_sa.txt -Sx.28 -W2.5,black -B -O -P -K >> $PS
# gmt6 psxy $Rmap $J mine_sa.txt -Sx.22 -W1.5,SlateGray -B -O -P -K >> $PS

gmt6 psscale -Dx12.2c/2.2c+w4c/.24c+e -O -K -Cdark_reds.cpt -Bx3 -Bx+l"Depth (km)" >> $PS
# gmt6 psscale -Dx11.2c/1.4c+w2.2c/.22c+e -O -K -G-100/700 -Cafrikakarte.cpt -Bx200 -By+l"m"  >> $PS #-G-100/1000
gmt6 gmtset FONT_ANNOT_PRIMARY 9.5p FONT_LABEL 13.5p

gmt6 pslegend -Dx11c/15c+w1.85c/2.05c+o-1c/-.5c -F+gwhite+p.1 -O $J $Rmap << EOF >> $PS
#S 0.2c - 0.3c - 1.5p,indianred 0.2i P axis
S 0.2c - 0.3c - 1.5p,Black 0.2i Faults
#S 0.2c s 0.22c Black - 0.2i ANSN
S 0.2c c 0.3c Firebrick@35 - 0.2i Eq
S 0.2c a 0.28c royalblue@15 - 0.2i Springs
S 0.2c s 0.3c grey@10 - 0.2i DPI
S 0.2c s 0.3c darkseagreen@10 - 0.2i NFZ
EOF


# gmt6 psconvert -A -Tf -P -Z $PS

gmt6 ps2raster -A -Tj -E720 -P -Z -Vq $PS # f is for Pdf; j JPEG
open map_ss.jpg
# cp map.pdf ~/Dropbox/map.pdf
###
