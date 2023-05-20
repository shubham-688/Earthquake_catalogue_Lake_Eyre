#!/bin/csh
#ogr2ogr -f "GMT" aus5wgd_l.gmt  aus5wgd_l.shp -t_srs WGS84
set PS=map_shm.ps
#set Rmap=-Rg
#set Rmap=-R-135/40/-90/90
set Rmap=-R130/142/-32.5/-25 # Lake_eyre+ Flinders
set Rmap=-R118/150/-35/-10 # Lake_eyre+ Flinders

set Rmap=-R108/-42/152/-9r

#set J=-JM8/22
set J=-JM12/11 #sa
set J=-JS130/-30/4i #oz
#set J=-Jx0.023i/0.029i

# set topogrd = /Users/shubham/Research/etopo/etopo_sa_15s.grd
set topogrd = ~/Research/etopo/ETOPO1_Ice_g_gmt5.grd
# set topoxyz = ~/Research/etopo/etopo1.xyz.
set GAB = ~/Research/geological_data_straya/artesian_basin/81672/GAB_Hydrological_Boundary.gmt

set oleron = ~/Documents/ScientificColourMaps5/oleron/oleron.cpt
set tpushuf = ~/Documents/cpt-city/tp/tpushuf
set afrikakarte = ~/Documents/cpt-city/wkp/lilleskut/afrikakarte

gmt6 gmtset MAP_FRAME_TYPE plain MAP_FRAME_WIDTH 4p FONT_ANNOT_PRIMARY 9.5p MAP_TICK_LENGTH_PRIMARY 7p MAP_FRAME_PEN 2.2p
gmt6 gmtset COLOR_BACKGROUND lightsteelblue MAP_ANNOT_OBLIQUE 6 # MAP_ANNOT_OBLIQUE lat_horizontal

gmt6 psbasemap -BneWS -Bxa10f5 -Bya10f5 $J $Rmap -K >! $PS
gmt6 makecpt -C$afrikakarte -M -T-1620/890/20 > afrikakarte.cpt #  for transparency

###########
gmt6 grdimage $topogrd $J $Rmap -Bx -By -FRgray -Cafrikakarte.cpt -I+nt.85 -K -O >> $PS # original color -I+nt.5
# gmt6 grdimage $topogrd $J $Rmap -Bx -By -Cfes.cpt -I+nt.85 -K -O >> $PS # gray scale

gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -K -Di -O -W.1p >> $PS #-A10+l -Ia
########
gmt6 psxy $GAB $Rmap $J -W0.3p,navy -Gdarkgrey@70 -O -K >> $PS ###

gmt6 psxy $Rmap $J -W1p,white,- -K -P -O << EOF >> $PS #131/141.5/-32.5/-24.5
131 -32.5
131 -24.5
141.5 -24.5
141.5 -32.5
131 -32.5
EOF

#awk '{if ($6 != "D" && $6 != "E" ) {print $0}} ' wsm_aus.txt
# awk '{if ($6 == "D" || $6 == "E") print $2,$1,$3,.3c}' wsm_aus.txt | gmt6 psxy $Rmap $J  -SV.3c -W.65,dimgray -V -B -O -P -K >> $PS # eq
# awk '{if ($6 == "D" || $6 == "E") print $2,$1,$3+180,.3c}' wsm_aus.txt | gmt6 psxy $Rmap $J  -SV.3c -W.65,dimgray -V -B -O -P -K >> $PS # eq
#
echo 121.5 -28 7cm/yr | gmt6 pstext $Rmap $J -F+f7p,darkred+a-10 -Gwhite@50 -V -O -P -K  >> $PS

# APM wrt to africa
gmt6 psxy $J $Rmap -SV.125i+ea -P -O -K -Gmaroon -W1.5p,maroon << EOF >> $PS
122 -27 15 .45i
EOF
gmt6 psxy $J $Rmap -SV.125i+ea -P -O -K -Gmaroon -W1.5p,maroon << EOF >> $PS
135 -41 12 .45i
EOF
gmt6 psxy $J $Rmap -SV.125i+ea -P -O -K -Gmaroon -W1.5p,maroon << EOF >> $PS
138 -16 12 .45i
EOF
#officer
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
135 -28 95 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
135 -28 -95 .25i
EOF
#flinders
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
138 -32.5 91 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
138 -32.5 271 .25i
EOF
#
#darling
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
143 -30.5 91 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
143 -30.5 271 .25i
EOF
#cooper
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
140.9 -26.5 100 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
140.9 -26.5 280 .25i
EOF
#cooper
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
140.9 -22.25 87 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
140.9 -22.25 267 .25i
EOF
# long lat azi length size (head)
#amadeus
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
131 -24 27 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
131 -24 207 .25i
EOF
#McA
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
134 -16.5 26 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
134 -16.5 206 .25i
EOF
#
#SE-gipps
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
149 -38 135 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
149 -38 315 .25i
EOF
#SE-otway
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
142 -38 133 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
142 -38 313 .25i
EOF
#syd
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
150 -32 57 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
150 -32 237 .25i
EOF
#bowen-surat2
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
150 -25.5 65 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
150 -25.5 245 .25i
EOF
#bowen-surat1
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
149 -21.5 32 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
149 -21.5 212 .25i
EOF
#### Western A
#perth
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
115 -30.8 88 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
115 -30.8 268 .25i
EOF
#canning
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
124.5 -16.8 53 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
124.5 -16.8 233 .25i
EOF
#N Carnanvon
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
117 -20 108 .25i
EOF
gmt6 psxy $J $Rmap -SV.125i+be -P -O -K -Gblack -W1.5p << EOF >> $PS
117 -20 288 .25i
EOF
#
echo 139.5 -28 Lake | gmt6 pstext $Rmap $J -F+f6p,Courier -O -P -K  >> $PS
echo 139.5 -28.75 Eyre | gmt6 pstext $Rmap $J -F+f6p,Courier -O -P -K  >> $PS


gmt6 psbasemap -B -V $J $Rmap -O -P -K >> $PS

 #-A10+l -Ia
gmt6 gmtset FONT_ANNOT_PRIMARY 8p MAP_FRAME_PEN .8p FONT_LABEL 7.5p
gmt6 psscale -Dx1c/1c+w2.6c/.24c+e+h -O -G-100/880 -F+gwhite+p.1 -Cafrikakarte.cpt -Bx250 -By+l"(m)" -K >> $PS #-G-100/1000

# gmt6 psscale -Dx11.2c/2.2c+w4c/.24c+e -O -K -Cdark_reds.cpt -Bx3 -Bx+l"Depth (km)" >> $PS
#gmt6 pslegend -Dx1.2c/4c+w2.4c/1.6c+o-1c/-.5c -F+gwhite+p.1 -O -K $J $Rmap << EOF >> $PS
gmt6 pslegend -Dx1.2c/8c+w2.2c/1.2c+o-1c/-.75c -F+gwhite+p.1 -O $J $Rmap << EOF >> $PS
S 0.2c v 0.5c black 1.0p,black 0.2i mean S@-hmax@-
S 0.2c v 0.5c maroon 1.0p,maroon 0.2i APM
S 0.2c s 0.22c grey - 0.2i GAB
# S 0.2c - 0.3c - 1.25p,dimgray 0.2i S@-hmax@- (D-E)
EOF

####LEGEND
# gmt6 pslegend -DJLT+w2.15c/.27c+o-2.15c/-.27c -F+gwhite@10+p.05 -O $J $Rmap << EOF >> $PS
# # #S 0.2c s 0.3c DarkOrange - 0.2i 1Ds-1Dp
# EOF

#pstext $J $Rmap -P -O -K <<End >> $PS
#-40 -7 18 0 31 1 (b)
#End

gmt6 ps2raster -A -Tj -E720 -P -Z -Vq $PS
open map_shm.jpg
# cp map.pdf ~/Dropbox/map_na_aust_ml1.pdf
###
