#!/usr/bin/perl -w
# see https://github.com/Dal-mzhang/REAL for details about REAL. This script find events using P-S picks from STEP I, month-by-month.
$year = "2021";
$mon = "1"; # change here for month
$Jday = "1"; # change line below for number of days in a month
for ($day = "1"; $day < 32; $day++){
    if ($Jday <100) {$Jday= sprintf("%03d", $Jday);}
    if ($day <10) {$day= sprintf("%02d", $day);}
    if ($mon <10) {$mon= sprintf("%02d", $mon);}
    print"Day/Month: $day/$mon \n";
    print"Doing day $Jday\n";

    $D = "$year/$mon/$day";
    $R = "4.5/10/0.15/1000/90/360"; #deg/km/deg/km/sec; search range in deg/km search grid deg/km iint in sec
    $G = "5.5/10/0.05/10"; #  (deg/km/deg/km) h-range in tt table; grid size in tt
    # $G = "1.4/20/0.01/1"; # original REAL parameters

    $V = "6.35/3.7"; # average P and S velocity ori 6.2/3.5# for the 1d model AUREM 6.35/3.7
    $S = "7/0/9/1/2/2/2.25/1.5"; # here to change no of P and S picks
    #(np0/ns0/nps0/nppp(1)/std0/dtps/nrt[/rsel/ires])

    # std0 and nrt affects the first pick a lot
    # rsel doesn't affect a lot..2 and 4 were kidna same..just some less S picks in rsel=2

    $dir = "../Pick/2021_1/$year$mon$day"; # change here yearly
    $station = "../Data/stations_19.txt";
    # $station = "../Data/station.dat";
    # $ttime = "./tt_db/ttdb.txt";
    $ttime = "./tt_db/ttdb_moho45.txt";
    # $ttime = "./tt_db/itvel_moho35.txt";

    system("REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime");
    print"REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime\n";
    print"-----------------------------------------------------------------------\n";


    $cat_name="./catalog/jun/catalog_sel_$Jday.txt"; #change here
    $cat="catalog_sel.txt";
    rename $cat,$cat_name;
    $phase_name="./catalog/jun/phase_sel_$Jday.txt"; #change here
    $phase="phase_sel.txt";
    rename $phase,$phase_name;
    $Jday++;
  }

# np0: threshold for number of P picks
# ns0: threshold for number of S picks
#
# nps0: threshold for total number of picks (P&S)
#
# nppp: effective number of grids that meet the thresholds (np0, ns0 and nps0)
# [The recommend value is 1. nppp>1 means more critical threshold (eligible for fine grid size)]
#
# std0: standard deviation threshold (residual).[Here residual is defined as the deviation of the
# origin times from different picks (for the same event). It is not the traditional RMS residual]
#
# dtps: time threshold for S and P separation. [dtps is used to remove some false S picks.
# P picks appear in S pick pool in real data when applying STA/LTA pickers. Set dtps = 0 to turn this constraint off]
#
# nrt: nrt*default time window (>1) [i.e., nrt*sqrt(tdx**2+tdh**2)/vp0 for P time
# window or nrt*sqrt(tdx**2+tdh**2)/vs0 for S time window. It accommodates the inaccuracy
# of velocity model (as well as pick uncertainty). Use larger nrt if velocity model is insufficient]
#
# rsel: tolerance multiplier; keep picks in phase_sel.txt with residuals less than rsel*STD
# [rsel is used to remove picks with large residuals. Default rsel is 5.0
# (a large value, i.e., approximately turn this constraint off).
# For example, if your final standard deviation is 0.1 sec,
