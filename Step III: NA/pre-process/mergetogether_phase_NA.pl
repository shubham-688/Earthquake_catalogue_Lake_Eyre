#!usr/bin/perl -w
$mag = "2";
$num = 0;

# $dir = "../catalog";

$together = "phase_all_NA_format_mine.txt";
#Actually, it is the hypoDD input format
open(OT,">$together");
	$file = "phase_all/phase_sel_mine_2.txt";
	open(JK,"<$file");
	@par = <JK>;
	close(JK);
	foreach $file(@par){
		($test,$jk) = split(' ',$file);
		if($test =~ /^\d+$/){
			($jk,$year1,$mon1,$dd1,$time,$ot,$std,$lat,$lon,$dep,$mg,$mg_var,$P,$S,$PS,$st_gap) = split(' ',,$file);
			($hour,$min,$sec) = split('\:',$time);
			$num++;
			print OT "# $year1  $time $mon1  $dd1   $hour    $min    $sec    $std    $lat    $lon    $dep     $mg     $P     $S     $PS     $st_gap   $num\n";
		}else{
			($net,$station,$phase,$tt_abs,$tt_rel,$amplitude,$err,$weight,$azi) = split(' ',$file);
			print OT "$station $phase $tt_rel $err $azi \n";
		}
	}

close(OT);
