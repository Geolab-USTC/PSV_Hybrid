#!/usr/bin/perl -w
use strict;
use warnings;

@ARGV==1 or die "  Usage: perl $0 snapfile\n";

my ($Snapfile) = @ARGV;
my ($xmin, $xmax, $ymin, $ymax, $zmin, $zmax) = split " ", `minmax -C $Snapfile`;
my $CPT="div.pal";
my $PS=$Snapfile.".ps";
my $GRD=$Snapfile.".grd";
my $J="X21c/-7c";
my $R="$xmin/$xmax/$ymin/$ymax";

my $ABSMAX=max(abs($zmin),abs($zmax));
my $B="a1000/a500:.'Normlized By $ABSMAX':WSEN";

open(SNAP,">snapnew");
open(IN,"$Snapfile");
while(my $line = <IN>){
	chomp($line);
	my ($i,$j,$val) = split(" ",$line);
	$val /= $ABSMAX;
	print SNAP "$i $j $val\n"
}
close(IN);
close(SNAP);

system("xyz2grd  snapnew -G$GRD -R$R -I5/5 -V");
system("grdimage $GRD -R$R -J$J -C$CPT -B$B -K > $PS");
system("grdcontour $GRD  -R$R -J$J -C0.1 -L0.01/1 -W0.5p/0/0/0 -K -O -V >> $PS");
system("grdcontour $GRD  -R$R -J$J -C0.1 -L-1/-0.01 -W0.5p/0/0/0 -K -O -V >> $PS");
system("psscale -D3.5i/-0.5i/6i/0.1ih -C$CPT -Ba5g5 -O >> $PS");

unlink $GRD;
unlink "snapnew";
unlink glob ".gmt*";

sub max{
	my ($v1,$v2) = @_;
	if($v1>$v2){return $v1;}
	else {return $v2;}
}
