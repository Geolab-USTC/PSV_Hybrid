#!/usr/bin/env perl
#
# Plot velocity structure of FD region
#
use strict;
use warnings;

@ARGV==1 or die "  Usage: perl $0 modelfile\n";

my ($file) = @ARGV;

die "file not found!\n" if !-e $file;

my ($xmin, $xmax, $ymin, $ymax, $zmin, $zmax) = split " ", `minmax -C $file`;
my $PS=$file.".ps";
my $GRD=$file.".grd";
my $J="X21c/-7c";
my $R="$xmin/$xmax/$ymin/$ymax";
my $CPT = "vel.cpt";
my $B="a50f10/a50f10";

system("xyz2grd $file -G$GRD -R$R -I1/1 -V");
system("grd2cpt $GRD -Crainbow -Z > $CPT");
system("grdimage $GRD -R$R -J$J -C$CPT -B$B -K -Y2i> $PS");
system("psscale -D3.5i/-0.5i/6i/0.1ih -C$CPT -Ba500 -O >> $PS");

unlink $GRD;
unlink glob ".gmt*";
unlink $CPT;
