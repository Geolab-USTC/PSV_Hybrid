#!/usr/bin/env perl
#
# The Standard PREM model prem.nd comes from the TauP package,
# this file gives parameters of specified depth.
#
# This script convert prem.nd to layered homogeneous model.
#
use strict;
use warnings;

my $in = "prem.nd";
my $out = "prem";

open(IN, "< $in") or die "Error in opening file $in\n";
my $nlay = 0;               # total number of layers
my ($moho, $cmb, $icb);     # major boundaries
my $DH;                     # maximum thickness for each layer
my (@dep, @th, @vp, @vs, @rh, @Qp, @Qs);

# read prem model in
foreach (<IN>) {
    my @lines = split ' ', $_;
    if ( @lines == 1 ) {    # mantle | outer-core | inner-core
        if ($lines[0] =~ 'mantle') {
            $moho = $nlay;
        } elsif ($lines[0] =~ 'outer-core') {
            $cmb = $nlay;
        } elsif ($lines[0] =~ 'inner-core') {
            $icb = $nlay;
        } else {
            die "Error in parsing major boundary!\n";
        }
        next;
    }
    ($dep[$nlay], $vp[$nlay], $vs[$nlay], $rh[$nlay], $Qp[$nlay], $Qs[$nlay]) = @lines;
    $nlay++;
}
close(IN);

print STDERR "Moho=$moho CMB=$cmb ICB=$icb\n";

open(OUT, "> $out") or die "Error in opening $out to write.\n";

# Free surface => Moho
$DH = 20;
for (my $i=1; $i<$moho; $i++) {
    &discret($i, $DH);
}

# Moho => CMB
for (my $i=$moho; $i<$cmb; $i++) {
    &discret($i, $DH);
}

# CMB => ICB
$DH = 10;
for (my $i=$cmb; $i<$icb; $i++) {
    &discret($i, $DH);
}

# ICB => Center of the Earth
$DH = 10;
#for (my $i=$icb; $i<$nlay; $i++) {
#    &discret($i, $DH);
#}
$DH = 1;
for (my $i=$icb; $i<$icb+3; $i++) {
    &discret($i, $DH);
}
$DH = 10;
for (my $i=$icb+3; $i<$nlay; $i++) {
    &discret($i, $DH);
}
close(OUT);

sub discret{
    my ($lay, $dh) = @_;
    my $th = $dep[$lay] - $dep[$lay-1];     # thickness
    return if $th == 0;                     # This is a sharp boundary

    my $nn = int($th/$dh);                  # number of sublayer with thickness=$dh

    my $i;
    for ($i=0; $i<$nn; $i++) {
        my $h = $dh/2.0 + $i*$dh;           # depth relative to top interface
        my $vp = &linear($vp[$lay], $vp[$lay-1], $th, $h);
        my $vs = &linear($vs[$lay], $vs[$lay-1], $th, $h);
        my $rh = &linear($rh[$lay], $rh[$lay-1], $th, $h);
        my $Qp = &linear($Qp[$lay], $Qp[$lay-1], $th, $h);
        my $Qs = &linear($Qs[$lay], $Qs[$lay-1], $th, $h);

        $vs = 0.001 if $vs < 1.0e-6;        # vs has a tiny value
        printf OUT "%15.6f %15.6f %15.6f %15.6f %10.1f %10.1f\n", $dh, $vp, $vs, $rh, $Qp, $Qs;
    }

    # last sublayer
    return if $th-($i*$dh) <= 1.0e-6;       # thickness of last sublayer is 0
    my $dep = ($i*$dh+$th)/2.0;
    my $vp = &linear($vp[$lay], $vp[$lay-1], $th, $dep);
    my $vs = &linear($vs[$lay], $vs[$lay-1], $th, $dep);
    my $rh = &linear($rh[$lay], $rh[$lay-1], $th, $dep);
    my $Qp = &linear($Qp[$lay], $Qp[$lay-1], $th, $dep);
    my $Qs = &linear($Qs[$lay], $Qs[$lay-1], $th, $dep);
    $th -= ($i*$dh);
    $vs = 0.001 if $vs < 1.0e-6;        # vs has a tiny value
    printf OUT "%15.6f %15.6f %15.6f %15.6f %10.1f %10.1f\n", $th, $vp, $vs, $rh, $Qp, $Qs;
}

sub linear{
    my ($v2, $v1, $th, $h) = @_;        # bottom interface, top interface, dh, depth
    my $slop = ($v2-$v1)/$th;
    my $value = $slop*$h + $v1;         # y = (v2-v1)/h*x + v1
    return $value;
}
