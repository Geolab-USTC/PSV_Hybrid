#!/usr/bin/env perl
#
# Perl script to setup 2D FD model with dipping edge
#                ______________
#               /              \
#              /                \
#             /                  \
#  __________/                    \___________
# x1        x2  x3             x4 x5         x6
#
use strict;
use warnings;

my $model = "fdmodel";  # filename for output fd model
my $np  =   6;          # six points to define a 2D structure

my $xstart1 = -1;       # start from -1
my $xstart2 = 2400;     # start point of left-side dipping edge
my $xstart3 = 2420;     # end point of left-side dipping edge
my $xstart4 = 2680;     # start point of right-side dipping edge
my $xstart5 = 2700;     # end point of right-side dipping edge
my $xstart6 = 80000;    # go to infinite distance

my $thick   = 73;       # thickness of dipping structure (flattened)

my $ystart1 = 445;      # depth relative to GRT-FD interface
my $ystart2 = $ystart1;
my $ystart3 = $ystart2-$thick;
my $ystart4 = $ystart3;
my $ystart5 = $ystart1;
my $ystart6 = $ystart1;

# 2D velocity structure
my $vptop  =   52.112;
my $vstop  =   0.005;
my $dentop =   61.227;
my $vpbottom = 57.545;
my $vsbottom = 18.286;
my $denbottom = 66.6;

my $layer = 3;     # total layer number

open(Model, "> $model") or die "Error in opening file $model\n";
#print Model "-$layer\n";
print Model "$layer\n";

##set up the left-side region
# layer 0
print Model "$np $vptop $vstop $dentop\n";
print Model "$xstart1 $xstart2 $xstart3 $xstart4 $xstart5 $xstart6\n";
print Model "$ystart1 $ystart2 $ystart3 $ystart4 $ystart5 $ystart6\n";

# layer 1-19
for(my $tt=1; $tt<$layer-1; $tt++){
    my $x3=$xstart3-($xstart3-$xstart2)*$tt/($layer-1);
    my $y3=$ystart3-($ystart3-$ystart2)*$tt/($layer-1);

    my $x4=$xstart4+($xstart5-$xstart4)*$tt/($layer-1);
    my $y4=$y3;
    my $vp=$vptop-($vptop-$vpbottom)*($tt-1)/($layer-1);
    my $vs=$vstop-($vstop-$vsbottom)*($tt-1)/($layer-1);
    my $den=$dentop-($dentop-$denbottom)*($tt-1)/($layer-1);

    print Model "$np $vp $vs $den\n";
    print Model "$xstart1 $xstart2 $x3 $x4 $xstart5 $xstart6\n";
    print Model "$ystart1 $ystart2 $y3 $y4 $ystart5 $ystart6\n";
}

# layer 20
print Model "2 $vpbottom $vsbottom $denbottom\n";
print Model "$xstart1 $xstart6\n";
print Model "$ystart1 $ystart6\n";

close(Model);
