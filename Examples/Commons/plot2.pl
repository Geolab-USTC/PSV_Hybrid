#!/usr/bin/env perl
use strict;
use warnings;

@ARGV == 2 or die "perl plot2.pl dir1 dir2\n";

my ($dir1, $dir2) = @ARGV;

open(SAC, "|sac");
my $b0 = 899;
for (my $i=35; $i<=55; $i+=0.2) {
    $b = $b0 + $i;
    my $e = $b + 12;
    my $file = sprintf("%.2f.accz", $i);
    my $pdf = sprintf("%.1f.pdf", $i);
    print SAC "r $dir1/$file $dir2/$file \n";
    print SAC "bp c 1 3\n";
    print SAC "xlim $b $e\n";
    print SAC "div &1,depmax &2,depmax \n";
    print SAC "color black incre list black blue \n";
    print SAC "p2 \n";
    print SAC "saveimg $pdf\n";
}
print SAC "q\n";
close(SAC);
