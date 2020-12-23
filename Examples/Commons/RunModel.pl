#!/usr/bin/env perl
use warnings;
use strict;
#
#   perl shell for running the hybrid method and Kirchhoff
#

my $BIN      = '../../bin';
my $PARFILE  = "testprem.par";     # parameters file

my $debug = 0;      # run check_* to debug

my $GRTWKMrun            = "yes";
my $GRTrun               = "yes";
my $FDrun                = "yes";
my $FDdemult2kir         = "yes";
my $KIRrun               = "yes";

# epicentral distance in degree
my @Recivers;
my $nrcv = 0;
for (my $i=35; $i<=55; $i+=1.0) {
    $Recivers[$nrcv] = $i;
    $nrcv++;
}

if($GRTWKMrun eq "yes"){
    printf STDERR "###########################################################\n";
    printf STDERR "# Create GRT library from receiver to Kirchhoff interface #\n";
    printf STDERR "###########################################################\n";
    my $REFMODEL = "prem.grt";
    system("$BIN/aseriesCMB par=$PARFILE model=$REFMODEL");

    if ($debug) {
        printf STDERR "##############################################################\n";
        printf STDERR "# Debug for GRT library from receiver to Kirchhoff interface #\n";
        printf STDERR "##############################################################\n";
        mkdir "WKM", 0755 if !-e "WKM";
        system "$BIN/check_wkm par=$PARFILE";
    }
}

if($GRTrun eq "yes"){
    printf STDERR "####################################################\n";
    printf STDERR "# Running GRT calculations for FD input wavefields #\n";
    printf STDERR "####################################################\n";
    system("$BIN/aserpsvfd par=$PARFILE");

    if ($debug) {
        printf STDERR "#########################################################\n";
        printf STDERR "# Debug for GRT library from source to GRT-FD interface #\n";
        printf STDERR "#########################################################\n";
        mkdir "GRT", 0755 if !-e "GRT";
        system "$BIN/check_grt par=$PARFILE";
    }

    printf STDERR "####################################################\n";
    printf STDERR "# Demult the wavefield from GRT to input of FD     #\n";
    printf STDERR "####################################################\n";
    system ("$BIN/demult par=$PARFILE");
}

if($FDrun eq "yes"){
    my $greenfile;
    my $diffile =  &findpar($PARFILE, "diffile");
    if($diffile == 1){
       $greenfile = &findpar($PARFILE, "greenfile1");
    } else {
       $greenfile = &findpar($PARFILE, "greenfile");
    }
    print STDERR "greenfile $greenfile\n";
    my $readpars =  &findpar($PARFILE, "readpars");
    my $raymodel =  &findpar($PARFILE, "raymodel");
    if($readpars != 0){
        system "$BIN/genmodel par=$PARFILE raymodel=$raymodel > ranmodel.log";
    }
    printf STDERR "####################################################\n";
    printf STDERR "# Run FD                                           #\n";
    printf STDERR "####################################################\n";
    system "$BIN/psvfd par=$PARFILE greenfile=$greenfile";
}

if ($FDdemult2kir eq 'yes') {
    printf STDERR "######################################################\n";
    printf STDERR "# demult the wavefield from FD to input of Kirchhoff #\n";
    printf STDERR "######################################################\n";
    system "$BIN/demult2kir par=$PARFILE rcore=0";

    if ($debug) {
        my $kirfile_x = &findpar($PARFILE, "kirfile_x");
        my $kirfile_z = &findpar($PARFILE, "kirfile_z");
        printf STDERR "######################################################\n";
        printf STDERR "# Debug for FD output to Kirchhoff                   #\n";
        printf STDERR "######################################################\n";
        mkdir "FDX", 0755 if !-e "FDX";
        mkdir "FDZ", 0755 if !-e "FDZ";
        system "$BIN/check_fdkir par=$PARFILE kirfile=$kirfile_x dir=FDX";
        system "$BIN/check_fdkir par=$PARFILE kirfile=$kirfile_z dir=FDZ";
    }
}

my $RAD   = 111.194924748;
if($KIRrun eq "yes"){

    my $nt_kir = &findpar($PARFILE, "nt_kir");
    my $raymodel = &findpar($PARFILE, "raymodel");
    $nt_kir *= 3;

    for(my $idist = 0; $idist < @Recivers; $idist++){
        #Reading out the Green's functions from library
        my $phase = "PKiKP";
        my $sdepth    = &SDEPTH($raymodel);

        # xcenter seems to be the distance from receiver to pierce point
        open(TAUP, "taup_pierce -ph $phase -mod prem -h $sdepth -deg $Recivers[$idist] -turn|");
        my $junk = <TAUP>;
        my ($pierce, $icbdep, $ttt) = split " ", <TAUP>;
        close(TAUP);
        my $xcenter = ($Recivers[$idist]-$pierce) * $RAD;
        print STDERR "gcarc=$Recivers[$idist]\n";

        my $xdist  = 2000;
        my $xleft  = $xcenter + 0.5*$xdist;
        my $xright = $xcenter - 0.5*$xdist;
        my $xtaper = 500;

	    my $dist  = $Recivers[$idist] * $RAD;
	    system("$BIN/Read_kirgreen par=$PARFILE dist=$dist taper=1 xcenter=$xcenter xleft=$xleft xright=$xright xtaper=$xtaper");
        if ($debug) {
            my $dir = sprintf("GRTKIR_%.2f", $Recivers[$idist]);
            mkdir $dir,0755 if !-e $dir;
            system "$BIN/check_grtkir par=$PARFILE dir=$dir";
        }

        #doing Kirchhoff
	    for(my $comp = 1; $comp < 2; $comp++){      # Only Z component are calculated!
            my ($fdkirfile, $sacfile);
	        if($comp == 0){
	            $fdkirfile = &findpar($PARFILE, "kirfile_x");
	            $sacfile   = sprintf("%.2f.accx", $Recivers[$idist]);
	        } else{
	            $fdkirfile = &findpar($PARFILE, "kirfile_z");
	            $sacfile   = sprintf("%.2f.accz", $Recivers[$idist]);
            }
            my $debugfile = sprintf("kirdebug.%.2f", $Recivers[$idist]);
	        my $tstart = `$BIN/kirch par=$PARFILE accfile=acc_junk dist=$dist fdkirfile=$fdkirfile debug=$debug debugfile=$debugfile`;
            if ($debug) {
                my $dir = sprintf("KIRCONV_%.2f", $Recivers[$idist]);
                mkdir $dir, 0755 if !-e $dir;
                system "$BIN/check_kirconv debugfile=$debugfile dir=$dir";
            }
            chomp($tstart);

	        system("$BIN/line2point par=$PARFILE accfile=acc_junk nt_kir=$nt_kir ");
	        system "mv syn.SAC.1 $sacfile";
	        &ChgSacHead($sacfile,"b",$tstart);
	        &ChgSacHead($sacfile,"gcarc",$Recivers[$idist]);

            unlink "acc_junk";
        }
    }
    my $greenfile_kir    = &findpar($PARFILE, "greenfile_kir");
    unlink $greenfile_kir;
}

sub ChgSacHead{
    $ENV{SAC_DISPLAY_COPYRIGHT}=0;
    my($file,$head,$val) = @_;
    open(SAC, "|sac");
    print SAC "r $file\n";
    print SAC "ch $head $val\n";
    print SAC "w over\n";
    print SAC "quit\n";
    close(SAC);
}

#
# subroutine for finding the value of specific par parameter from par source file
# USAGE:
#   $val = findpar("parfile","par");
# Revision History:
#         Feb 11 1997    Lianxing Wen Initial Revision
#         Feb 23 1997    Lianxing Wen add split(/ /, $c[1])
#
sub findpar{
    my ($parfile, $par) = @_;
    my $par_val;
    die "No $parfile found!\n" unless -f $parfile;

    open(PT,$parfile);
    foreach (<PT>) {
        chomp($_);
        my @c = split('=');
        next if @c==0;
        if ($c[0] eq $par) {
            my @val = split(/ /, $c[1]);
            $par_val = $val[0];
        }
    }
    close(PT);
    return ($par_val);
}

sub SDEPTH{
    # find source depth from input model file
    my($model) = @_;
    open(IN,"< $model") or die "cannot open file $model in SDEPTH";
    my $l = <IN>;
    chomp($l);
    my ($nlayer,$nb) = split(/ /,$l);

    my $sss = 0.0;
    for(my $i = 0; $i < $nb; $i++){
        $l = <IN>;
        chomp($l);
        $l =~ s/\s+/\ /g;
        $l =~ s/^\s+//;
        my($vpdum, $vsdum, $dddum, $th) = split(/ /,$l);
        $sss += $th;
    }

    close(IN);
    return $sss;
}
