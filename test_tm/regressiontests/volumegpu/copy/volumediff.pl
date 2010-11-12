#!/usr/bin/env perl

# This is a repeat of mridiff.pl, but with extra flags
# since the volume copy doesn't do everything

use warnings;
use strict;
use Getopt::Long;

my( $outFile );
my( $cmp1File, $cmp2File );
my( $threshold );


# -----------------------------------
# Check command line options

$outFile = "";
$cmp1File = "";
$cmp2File = "";
$threshold = 0;

GetOptions( 'results=s' => \$outFile,
	    'gold=s' => \$cmp1File,
	    'test=s' => \$cmp2File,
	    'threshold=f' => \$threshold );


# ----------------------------------
# Open output file

open( resOut, ">".$outFile ) || die "Couldn't open output file!\n";

my $dateString = `date`;
chomp $dateString;

print resOut <<"EOT";
-- -*- lua -*-
-- created: $dateString --
    
EOT
;

# ----------------------------------
# Run mri_diff


my $binPath = $ENV{'TM_BIN_DIR'};
my $cmdline = "$binPath/mri_diff --notallow-res --notallow-geo --verbose $cmp1File $cmp2File --thresh $threshold --diff diff.mgz";
#print "$cmdline\n";

system( $cmdline );
if ($? == -1) {
    die "Failed to execute: $!\n";
} elsif ($? & 127) {
    printf "Child died with signal %d, %s coredump\n",
    ($? & 127), ($? & 128) ? 'with' : 'without';
    die "Exiting\n"
}

my $exitVal = ( $? >> 8 );

#print "$exitVal\n";


# ------------------------------
# Write results

my $resString;

if( $exitVal == 0 ) {
    $resString = "passed";
} elsif( $exitVal == 106 ) {
    $resString = "diff";
} else {
    $resString = "failed";
    print "Unrecognised error code from mri_diff: $exitVal\n";
}

print resOut <<"EOT"
myTbl = {
  {
    [\"result\"] = \"$resString\",
  },
}

EOT
;
close( resOut );
