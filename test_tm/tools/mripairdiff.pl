#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

my( $outFile, $threshold );
my( $goldPattern, $cmpPattern );
my( $start, $stop );

# ------------------------------------

# Write a welcome

my( $rcs, $rev, $rcsdate );
$rcs = '$RCSfile: mripairdiff.pl,v $';
$rev = '$Revision: 1.1 $';
$rcsdate = '$Date: 2010/12/21 20:03:19 $';

print "MRI Pair Diff\n";
print "=============\n\n";
print "Original Author: Richard G. Edgar\n";
print "$rcs\n";
print "$rev\n";
print "$rcsdate\n\n";


# -----------------------------------
# Check command line options

$outFile = "results.lua";
$threshold = $start = $stop = 0;
$goldPattern = $cmpPattern = "";

GetOptions( 'results=s' => \$outFile,
	    'threshold=f' => \$threshold,
	    'start=i' => \$start,
	    'stop=i' => \$stop,
	    'gold=s' => \$goldPattern,
	    'cmp=s' => \$cmpPattern );


# --------------------------------
# Open the output file

open( resOut, ">".$outFile ) || die "Couldn't open output file!\n";

my $dateString = `date`;
chomp $dateString;

print resOut <<"EOT";
-- -*- lua -*-
-- created: $dateString --

myTbl = {
EOT
;

# ---------------------------------



my $i;
for( $i=$start; $i<=$stop; $i++ ) {
    # Get the names
    my $goldName = sprintf( $goldPattern, $i );
    my $cmpName = sprintf( $cmpPattern, $i );
    my $diffName = sprintf( "diff%04i.mgz", $i );

    # Run the command
    my $binPath = $ENV{'TM_BIN_DIR'};
    my $cmdline = "$binPath/mri_diff --verbose $goldName $cmpName --thresh $threshold --diff $diffName";

    system( $cmdline );
    if ($? == -1) {
	die "Failed to execute: $!\n";
    } elsif ($? & 127) {
	printf "Child died with signal %d, %s coredump\n",
	($? & 127), ($? & 128) ? 'with' : 'without';
	die "Exiting\n"
	}

    # Find the result
    my $resString;
    my $exitVal = ( $? >> 8 );
    if( $exitVal == 0 ) {
	$resString = "passed";
    } elsif( $exitVal == 106 ) {
	$resString = "diff";
    } else {
	$resString = "failed";
	print "Unrecognised error code from mri_diff: $exitVal\n";
    }

    # Write out the result
    print resOut << "EOT"
{
   [\"result\"] = \"$resString\",
   [\"comment\"] = \"Comparison of $goldName $cmpName\",
},
EOT
    ;
}



# Finish off the lua file
print resOut <<"EOT";
}
EOT
    ;
close( resOut );

exit;

