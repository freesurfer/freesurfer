#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

my( $outFile );
my( $cmp1File, $cmp2File );
my( $diffThres, $matchThres );

# -----------------------------------
# Check command line options

$outFile = "";
$cmp1File = "";
$cmp2File = "";
$diffThres = 1e-4;
$matchThres = 1e-5;

GetOptions( 'results=s' => \$outFile,
            'gold=s' => \$cmp1File,
            'test=s' => \$cmp2File,
	    'match=f' => \$matchThres,
	    'diff=f' => \$diffThres );


# -------------------------------
# Get the values
my( $goldVal, $cmpVal );

$goldVal = &GetValue( $cmp1File );
$cmpVal = &GetValue( $cmp2File );

# Check for match
my $diffVal;
if( $goldVal != 0 ) {
    $diffVal = abs( 1.0 - ($cmpVal/$goldVal) );
} else {
    $diffVal = abs( $cmpVal );
}

my $resString;

print "diffVal = $diffVal\n";

if( $diffVal <= $matchThres ) {
    $resString = "passed";
} elsif ( $diffVal <= $diffThres ) {
    $resString = "diff";
} else {
    $resString = "failed";
}

# ----------------------------------
# Write output


open( resOut, ">".$outFile ) || die "Couldn't open output file!\n";

my $dateString = `date`;
chomp $dateString;

print resOut <<"EOT";
-- -*- lua -*-
-- created: $dateString --

myTbl = {
    {
	[\"result\"] = \"$resString\",
    },
}
EOT
;


exit;


# ===========================================

sub GetValue{
    my( $filename ) = @_;
    
    open( myFile, "<".$filename ) || die "Couldn't open $filename\n";
    binmode myFile;
    my $data;
    # Read 8 bytes
    read( myFile, $data, 8 );
    close( myFile );
    # Unpack into a double
    return( unpack( "d", $data ) );
}
