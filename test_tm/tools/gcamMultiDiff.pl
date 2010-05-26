#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;


my( $configFile );
my( $outFile );
my( $goldFile, $cmpFile );

# ------------------------------------

# Write a welcome

my( $rcs, $rev, $rcsdate );
$rcs = '$RCSfile: gcamMultiDiff.pl,v $';
$rev = '$Revision: 1.3 $';
$rcsdate = '$Date: 2010/05/26 19:20:00 $';

print "GCAM Multi-Diff\n";
print "===============\n\n";
print "Original Author: Richard G. Edgar\n";
print "$rcs\n";
print "$rev\n";
print "$rcsdate\n\n";


# ------------------------------
# Check command line options

$configFile = "";
$outFile = "";
$goldFile = "";
$cmpFile = "";

GetOptions( 'config=s' => \$configFile,
	    'results=s' => \$outFile,
	    'gold=s' => \$goldFile,
	    'cmp=s' => \$cmpFile );

# -------------------------------
# Load the configuration file

my( @testVars );

@testVars = &LoadConfigFile( $configFile );



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

# --------------------------------
# Run comparisons

foreach my $testVar ( @testVars ) {
    my $resString = &RunCompare( $testVar->{"field"},
				 $testVar->{"mode"},
				 $testVar->{"match"},
				 $testVar->{"diff"} );
    my $name = $testVar->{"field"};

    print resOut <<"EOT"
{
    [\"result\"] = \"$resString\",
    [\"comment\"] = \"Comparison of $name field\",
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



# =====================================

sub LoadConfigFile{

    my( $configFile ) = @_;
    my( @configLines, $line );

    my( %testItem, @tvs );

    # Slurp in the config file
    open( myConfig, "<".$configFile ) || die "Can't open $configFile\n";
    @configLines = <myConfig>;
    close( myConfig );
    
    # Loop over the lines
    foreach $line ( @configLines ) {
	if( $line =~ /^\#/ ) {
	    # Comment lines start with "#"
	}
	elsif( $line =~ /^(\w+)\s+(\w+)\s+(\S+)\s+(\S+)/ ) {
	    %testItem = ( "field"=>$1,
			  "mode"=>$2,
			  "match"=>$3,
			  "diff"=>$4 );
	    push( @tvs, {%testItem} );
	}
	
    }
    
    return( @tvs );
}


# =====================================
sub RunCompare{
    my( $varName, $mode, $matchTol, $diffTol ) = @_;

    my $binPath = $ENV{'TM_BIN_DIR'};

    my $cmdline = "$binPath/gcamCompare $varName $goldFile $cmpFile --$mode -m $matchTol -d $diffTol";

    print "$cmdline\n";

    system( $cmdline );
    if( $? == -1 ) {
	print "Failed to execute: $!\n";
	return( "failed" );
    } elsif( $? & 127 ) {
	printf "Child died with signal %d, %s coredump\n",
	($? & 127), ($? & 128) ? 'with' : 'without';
	return( "failed" );
    }

    my $exitVal = ( $? >> 8 );
    
    my $result;

    if( $exitVal == 0 ) {
	$result = "passed";
    } elsif( $exitVal == 100 ) {
	$result = "diff";
    } else {
	$result = "failed";
    }

    return( $result );
}
