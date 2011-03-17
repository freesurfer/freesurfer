#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;


my( $configFile, $outFile );

# ---------------------------------

# Write a welcome

my( $rcs, $rev, $rcsdate );
$rcs = '$RCSfile: ltaMultiDiff.pl,v $';
$rev = '$Revision: 1.3 $';
$rcsdate = '$Date: 2011/03/17 18:28:35 $';

print "LTA Multi-Diff\n";
print "==============\n\n";
print "Original Author: Richard G. Edgar\n";
print "$rcs\n";
print "$rev\n";
print "$rcsdate\n\n";


# ------------------------------
# Check command line options

$configFile = "";
$outFile = "";

GetOptions( 'config=s' => \$configFile,
	    'results=s' => \$outFile );


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
    my $resString = &RunCompare( $testVar->{"source"},
				 $testVar->{"cpu"},
				 $testVar->{"gpu"},
				 $testVar->{"match"},
				 $testVar->{"diff"} );
    my $name = $testVar->{"source"};

       print resOut <<"EOT"
{
    [\"result\"] = \"$resString\",
    [\"comment\"] = \"Comparison of $name MRI\",
},
EOT
    ;
}



# --------------------------------
# Finish off the lua file
print resOut <<"EOT";
}
EOT
    ;
close( resOut );

print "\nScript complete\n";

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
	elsif( $line =~ /^(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
	    %testItem = ( "source"=>$1,
			  "cpu"=>$2,
			  "gpu"=>$3,
			  "match"=>$4,
			  "diff"=>$5 );
	    push( @tvs, {%testItem} );
	}
	
    }
    
    return( @tvs );
}



# =====================================
sub RunCompare{
    my( $source, $cpu, $gpu, $matchTol, $diffTol ) = @_;

    my $binPath = $ENV{'TM_BIN_DIR'};

    my $cmdLine = "$binPath/lta_diff $cpu $gpu";

    my $cmdOut = `$cmdLine`;

     if( $? == -1 ) {
	print "Failed to execute: $!\n";
	return( "failed" );
    } elsif( $? & 127 ) {
	printf "Child died with signal %d, %s coredump\n",
	($? & 127), ($? & 128) ? 'with' : 'without';
	return( "failed" );
    }
    my $exitVal = ( $? >> 8 );
    if( $exitVal != 0 ) {
        return( "failed" );
    }

    my $result;
    if( $cmdOut <= $matchTol ) {
	$result = "passed";
    } elsif( $cmdOut <= $diffTol ) {
	$result = "diff";
    } else {
	$result = "failed";
    }

    print "Comparison $cpu $gpu ->  $cmdOut\n";

    return( $result );
}
