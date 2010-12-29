#! /usr/bin/perl -w

#use strict;

use lib "$ENV{'FREESURFER_HOME_BIN'}";

BEGIN {
   push @INC,"$ENV{'FREESURFER_HOME_BIN'}";
}

use Getopt::Long;
use getConfig;

############### 

$sdir = "$ENV{'SUBJECTS_DIR'}";
$tdir = "$ENV{'TEMPLATE_DIR'}";
$annotFile = "$ENV{'ANNOTFILE'}";
$volType = "$ENV{'VOLTYPE'}";
$hemi = "$ENV{'HEMI'}";

print " =====================\n settings: $sdir      \n========================\n";
print " =====================\n settings: $annotFile \n========================\n";
print " =====================\n settings: $volType   \n========================\n";

###############

my $exeFile = "surf2vol";
my $vol = "";
my $refVol = "";
my $settingsFile = "";
my $dbgOut = 0;
my $usePial = 0;
my $useMpi = 0;

GetOptions ( "moving=s" => \$vol,
	     "fixed=s" => \$refVol,
	     "settings=s" => \$settingsFile,
	     "pial" => \$usePial,
	     "dbgOut" => \$dbgOut,
	     "exe=s" => \$exeFile,
	     "mpi" => \$useMpi,
	   );

if ( length $vol == 0 or
     length $refVol == 0 or
     length $settingsFile == 0 )
  {
    print " Please provide moving, fixed and settins arguments\n Optionally, dbgOut will write intermediate morph files\n";
    print " Other options \n" .
      " exe <s> - executable to use\n" .
	" pial - use pial surfaces\n" .
	  " cpus <n> - number of CPUs to use\n" .
	    " no_mpi - option NOT to use MPI\n";
    exit 1;B
  }


#------------------------------------------------------------
# STANDARD freesurfer convention
$volData    = "$sdir/$vol/mri/$volType.mgz";
$refVolData = "$tdir/$refVol/mri/$volType.mgz";

$refSurf_white = "$tdir/$refVol/surf/$hemi.white";
$refSurf_pial  = "$tdir/$refVol/surf/$hemi.pial";
#$refSurf_lh_white = "$sdir/$refVol/surf/lh.white";
#$refSurf_rh_white = "$sdir/$refVol/surf/rh.white";
#$refSurf_lh_pial  = "$sdir/$refVol/surf/lh.pial";
#$refSurf_rh_pial  = "$sdir/$refVol/surf/rh.pial";

#------------------------------------------------------------

my $cmdLine;

print " =====================\n processing brain $vol\n========================\n";
print " =====================\n template brain $refVol\n========================\n";

#my $outPath = "$sdir/$vol";
## make sure out dir exists
#if (not -d $outPath)
#  {
#    mkdir $vol;
#  }

$outDir = "$ENV{'OUTDIR'}";

# NOTE: hard-coded convention for the naming convention!
$surfResample = "resample";

$surf_white = "$outDir/$hemi.$surfResample.white";
$surf_pial = "$outDir/$hemi.$surfResample.pial";

#$surf_lh_white = "$outDir/lh.$surfResample.white";
#$surf_rh_white = "$outDir/rh.$surfResample.white";
#$surf_lh_pial = "$outDir/lh.$surfResample.pial";
#$surf_rh_pial = "$outDir/rh.$surfResample.pial";

# static elastic morph

##################################
#
# elastic morph
# process settings file
#
#################################

# this will be a hash reference
my $hash = &getConfig( conf_file => "$settingsFile" );

# populate values
if ( exists $$hash{ksp_rtol} ) { $kspRtol = $$hash{ksp_rtol}; }
else { $kspRtol = 10; }

if ( exists $$hash{weight} ) { $weight = $$hash{weight}; }
else { $weight = 1; }

if ( exists $$hash{out_root} ) { $outRoot = $$hash{out_root}; }
else { $outRoot = "out_root"; }

if ( exists $$hash{surf_root} ) { $surfRoot = $$hash{surf_root}; }
else { $surfRoot = "out_surf"; }

if ( exists $$hash{overwrite} ) { $overwrite = $$hash{overwrite}; }
else { $overwrite = 0; }

if ( exists $$hash{options} ) { $otherOptions = $$hash{options}; }
else { $otherOptions = " "; }

$outElastic = "$outDir/${outRoot}_to${refVol}.mgz";
if ( (not -e "$outElastic") or $overwrite )
  {
    $cmdVols = " -fixed_mri $refVolData -moving_mri $volData";
#    $cmdAparc = " -aparc $sdir/$refVol/label/lh.$annotFile -aparc_2 $sdir/$refVol/label/rh.$annotFile";
    $cmdAparc = " -aparc $tdir/$refVol/label/$hemi.$annotFile ";
#    $cmdSurfWhite_lh = "-fixed_surf $refSurf_lh_white   -moving_surf $surf_lh_white";
#    $cmdSurfWhite_rh = "-fixed_surf_2 $refSurf_rh_white -moving_surf_2 $surf_rh_white";
    $cmdSurfWhite = "-fixed_surf $refSurf_white   -moving_surf $surf_white";

    if ( $useMpi )
      {
        $cmdOptions = "-lin_res 20 -ksp_rtol 1.0e-$kspRtol -cache_transform $outDir/transform.txt -penalty_weight $weight $mpiOtherOptions";
      }
    else
      {
        $cmdOptions = "-lin_res 20 -ksp_rtol 1.0e-$kspRtol -cache_transform $outDir/transform.txt -penalty_weight $weight $otherOptions";
      }

    $cmdOut = "-out $outElastic -out_surf $outDir/${surfRoot}_to${refVol} -out_mesh $outDir/${outRoot}_to${refVol}";
    if ( $dbgOut )
      {
	$cmdOut = $cmdOut . " -dbg_output $outDir/${outRoot}_to${refVol}-dbg ";
      }
#    $cmdSurf = " $cmdSurfWhite_lh $cmdSurfWhite_rh ";
    $cmdSurf = " $cmdSurfWhite ";
    if ( $usePial )
      {
#	$cmdSurfPial_lh = " -fixed_surf_3 $refSurf_lh_pial -moving_surf_3 $surf_lh_pial";
#	$cmdSurfPial_rh = " -fixed_surf_4 $refSurf_rh_pial -moving_surf_4 $surf_rh_pial";
	$cmdSurfPial = " -fixed_surf_2 $refSurf_pial -moving_surf_2 $surf_pial";
#	$cmdSurf = $cmdSurf . " $cmdSurfPial_lh $cmdSurfPial_rh ";
	$cmdSurf = $cmdSurf . " $cmdSurfPial ";

	# option 2 - same aparcs as for the white surface
#	$cmdAparc = $cmdAparc . " -aparc_3 $sdir/$refVol/label/lh.$annotFile -aparc_4 $sdir/$refVol/label/rh.$annotFile";
	$cmdAparc = $cmdAparc . " -aparc_2 $tdir/$refVol/label/$hemi.$annotFile ";
      }
    
    $cmdMain = "surf2vol $cmdVols $cmdSurf $cmdAparc $cmdOptions $cmdOut";
    $cmdLine = "$cmdMain &> $outDir/trace_${outRoot}_to${refVol}_$vol.txt";

    print ("trying to execute elastic registration -> cmd line = \n $cmdLine\n");
    system($cmdLine)==0 or die "error executing surf2vol";
    
  }
else 
  {
    print " output volume $outElastic exists and the noOverwrite option has been activated\n";
  }
