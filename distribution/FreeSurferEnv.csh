#############################################################################
# Name:    FreeSurferEnv.csh
# Purpose: Setup the environment to run FreeSurfer/FS-FAST (and FSL)
# Usage:   See help section below
# Note:    The bash equivalent script is FreeSurferEnv.sh, and should
#          be maintained to operate the same way.
#
#############################################################################

set VERSION = 'FreeSurferEnv.csh @FS_VERSION@';

## Print help if --help or -help is specified
if (("$1" == "--help") || ("$1" == "-help")) then
    echo "FreeSurferEnv.csh"
    echo ""
    echo "Purpose: Setup the environment to run FreeSurfer and FS-FAST"
    echo ""
    echo "Usage:"
    echo ""
    echo "1. Create an environment variable called FREESURFER_HOME and"
    echo "   set it to the directory in which FreeSurfer is installed."
    echo "2. From a csh or tcsh shell or (.login): "
    echo '       source $FREESURFER_HOME/FreeSurferEnv.csh'
    echo "3. There are environment variables that should point to locations"
    echo "   of software or data used by FreeSurfer. If set prior to"
    echo "   sourcing, they will not be changed, but will otherwise be"
    echo "   set to default locations:"
    echo "       FSFAST_HOME"
    echo "       SUBJECTS_DIR"
    echo "       MINC_BIN_DIR"
    echo "       MINC_LIB_DIR"
    echo "       FSL_DIR"
    echo "4. If NO_MINC is set (to anything), "
    echo "   then all the MINC stuff is ignored."
    echo "5. If NO_FSFAST is set (to anything), "
    echo "   then the startup.m stuff is ignored."
    echo "6. The script will print the final settings for the above "
    echo "   variables as well as any warnings about missing directories."
    echo "   If FS_FREESURFERENV_NO_OUTPUT is set, then no normal output"
    echo "   will be made (only error messages)."
    echo ""
    echo "The most convenient way to use this script is to write another"
    echo "script that sets FREESURFER_HOME and possibly SUBJECTS_DIR for"
    echo "your set-up, as well as NO_MINC, NO_FSFAST, or"
    echo "FS_FREESURFERENV_NO_OUTPUT as appropriate, and then source this"
    echo "script.  See SetUpFreeSurfer.csh for an example."
    exit 0;
endif

## Get the name of the operating system
set os = `uname -s`
setenv OS $os

## Set this environment variable to suppress the output.
if( $?FS_FREESURFERENV_NO_OUTPUT ) then
    set output = 0
else
    set output = 1
endif
if($?USER == 0 || $?prompt == 0) then
    set output = 0
endif

## Check if FREESURFER_HOME variable exists, then check if the actual
## directory exists.
if(! $?FREESURFER_HOME) then
    echo "ERROR: environment variable FREESURFER_HOME is not defined"
    echo "       Run the command 'setenv FREESURFER_HOME <FreeSurferHome>'"
    echo "       where <FreeSurferHome> is the directory where FreeSurfer"
    echo "       is installed."
    exit 1;
endif

if(! -e $FREESURFER_HOME) then
    echo "ERROR: $FREESURFER_HOME "
    echo "       does not exist. Check that this value is correct.";
    exit 1;
endif

if( $output ) then
    if (-e $FREESURFER_HOME/build-stamp.txt) then
        echo "-------- `cat $FREESURFER_HOME/build-stamp.txt` --------"
    endif
    echo "Setting up environment for FreeSurfer/FS-FAST (and FSL)"
    if (("$1" == "--version") || \
        ("$1" == "--V") || \
        ("$1" == "-V") || \
        ("$1" == "-v")) then
        echo $VERSION
    endif
endif

## Now we'll set directory locations based on FREESURFER_HOME for use
## by other programs and scripts.

## Set up the path. They should probably already have one, but set a
## basic one just in case they don't. Then add one with all the
## directories we just set.  Additions are made along the way in this
## script.
if(! $?path ) then
    set path = ( ~/bin /bin /usr/bin /usr/local/bin )
endif

## If FS_OVERRIDE is set, automatically assign the standard freesurfer
## defaults to all locations.  Otherwise, it will only do so if the
## variable isn't already set.

if(! $?FS_OVERRIDE) then
    setenv FS_OVERRIDE 0
endif

if(! $?FSFAST_HOME || $FS_OVERRIDE) then
    setenv FSFAST_HOME $FREESURFER_HOME/fsfast
endif

if(! $?SUBJECTS_DIR  || $FS_OVERRIDE) then
    setenv SUBJECTS_DIR $FREESURFER_HOME/subjects
endif

if(! $?FUNCTIONALS_DIR  || $FS_OVERRIDE) then
    setenv FUNCTIONALS_DIR $FREESURFER_HOME/sessions
endif

if((! $?NO_MINC) && (! $?MINC_BIN_DIR  || $FS_OVERRIDE)) then
    # try to find minc toolkit binaries
    if ( $?MNI_INSTALL_DIR) then
        setenv MINC_BIN_DIR $MNI_INSTALL_DIR/bin
        setenv MNI_DIR $MNI_INSTALL_DIR
    else if ( -e $FREESURFER_HOME/mni/bin) then
        setenv MINC_BIN_DIR $FREESURFER_HOME/mni/bin
        setenv MNI_DIR $FREESURFER_HOME/mni
    else if ( -e /usr/pubsw/packages/mni/current/bin) then
        setenv MINC_BIN_DIR /usr/pubsw/packages/mni/current/bin
        setenv MNI_DIR /usr/pubsw/packages/mni/current
    else if ( -e /usr/local/pubsw/packages/mni/current/bin) then
        setenv MINC_BIN_DIR /usr/local/pubsw/packages/mni/current/bin
        setenv MNI_DIR /usr/local/pubsw/packages/mni/current
    else if ( -e /usr/local/mni/bin) then
        setenv MINC_BIN_DIR /usr/local/mni/bin
        setenv MNI_DIR /usr/local/mni
    endif
endif
if((! $?NO_MINC) && (! $?MINC_LIB_DIR  || $FS_OVERRIDE)) then
    # try to find minc toolkit libraries
    if ( $?MNI_INSTALL_DIR) then
        setenv MINC_LIB_DIR $MNI_INSTALL_DIR/lib
    else if ( -e $FREESURFER_HOME/mni/lib) then
        setenv MINC_LIB_DIR $FREESURFER_HOME/mni/lib
    else if ( -e /usr/pubsw/packages/mni/current/lib) then
        setenv MINC_LIB_DIR /usr/pubsw/packages/mni/current/lib
    else if ( -e /usr/local/pubsw/packages/mni/current/lib) then
        setenv MINC_LIB_DIR /usr/local/pubsw/packages/mni/current/lib
    else if ( -e /usr/local/mni/lib) then
        setenv MINC_LIB_DIR /usr/local/mni/lib
    endif
endif
if((! $?NO_MINC) && (! $?MNI_DATAPATH  || $FS_OVERRIDE)) then
    # try to find minc toolkit data (MNI::DataDir)
    if ( $?MNI_INSTALL_DIR) then
        setenv MNI_DATAPATH $MNI_INSTALL_DIR/data
    else if ( -e $FREESURFER_HOME/mni/data) then
        setenv MNI_DATAPATH $FREESURFER_HOME/mni/data
    else if ( -e /usr/pubsw/packages/mni/current/data) then
        setenv MNI_DATAPATH /usr/pubsw/packages/mni/current/data
    else if ( -e /usr/local/pubsw/packages/mni/current/data) then
        setenv MNI_DATAPATH /usr/local/pubsw/packages/mni/current/data
    else if ( -e /usr/local/mni/data) then
        setenv MNI_DATAPATH /usr/local/mni/data
    endif
endif

if(! $?FSL_DIR  || $FS_OVERRIDE) then
    # FSLDIR is the FSL declared location, use that.
    # else try to find an installation.
    if ( $?FSLDIR ) then
        setenv FSL_DIR $FSLDIR
    else if ( -e $FREESURFER_HOME/fsl) then
        setenv FSL_DIR $FREESURFER_HOME/fsl
    else if ( -e /usr/pubsw/packages/fsl/current) then
        setenv FSL_DIR /usr/pubsw/packages/fsl/current
    else if ( -e /usr/local/pubsw/packages/fsl/current) then
        setenv FSL_DIR /usr/local/pubsw/packages/fsl/current
    else if ( -e $HOME/fsl); then
        setenv FSL_DIR $HOME/fsl
    else if ( -e /usr/local/fsl) then
        setenv FSL_DIR /usr/local/fsl
    endif
endif

setenv FREESURFER_HOME  $FREESURFER_HOME
setenv LOCAL_DIR        $FREESURFER_HOME/local

## Make sure these directories exist.
foreach d ($FSFAST_HOME $SUBJECTS_DIR)
    if(! -e $d ) then
        if( $output ) then
            echo "WARNING: $d does not exist"
        endif
    endif
end

if( ! $?FSF_OUTPUT_FORMAT) setenv FSF_OUTPUT_FORMAT nii.gz
if( $output ) then
    echo "FREESURFER_HOME   $FREESURFER_HOME"
    echo "FSFAST_HOME       $FSFAST_HOME"
    echo "FSF_OUTPUT_FORMAT $FSF_OUTPUT_FORMAT"
    echo "SUBJECTS_DIR      $SUBJECTS_DIR"
endif

if ( $output && $?TUTORIAL_DATA ) then
    if ( -d $TUTORIAL_DATA) then
    echo "TUTORIAL_DATA     $TUTORIAL_DATA"
    endif
endif

######## --------- Functional Analysis Stuff ----------- #######
if( ! $?NO_FSFAST) then
    setenv FMRI_ANALYSIS_DIR $FSFAST_HOME # backwards compatability
    set SUF = ~/matlab/startup.m
    if(! -e $SUF) then
        echo "INFO: $SUF does not exist ... creating"
        mkdir -p ~/matlab
        touch $SUF

        echo "%------------ FreeSurfer -----------------------------%" >> $SUF
        echo "fshome = getenv('FREESURFER_HOME');"                     >> $SUF
        echo "fsmatlab = sprintf('%s/matlab',fshome);"                 >> $SUF
        echo "if (exist(fsmatlab) == 7)"                               >> $SUF
        echo "    addpath(genpath(fsmatlab));"                         >> $SUF
        echo "end"                                                     >> $SUF
        echo "clear fshome fsmatlab;"                                  >> $SUF
        echo "%-----------------------------------------------------%" >> $SUF
        echo "" >> $SUF
        echo "%------------ FreeSurfer FAST ------------------------%" >> $SUF
        echo "fsfasthome = getenv('FSFAST_HOME');"                     >> $SUF
        echo "fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);"       >> $SUF
        echo "if (exist(fsfasttoolbox) == 7)"                          >> $SUF
        echo "    path(path,fsfasttoolbox);"                           >> $SUF
        echo "end"                                                     >> $SUF
        echo "clear fsfasthome fsfasttoolbox;"                         >> $SUF
        echo "%-----------------------------------------------------%" >> $SUF
    endif

    set tmp1 = `grep FSFAST_HOME $SUF       | wc -l`;
    set tmp2 = `grep FMRI_ANALYSIS_DIR $SUF | wc -l`;
    set tmp3 = `grep FREESURFER_HOME $SUF   | wc -l`;

    if($tmp1 == 0 && $tmp2 == 0 && $tmp3 == 0) then
            if( $output ) then
            echo ""
            echo "WARNING: The $SUF file does not appear to be";
            echo "         configured correctly. You may not be able"
            echo "         to run the FS-FAST programs";
            echo "Try adding the following lines to $SUF"
            echo "-----------------cut---------------------"
            echo "fshome = getenv('FREESURFER_HOME');"
            echo "fsmatlab = sprintf('%s/matlab',fshome);"
            echo "if (exist(fsmatlab) == 7)"
            echo "    addpath(genpath(fsmatlab));"
            echo "end"
            echo "clear fshome fsmatlab;"
            echo "fsfasthome = getenv('FSFAST_HOME');"
            echo "fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);"
            echo "if (exist(fsfasttoolbox) == 7)"
            echo "    path(path,fsfasttoolbox);"
            echo "end"
            echo "clear fsfasthome fsfasttoolbox;"
            echo "-----------------cut---------------------"
            echo ""
            endif
    endif
endif

### ----------- MINC Stuff -------------- ####
if( $output && $?MNI_DIR ) then
    echo "MNI_DIR           $MNI_DIR"
endif
#if( $output && $?MINC_BIN_DIR ) then
#    echo "MINC_BIN_DIR    $MINC_BIN_DIR"
#endif
#if( $output && $?MINC_LIB_DIR ) then
#    echo "MINC_LIB_DIR    $MINC_LIB_DIR"
#endif
if(! $?NO_MINC) then
    if( $?MINC_BIN_DIR) then
        if (! -d $MINC_BIN_DIR) then
            if( $output ) then
                echo "WARNING: MINC_BIN_DIR '$MINC_BIN_DIR' does not exist.";
            endif
        endif
    else
        if( $output ) then
            echo "WARNING: MINC_BIN_DIR not defined."
            echo "         'nu_correct' and other MINC tools"
            echo "         are used by some Freesurfer utilities."
            echo "         Set NO_MINC to suppress this warning."
        endif
    endif
    if( $?MINC_LIB_DIR) then
        if (! -d $MINC_LIB_DIR) then
            if( $output ) then
                echo "WARNING: MINC_LIB_DIR '$MINC_LIB_DIR' does not exist.";
            endif
        endif
    else
        if( $output ) then
            echo "WARNING: MINC_LIB_DIR not defined."
            echo "         Some Freesurfer utilities rely on the"
            echo "         MINC toolkit libraries."
            echo "         Set NO_MINC to suppress this warning."
        endif
    endif
    ## nu_correct and other MINC tools require a path to mni perl scripts
    if ((! $?MNI_PERL5LIB || $FS_OVERRIDE )) then
        if ( -e $FREESURFER_HOME/mni/share/perl5) then
            # Linux CentOS 6 w/ mni/1.5 build:
            setenv MNI_PERL5LIB       "$FREESURFER_HOME/mni/share/perl5"
        else if ( -e $MINC_LIB_DIR/perl5/5.8.8) then
            # Linux CentOS5:
            setenv MNI_PERL5LIB       "$MINC_LIB_DIR/perl5/5.8.8"
        else if ( -e $MINC_LIB_DIR/perl5/5.8.5) then
            # Linux CentOS4:
            setenv MNI_PERL5LIB       "$MINC_LIB_DIR/perl5/5.8.5"
        else if ( -e $MINC_LIB_DIR/perl5/5.8.3) then
            # Linux FC2:
            setenv MNI_PERL5LIB       "$MINC_LIB_DIR/perl5/5.8.3"
        else if ( -e $MINC_LIB_DIR/perl5/site_perl/5.8.3) then
            # Linux:
            setenv MNI_PERL5LIB       "$MINC_LIB_DIR/perl5/site_perl/5.8.3"
        else if ( -e $MINC_LIB_DIR/perl5/5.8.0) then
            # Linux RH9:
            setenv MNI_PERL5LIB       "$MINC_LIB_DIR/perl5/5.8.0"
        else if ( -e $MINC_LIB_DIR/5.6.0) then
            # Linux RH7 and RH9:
            setenv MNI_PERL5LIB       "$MINC_LIB_DIR/5.6.0"
        else if ( -e $MINC_LIB_DIR/../Library/Perl/Updates/5.12.3 ) then
            # Max OS X Mavericks build default:
            setenv MNI_PERL5LIB "$MINC_LIB_DIR/../Library/Perl/Updates/5.12.3"
        else if ( -e $MINC_LIB_DIR/../Library/Perl/Updates/5.10.0 ) then
            # Max OS X Snow Leopard default:
            setenv MNI_PERL5LIB "$MINC_LIB_DIR/../Library/Perl/Updates/5.10.0"
        else if ( -e $MINC_LIB_DIR/../System/Library/Perl/5.8.6 ) then
            # Max OS X Tiger default:
            setenv MNI_PERL5LIB "$MINC_LIB_DIR/../System/Library/Perl/5.8.6"
        else if ( -e $MINC_LIB_DIR/../System/Library/Perl/5.8.1 ) then
            # Max OS X Panther default:
            setenv MNI_PERL5LIB "$MINC_LIB_DIR/../System/Library/Perl/5.8.1"
        else if ( -e $MINC_LIB_DIR/MNI) then
            # Solaris:
            setenv MNI_PERL5LIB       "$MINC_LIB_DIR"
        else if ( -e $MINC_LIB_DIR/perl5/5.8) then
            # Cygwin:
            setenv MNI_PERL5LIB       "$MINC_LIB_DIR/perl5/5.8"
        else
            setenv MNI_PERL5LIB       ""
        endif
    endif
    if ((! $?PERL5LIB || $FS_OVERRIDE )) then
        setenv PERL5LIB       $MNI_PERL5LIB
    else if ( "$PERL5LIB" != "$MNI_PERL5LIB" ) then
        setenv PERL5LIB      "$MNI_PERL5LIB":"$PERL5LIB"
    endif
#    if( $output && $?PERL5LIB ) then
#        echo "PERL5LIB        $PERL5LIB"
#    endif
endif
if(! $?NO_MINC) then
    if ( $?MINC_BIN_DIR) then
        set path = ( $MINC_BIN_DIR $path )
    endif
endif


### - Miscellaneous support libraries (tiff/jpg/glut - Mac OS only) - ###
if ( -e $FREESURFER_HOME/lib/misc/bin ) then
    set path = ( $FREESURFER_HOME/lib/misc/bin \
                 $path \
                )
endif
if ( -e $FREESURFER_HOME/lib/misc/lib ) then
    setenv MISC_LIB  $FREESURFER_HOME/lib/misc/lib
    if(! $?LD_LIBRARY_PATH ) then
        setenv LD_LIBRARY_PATH $MISC_LIB
    else
        setenv LD_LIBRARY_PATH "$MISC_LIB":"$LD_LIBRARY_PATH"
    endif
    if(! $?DYLD_LIBRARY_PATH ) then
        setenv DYLD_LIBRARY_PATH $MISC_LIB
    else
        setenv DYLD_LIBRARY_PATH "$MISC_LIB":"$DYLD_LIBRARY_PATH"
    endif
endif
if( $output && $?MISC_LIB ) then
    echo "MISC_LIB        $MISC_LIB"
endif

### - freeview binary should be in the path - Mac OS only - ### 
if ( -e $FREESURFER_HOME/bin/freeview.app ) then
    set path = ( $FREESURFER_HOME/bin/freeview.app/Contents/MacOS/ \
                 $path \
                )
endif

### - Add gcc libs to DYLD_LIBRARY_PATH - MacOS only - ###
if ( -e $FREESURFER_HOME/lib/gcc/lib ) then
    if(! $?DYLD_LIBRARY_PATH ) then
        setenv DYLD_LIBRARY_PATH $FREESURFER_HOME/lib/gcc/lib
    else
        setenv DYLD_LIBRARY_PATH "$FREESURFER_HOME/lib/gcc/lib":"$DYLD_LIBRARY_PATH"
    endif

endif

### ----------- VXL (shared lib support)  ------------ ####
if ( $?VXL_LIB ) then
    if(! $?LD_LIBRARY_PATH ) then
        setenv LD_LIBRARY_PATH $VXL_LIB
    else
        setenv LD_LIBRARY_PATH "$VXL_LIB":"$LD_LIBRARY_PATH"
    endif
    if(! $?DYLD_LIBRARY_PATH ) then
        setenv DYLD_LIBRARY_PATH $VXL_LIB
    else
        setenv DYLD_LIBRARY_PATH "$VXL_LIB":"$DYLD_LIBRARY_PATH"
    endif
endif
if( $output && $?VXL_LIB ) then
    echo "VXL_LIB         $VXL_LIB"
endif


### ----------- FSL ------------ ####
if ( $?FSL_DIR ) then
    setenv FSLDIR $FSL_DIR
    # FSL >= 6.0.6
    if ( -d $FSL_DIR/share/fsl/bin) then
        setenv FSL_BIN $FSL_DIR/share/fsl/bin
    # FSL <= 6.0.5.2
    else
        setenv FSL_BIN $FSL_DIR/bin
    endif
    if(! -d $FSL_BIN) then
        if( $output ) then
            echo "WARNING: $FSL_BIN does not exist.";
        endif
    endif
    if ( -e ${FSL_DIR}/etc/fslconf/fsl.csh ) then
        source ${FSL_DIR}/etc/fslconf/fsl.csh
    endif
    setenv FSLOUTPUTTYPE NIFTI_GZ
    # use local ImageMagick stuff
    if ( -e /usr/bin/display) setenv FSLDISPLAY /usr/bin/display
    if ( -e /usr/bin/convert) setenv FSLCONVERT /usr/bin/convert
endif
if ( $?FSL_BIN ) then
    # set path = ( $FSL_BIN $path )
    # avoid PATH conflicts with what might be in FSL_BIN for newer FSL releases
    set path = ( $path $FSL_BIN )
endif
if( $output && $?FSL_DIR ) then
    echo "FSL_DIR           $FSL_DIR"
endif


### ----------- Freesurfer Bin and Lib Paths  ------------ ####
if ( -e $FREESURFER_HOME/tktools ) then
    # tktools dir could be deleted to remove Cortech license dependency
    set path = ( $FREESURFER_HOME/tktools $path )
endif
set path = ( $FREESURFER_HOME/bin \
             $FSFAST_HOME/bin \
             $path \
            )

# set FREESURFER to match FREESURFER_HOME
setenv FREESURFER $FREESURFER_HOME

# This turns on "fixing" of group surface area. A group subject made
# with make_average_subject will have a surface area smaller than
# the average of the subjects. This makes it appear to have a surface
# area the same as that of the average of the input subjects. This
# affects surface smoothing (fwhm->niterations), computation of cluster 
# sizes, and cluster-wise correction for multiple comparisons.
# FIX_VERTEX_AREA does not need to be set to anything in particular.
# To turn off, the env variable should not exist at all.
setenv FIX_VERTEX_AREA 


# cause OS to build new bin path cache:
rehash;

exit 0;
####################################################################
