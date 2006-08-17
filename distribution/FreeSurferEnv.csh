#############################################################################
# Name:    FreeSurferEnv.csh
# Purpose: Setup the environment to run FreeSurfer/FS-FAST (and FSL)
# Usage:   See help section below
# Note:    The bash equivalent script is FreeSurferEnv.sh, and should
#          be maintained to operate the same way.
#
# $Id: FreeSurferEnv.csh,v 1.57 2006/08/17 17:00:13 nicks Exp $
#############################################################################

set VERSION = '$Id: FreeSurferEnv.csh,v 1.57 2006/08/17 17:00:13 nicks Exp $'

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
    echo "       FUNCTIONALS_DIR"
    echo "       MINC_BIN_DIR"
    echo "       MINC_LIB_DIR"
    echo "       GSL_DIR"
    echo "       VXL_DIR"
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

## If FS_OVERRIDE is set, this script will automatically assign
## defaults to all locations. Otherwise, it will only do so if the
## variable isn't already set
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
    else if ( -e $FREESURFER_HOME/mni/bin) then
        setenv MINC_BIN_DIR $FREESURFER_HOME/mni/bin
    else if ( -e /usr/pubsw/packages/mni/current/bin) then
        setenv MINC_BIN_DIR /usr/pubsw/packages/mni/current/bin
    else if ( -e /usr/local/mni/bin) then
        setenv MINC_BIN_DIR /usr/local/mni/bin
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

if( $output ) then
    echo "FREESURFER_HOME $FREESURFER_HOME"
    echo "FSFAST_HOME     $FSFAST_HOME"
    echo "SUBJECTS_DIR    $SUBJECTS_DIR"
endif
if( $output && $?FUNCTIONALS_DIR ) then
    echo "FUNCTIONALS_DIR $FUNCTIONALS_DIR"
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
        echo "path(path,fsmatlab);"                                    >> $SUF
        echo "clear fshome fsmatlab;"                                  >> $SUF
        echo "%-----------------------------------------------------%" >> $SUF
        echo "" >> $SUF
        echo "%------------ FreeSurfer FAST ------------------------%" >> $SUF
        echo "fsfasthome = getenv('FSFAST_HOME');"                     >> $SUF
        echo "fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);"       >> $SUF
        echo "path(path,fsfasttoolbox);"                               >> $SUF
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
            echo "fsmatlab = sprintf('%s/matlab',fsmatlab);"
            echo "path(path,fsmatlab);"
            echo "clear fshome fsmatlab;"
            echo "fsfasthome = getenv('FSFAST_HOME');"
            echo "fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);"
            echo "path(path,fsfasttoolbox);"
            echo "clear fsfasthome fsfasttoolbox;"
            echo "-----------------cut---------------------"
            echo ""
            endif
    endif
endif

### ----------- MINC Stuff -------------- ####
if( $output && $?MINC_BIN_DIR ) then
    echo "MINC_BIN_DIR    $MINC_BIN_DIR"
endif
if( $output && $?MINC_LIB_DIR ) then
    echo "MINC_LIB_DIR    $MINC_LIB_DIR"
endif
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
    ## Set Load library path ##
    if(! $?LD_LIBRARY_PATH ) then
        if ( $?MINC_LIB_DIR) then
            setenv LD_LIBRARY_PATH  $MINC_LIB_DIR
        endif
    else
        if ( $?MINC_LIB_DIR) then
            setenv LD_LIBRARY_PATH "$MINC_LIB_DIR":"$LD_LIBRARY_PATH"
        endif
    endif
    ## nu_correct and other MINC tools require a path to mni perl scripts
    if (! $?MNI_PERL5LIB) then
        if ( -e $MINC_LIB_DIR/perl5/5.8.5) then
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
        else if ( -e $MINC_LIB_DIR/../System/Library/Perl/5.8.6 ) then
            # Max OS X Tiger default:
            setenv MNI_PERL5LIB "$MINC_LIB_DIR/../System/Library/Perl/5.8.6"
        else if ( -e $MINC_LIB_DIR/../System/Library/Perl/5.8.1 ) then
            # Max OS X Panther default:
            setenv MNI_PERL5LIB "$MINC_LIB_DIR/../System/Library/Perl/5.8.1"
        else if ( -e $MINC_LIB_DIR/MNI) then
            # Solaris:
            setenv MNI_PERL5LIB       "$MINC_LIB_DIR"
        else
            setenv MNI_PERL5LIB       ""
        endif
    endif
    if (! $?PERL5LIB) then
        setenv PERL5LIB       $MNI_PERL5LIB
    else if ( "$PERL5LIB" != "$MNI_PERL5LIB" ) then
        setenv PERL5LIB      "$MNI_PERL5LIB":"$PERL5LIB"
    endif
    if( $output && $?PERL5LIB ) then
        echo "PERL5LIB        $PERL5LIB"
    endif
endif
if(! $?NO_MINC) then
    if ( $?MINC_BIN_DIR) then
        set path = ( $MINC_BIN_DIR $path )
    endif
endif


### ----------- GSL (Gnu Scientific Library)  ------------ ####
if ( -e $FREESURFER_HOME/lib/gsl) then
    setenv GSL_DIR    $FREESURFER_HOME/lib/gsl
else if ( -e /usr/pubsw/packages/gsl/current) then
    setenv GSL_DIR    /usr/pubsw/packages/gsl/current
endif
if ( $?GSL_DIR ) then
    if (! $?LD_LIBRARY_PATH) then
        setenv LD_LIBRARY_PATH  $GSL_DIR/lib
    else
        setenv LD_LIBRARY_PATH  "$GSL_DIR/lib":"$LD_LIBRARY_PATH"
    endif
    if (! $?DYLD_LIBRARY_PATH) then
        setenv DYLD_LIBRARY_PATH  $GSL_DIR/lib
    else
        setenv DYLD_LIBRARY_PATH  "$GSL_DIR/lib":"$DYLD_LIBRARY_PATH"
    endif
endif
if( $output && $?GSL_DIR ) then
    echo "GSL_DIR         $GSL_DIR"
endif


### ------- Qt (scuba2 and qdec support libraries) ------- ####
# look for Qt in common NMR locations, overriding any prior setting
# NJS: QT is no longer included in the lib search path,
# as having too many files to search slows the operation of any command.
#if ( -e $FREESURFER_HOME/lib/qt) then
#    setenv QTDIR    $FREESURFER_HOME/lib/qt
#else if ( -e /usr/pubsw/packages/qt/current) then
#    setenv QTDIR    /usr/pubsw/packages/qt/current
#endif
#if ( $?QTDIR ) then
#    setenv PATH     $QTDIR/bin:$PATH
#    if (! $?LD_LIBRARY_PATH) then
#        setenv LD_LIBRARY_PATH  $QTDIR/lib
#    else
#        setenv LD_LIBRARY_PATH  "$QTDIR/lib":"$LD_LIBRARY_PATH"
#    endif
#    if (! $?DYLD_LIBRARY_PATH) then
#        setenv DYLD_LIBRARY_PATH  $QTDIR/lib
#    else
#        setenv DYLD_LIBRARY_PATH  "$QTDIR/lib":"$DYLD_LIBRARY_PATH"
#    endif
#endif
#if( $output && $?QTDIR ) then
#    echo "QTDIR           $QTDIR"
#endif


### ----------- Tcl/Tk/Tix/BLT  ------------ ####
if ( -e $FREESURFER_HOME/lib/tcltktixblt/bin ) then
    set path = ( $FREESURFER_HOME/lib/tcltktixblt/bin \
                 $path \
                )
endif
if ( -e $FREESURFER_HOME/lib/tcltktixblt/lib ) then
    setenv TCL_LIB_DIR  $FREESURFER_HOME/lib/tcltktixblt/lib
    if ( $?SET_TCL_VARS ) then
        setenv TCLLIBPATH  $TCL_LIB_DIR
        setenv TCL_LIBRARY $TCL_LIB_DIR/tcl8.4
        setenv TK_LIBRARY  $TCL_LIB_DIR/tk8.4
        setenv TIX_LIBRARY $TCL_LIB_DIR/tix8.1
        setenv BLT_LIBRARY $TCL_LIB_DIR/blt2.4
    endif
    if(! $?LD_LIBRARY_PATH ) then
        setenv LD_LIBRARY_PATH $TCL_LIB_DIR
    else
        setenv LD_LIBRARY_PATH "$TCL_LIB_DIR":"$LD_LIBRARY_PATH"
    endif
    if(! $?DYLD_LIBRARY_PATH ) then
        setenv DYLD_LIBRARY_PATH $TCL_LIB_DIR
    else
        setenv DYLD_LIBRARY_PATH "$TCL_LIB_DIR":"$DYLD_LIBRARY_PATH"
    endif
endif
if( $output && $?TCL_LIB_DIR ) then
    echo "TCL_LIB_DIR     $TCL_LIB_DIR"
endif


### -------------- VTK ------------- ###
# NJS: VTK is no longer included in the lib search path,
# as having too many files to search slows the operation of any command.
#if ( -e $FREESURFER_HOME/lib/vtk) then
#    setenv VTK_DIR    $FREESURFER_HOME/lib/vtk
#else if ( -e /usr/pubsw/packages/vtk/current) then
#    setenv VTK_DIR    /usr/pubsw/packages/vtk/current
#endif
if ( $?VTK_DIR ) then
    setenv PATH     $VTK_DIR/bin:$PATH
    if (! $?LD_LIBRARY_PATH) then
        setenv LD_LIBRARY_PATH  $VTK_DIR/lib
    else
        setenv LD_LIBRARY_PATH  "$VTK_DIR/lib":"$LD_LIBRARY_PATH"
    endif
    if (! $?DYLD_LIBRARY_PATH) then
        setenv DYLD_LIBRARY_PATH  $VTK_DIR/lib
    else
        setenv DYLD_LIBRARY_PATH  "$VTK_DIR/lib":"$DYLD_LIBRARY_PATH"
    endif
endif
if( $output && $?VTK_DIR ) then
    echo "VTK_DIR         $VTK_DIR"
endif


### -------------- VXL ------------- ###
if ( -e $FREESURFER_HOME/lib/vxl) then
    setenv VXL_DIR    $FREESURFER_HOME/lib/vxl
else if ( -e /usr/pubsw/packages/vxl/current) then
    setenv VXL_DIR    /usr/pubsw/packages/vxl/current
endif
if ( $?VXL_DIR ) then
    if (! $?LD_LIBRARY_PATH) then
        setenv LD_LIBRARY_PATH  $VXL_DIR/lib
    else
        setenv LD_LIBRARY_PATH  "$VXL_DIR/lib":"$LD_LIBRARY_PATH"
    endif
    if (! $?DYLD_LIBRARY_PATH) then
        setenv DYLD_LIBRARY_PATH  $VXL_DIR/lib
    else
        setenv DYLD_LIBRARY_PATH  "$VXL_DIR/lib":"$DYLD_LIBRARY_PATH"
    endif
endif
if( $output && $?VXL_DIR ) then
    echo "VXL_DIR         $VXL_DIR"
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


### ----------- FSL ------------ ####
if ( $?FSL_DIR ) then
    setenv FSLDIR $FSL_DIR
    setenv FSL_BIN $FSL_DIR/bin
    if(! -d $FSL_BIN) then
        if( $output ) then
            echo "WARNING: $FSL_BIN does not exist.";
        endif
    endif
    if ( -e ${FSL_DIR}/etc/fslconf/fsl.csh ) then
        source ${FSL_DIR}/etc/fslconf/fsl.csh
    endif
    setenv FSLOUTPUTTYPE NIFTI_GZ
endif
if ( $?FSL_BIN ) then
    set path = ( $FSL_BIN $path )
endif
if( $output && $?FSL_DIR ) then
    echo "FSL_DIR         $FSL_DIR"
endif


### ----------- Freesurfer Bin and Lib Paths  ------------ ####
set path = ( $FREESURFER_HOME/bin/ \
             $FSFAST_HOME/bin \
             $path \
            )

# cause OS to build new bin path cache:
rehash;

exit 0;
####################################################################
