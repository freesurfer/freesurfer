#!/bin/bash -p

#############################################################################
# Name:    FreeSurferEnv.sh
# Purpose: Setup the environment to run FreeSurfer/FS-FAST (and FSL)
# Usage:   See help section below
# Note:    The csh/tcsh equivalent script is FreeSurferEnv.csh, and should
#          be maintained to operate the same way.
#
#############################################################################

VERSION='FreeSurferEnv.sh @FS_VERSION@'

## Print help if --help or -help is specified
if [ $# -gt 1 ]; then
  if [[ "$1" == "--help" || "$1" == "-help" ]]; then
    echo "FreeSurferEnv.sh"
    echo ""
    echo "Purpose: Setup the environment to run FreeSurfer and FS-FAST"
    echo ""
    echo "Usage:"
    echo ""
    echo "1. Create an environment variable called FREESURFER_HOME and"
    echo "   set it to the directory in which FreeSurfer is installed."
    echo "2. From a sh or bash shell or (.bash_login): "
    echo '       source $FREESURFER_HOME/FreeSurferEnv.sh'
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
    echo "script.  See SetUpFreeSurfer.sh for an example."
    return 0;
  fi
fi

## Get the name of the operating system
os=`uname -s`
export OS=$os

## Set this environment variable to suppress the output.
if [ -n "$FS_FREESURFERENV_NO_OUTPUT" ]; then
    output=0
else
    output=1
fi

if [[ -z "$USER" || -z "$PS1" ]]; then
    output=0
fi

## Check if FREESURFER_HOME variable exists, then check if the actual
## directory exists.
if [ -z "$FREESURFER_HOME" ]; then
    echo "ERROR: environment variable FREESURFER_HOME is not defined"
    echo "       Run the command 'export FREESURFER_HOME <FreeSurferHome>'"
    echo "       where <FreeSurferHome> is the directory where FreeSurfer"
    echo "       is installed."
    return 1;
fi

if [ ! -d $FREESURFER_HOME ]; then
    echo "ERROR: $FREESURFER_HOME "
    echo "       does not exist. Check that this value is correct.";
    return 1;
fi

if [[ $output == 1 ]]; then
    if [ -e $FREESURFER_HOME/build-stamp.txt ]; then
        echo "-------- `cat $FREESURFER_HOME/build-stamp.txt` --------"
    fi
    echo "Setting up environment for FreeSurfer/FS-FAST (and FSL)"
    if [ $# -gt 1 ]; then
      if [[ "$1" == "--version" || \
        "$1" == "--V" || \
        "$1" == "-V" || \
        "$1" == "-v" ]]; then
        echo $VERSION
      fi
    fi
fi

## Now we'll set directory locations based on FREESURFER_HOME for use
## by other programs and scripts.

## Set up the path. They should probably already have one, but set a
## basic one just in case they don't. Then add one with all the
## directories we just set.  Additions are made along the way in this
## script.
if [ -z "$PATH" ]; then
    PATH="~/bin:/bin:/usr/bin:/usr/local/bin"
fi

## If FS_OVERRIDE is set, automatically assign the standard freesurfer
## defaults to all locations.  Otherwise, it will only do so if the
## variable isn't already set.

if [ -z "$FS_OVERRIDE" ]; then
    export FS_OVERRIDE=0
fi

if [[ -z "$FSFAST_HOME" || $FS_OVERRIDE != 0 ]]; then
    export FSFAST_HOME=$FREESURFER_HOME/fsfast
fi

if [[ -z "$SUBJECTS_DIR" || $FS_OVERRIDE != 0 ]]; then
    export SUBJECTS_DIR=$FREESURFER_HOME/subjects
fi

if [[ -z $FUNCTIONALS_DIR || $FS_OVERRIDE != 0 ]]; then
   export FUNCTIONALS_DIR=$FREESURFER_HOME/sessions
fi

if [[ -z "$NO_MINC" && ( -z "$MINC_BIN_DIR" || $FS_OVERRIDE != 0 ) ]]; then
    # try to find minc toolkit binaries
    if [ -n "$MNI_INSTALL_DIR" ]; then
        export MINC_BIN_DIR=$MNI_INSTALL_DIR/bin
        export MNI_DIR=$MNI_INSTALL_DIR
    elif [ -e $FREESURFER_HOME/mni/bin ]; then
        export MINC_BIN_DIR=$FREESURFER_HOME/mni/bin
        export MNI_DIR=$FREESURFER_HOME/mni
    elif [ -e /usr/pubsw/packages/mni/current/bin ]; then
        export MINC_BIN_DIR=/usr/pubsw/packages/mni/current/bin
        export MNI_DIR=/usr/pubsw/packages/mni/current
    elif [ -e /usr/local/pubsw/packages/mni/current/bin ]; then
        export MINC_BIN_DIR=/usr/local/pubsw/packages/mni/current/bin
        export MNI_DIR=/usr/local/pubsw/packages/mni/current
    elif [ -e /usr/local/mni/bin ]; then
        export MINC_BIN_DIR=/usr/local/mni/bin
        export MNI_DIR=/usr/local/mni
    fi
fi

if [[ -z "$NO_MINC" && ( -z "$MINC_LIB_DIR" || $FS_OVERRIDE != 0 ) ]]; then
    # try to find minc toolkit libraries
    if [ -n "$MNI_INSTALL_DIR" ]; then
        export MINC_LIB_DIR=$MNI_INSTALL_DIR/lib
    elif [ -e $FREESURFER_HOME/mni/lib ]; then
        export MINC_LIB_DIR=$FREESURFER_HOME/mni/lib
    elif [ -e /usr/pubsw/packages/mni/current/lib ]; then
        export MINC_LIB_DIR=/usr/pubsw/packages/mni/current/lib
    elif [ -e /usr/local/pubsw/packages/mni/current/lib ]; then
        export MINC_LIB_DIR=/usr/local/pubsw/packages/mni/current/lib
    elif [ -e /usr/local/mni/lib ]; then
        export MINC_LIB_DIR=/usr/local/mni/lib
    fi
fi

if [[ -z "$NO_MINC" && ( -z "$MNI_DATAPATH" || $FS_OVERRIDE != 0 ) ]]; then
    # try to find minc toolkit data (MNI::DataDir)
    if [ -n "$MNI_INSTALL_DIR" ]; then
        export MNI_DATAPATH=$MNI_INSTALL_DIR/data
    elif [ -e $FREESURFER_HOME/mni/data ]; then
        export MNI_DATAPATH=$FREESURFER_HOME/mni/data
    elif [ -e /usr/pubsw/packages/mni/current/data ]; then
        export MNI_DATAPATH=/usr/pubsw/packages/mni/current/data
    elif [ -e /usr/local/pubsw/packages/mni/current/data ]; then
        export MNI_DATAPATH=/usr/local/pubsw/packages/mni/current/data
    elif [ -e /usr/local/mni/data ]; then
        export MNI_DATAPATH=/usr/local/mni/data
    fi
fi

if [[ -z "$FSL_DIR" || $FS_OVERRIDE != 0 ]]; then
    # FSLDIR is the FSL declared location, use that.
    # else try find an installation.
    if [ -n "$FSLDIR" ]; then
        export FSL_DIR=$FSLDIR
    elif [ -e $FREESURFER_HOME/fsl ]; then
        export FSL_DIR=$FREESURFER_HOME/fsl
    elif [ -e /usr/pubsw/packages/fsl/current ]; then
        export FSL_DIR=/usr/pubsw/packages/fsl/current
    elif [ -e /usr/local/pubsw/packages/fsl/current ]; then
        export FSL_DIR=/usr/local/pubsw/packages/fsl/current
    elif [ -e $HOME/fsl ]; then
        export FSL_DIR=$HOME/fsl
    elif [ -e /usr/local/fsl ]; then
        export FSL_DIR=/usr/local/fsl
    fi
fi

export FREESURFER_HOME=$FREESURFER_HOME
export       LOCAL_DIR=$FREESURFER_HOME/local

## Make sure these directories exist.
for d in "$FSFAST_HOME" "$SUBJECTS_DIR"; do
    if [ ! -d $d ]; then
        if [[ $output == 1 ]]; then
            echo "WARNING: $d does not exist"
        fi
    fi
done

if [[ -z "$FSF_OUTPUT_FORMAT" || $FS_OVERRIDE != 0 ]]; then
  export FSF_OUTPUT_FORMAT=nii.gz
fi

if [[ $output == 1 ]]; then
    echo "FREESURFER_HOME   $FREESURFER_HOME"
    echo "FSFAST_HOME       $FSFAST_HOME"
    echo "FSF_OUTPUT_FORMAT $FSF_OUTPUT_FORMAT"
    echo "SUBJECTS_DIR      $SUBJECTS_DIR"
fi

if [[ $output == 1 && ! -z $TUTORIAL_DATA ]]; then
    if [ -d $TUTORIAL_DATA ]; then
    echo "TUTORIAL_DATA     $TUTORIAL_DATA"
    fi
fi

######## --------- Functional Analysis Stuff ----------- #######
if [[ -z "$NO_FSFAST" ]]; then
    export FMRI_ANALYSIS_DIR=$FSFAST_HOME # backwards compatibility
    SUF=~/matlab/startup.m
    if [ ! -e $SUF ]; then
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
    fi

    tmp1=`grep FSFAST_HOME $SUF       | wc -l`;
    tmp2=`grep FMRI_ANALYSIS_DIR $SUF | wc -l`;
    tmp3=`grep FREESURFER_HOME $SUF   | wc -l`;

    if [[ $tmp1 == 0 && $tmp2 == 0 && $tmp3 == 0 ]]; then
        if [[ $output == 1 ]]; then
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
        fi
    fi
fi

### ----------- MINC Stuff -------------- ####
if [[ $output == 1 && -n "$MNI_DIR" ]]; then
    echo "MNI_DIR           $MNI_DIR"
fi
#if [[ $output == 1 && -n "$MINC_BIN_DIR" ]]; then
#    echo "MINC_BIN_DIR    $MINC_BIN_DIR"
#fi
#if [[ $output == 1 && -n "$MINC_LIB_DIR" ]]; then
#    echo "MINC_LIB_DIR    $MINC_LIB_DIR"
#fi
if [ -z "$NO_MINC" ]; then
    if [ -n "$MINC_BIN_DIR" ]; then
        if [ ! -d $MINC_BIN_DIR ]; then
            if [[ $output == 1 ]]; then
                echo "WARNING: MINC_BIN_DIR '$MINC_BIN_DIR' does not exist.";
            fi
        fi
    else
        if [[ $output == 1 ]]; then
            echo "WARNING: MINC_BIN_DIR not defined."
            echo "         'nu_correct' and other MINC tools"
            echo "         are used by some Freesurfer utilities."
            echo "         Set NO_MINC to suppress this warning."
        fi
    fi
    if [ -n "$MINC_LIB_DIR" ]; then
        if [ ! -d $MINC_LIB_DIR ]; then
            if [[ $output == 1 ]]; then
                echo "WARNING: MINC_LIB_DIR '$MINC_LIB_DIR' does not exist.";
            fi
        fi
    else
        if [[ $output == 1 ]]; then
            echo "WARNING: MINC_LIB_DIR not defined."
            echo "         Some Freesurfer utilities rely on the"
            echo "         MINC toolkit libraries."
            echo "         Set NO_MINC to suppress this warning."
        fi
    fi
    ## nu_correct and other MINC tools require a path to perl
    if [[ -z "$MNI_PERL5LIB" || $FS_OVERRIDE != 0 ]]; then
        if [ -e $FREESURFER_HOME/mni/share/perl5 ]; then
            # Linux CentOS6 w/ mni 1.5 build
            export MNI_PERL5LIB="$FREESURFER_HOME/mni/share/perl5"
        elif [ -e $MINC_LIB_DIR/perl5/5.8.8 ]; then
            # Linux CentOS5:
            export MNI_PERL5LIB="$MINC_LIB_DIR/perl5/5.8.8"
        elif [ -e $MINC_LIB_DIR/../Library/Perl/Updates/5.12.3 ]; then
            # Max OS X Mavericks default:
            export MNI_PERL5LIB="$MINC_LIB_DIR/../Library/Perl/Updates/5.12.3"
        elif [ -e $MINC_LIB_DIR/../Library/Perl/Updates/5.10.0 ]; then
            # Max OS X Snow Leopard default:
            export MNI_PERL5LIB="$MINC_LIB_DIR/../Library/Perl/Updates/5.10.0"
        elif [ -e $MINC_LIB_DIR/../System/Library/Perl/5.8.6 ]; then
            # Max OS X Tiger default:
            export MNI_PERL5LIB="$MINC_LIB_DIR/../System/Library/Perl/5.8.6"
        elif [ -e $MINC_LIB_DIR/../System/Library/Perl/5.8.1 ]; then
            # Max OS X Panther default:
            export MNI_PERL5LIB="$MINC_LIB_DIR/../System/Library/Perl/5.8.1"
        elif [ -e $MINC_LIB_DIR/perl5/5.8.5 ]; then
            # Linux CentOS4:
            export MNI_PERL5LIB="$MINC_LIB_DIR/perl5/5.8.5"
        elif [ -e $MINC_LIB_DIR/perl5/5.8.3 ]; then
            # Linux FC2:
            export MNI_PERL5LIB="$MINC_LIB_DIR/perl5/5.8.3"
        elif [ -e $MINC_LIB_DIR/perl5/5.8.0 ]; then
            # Linux RH9:
            export MNI_PERL5LIB="$MINC_LIB_DIR/perl5/5.8.0"
        elif [ -e $MINC_LIB_DIR/perl5/site_perl/5.8.3 ]; then
            # Linux:
            export MNI_PERL5LIB="$MINC_LIB_DIR/perl5/site_perl/5.8.3"
        elif [ -e $MINC_LIB_DIR/5.6.0 ]; then
            # Linux RH7 and RH9:
            export MNI_PERL5LIB="$MINC_LIB_DIR/5.6.0"
        elif [ -e $MINC_LIB_DIR/perl5/5.8 ]; then
            # Cygwin:
            export MNI_PERL5LIB="$MINC_LIB_DIR/perl5/5.8"
        elif [ -e $MINC_LIB_DIR/MNI ]; then
            # Solaris:
            export MNI_PERL5LIB="$MINC_LIB_DIR"
        else
            export MNI_PERL5LIB=""
        fi
    fi
    if [[ -z "$PERL5LIB" || $FS_OVERRIDE != 0 ]]; then
        export PERL5LIB="$MNI_PERL5LIB"
    elif [ "$PERL5LIB" != "$MNI_PERL5LIB" ]; then
        export PERL5LIB="$MNI_PERL5LIB":"$PERL5LIB"
    fi
    #if [[ $output == 1 && -n "$PERL5LIB" ]]; then
    #    echo "PERL5LIB        $PERL5LIB"
    #fi
fi
if [ -z "$NO_MINC" ]; then
    if [ -n "$MINC_BIN_DIR" ]; then
        PATH=$MINC_BIN_DIR:$PATH
    fi
fi


### - Miscellaneous support libraries TIFF/JPG/GLUT (MacOS only) - ####
if [ -d $FREESURFER_HOME/lib/misc/bin ]; then
    PATH=$FREESURFER_HOME/lib/misc/bin:$PATH
fi
if [ -d $FREESURFER_HOME/lib/misc/lib ]; then
    export MISC_LIB=$FREESURFER_HOME/lib/misc/lib
    if [ -z "$LD_LIBRARY_PATH" ]; then
        export LD_LIBRARY_PATH=$MISC_LIB
    else
        export LD_LIBRARY_PATH="$MISC_LIB":"$LD_LIBRARY_PATH"
    fi
    if [ -z "$DYLD_LIBRARY_PATH" ]; then
        export DYLD_LIBRARY_PATH=$MISC_LIB
    else
        export DYLD_LIBRARY_PATH="$MISC_LIB":"$DYLD_LIBRARY_PATH"
    fi
fi
if [[ $output == 1 && -n "$MISC_LIB" ]]; then
    echo "MISC_LIB        $MISC_LIB"
fi

### - Freeview directory should be in the path - MacOS only - ### 
if [ -d $FREESURFER_HOME/bin/freeview.app ]; then
    PATH=$FREESURFER_HOME/bin/freeview.app/Contents/MacOS/:$PATH
fi

### - Add gcc libs to DYLD_LIBRARY_PATH - MacOS only - ###
if [ -d $FREESURFER_HOME/lib/gcc/lib ]; then
    if [ -z "$DYLD_LIBRARY_PATH" ]; then
        export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib
    else
        export DYLD_LIBRARY_PATH="$FREESURFER_HOME/lib/gcc/lib":"$DYLD_LIBRARY_PATH"
    fi
fi

### ----------- VXL (shared lib support) ------------ ####
if [ -n "$VXL_LIB" ]; then
    if [ -z "$LD_LIBRARY_PATH" ]; then
        export LD_LIBRARY_PATH=$VXL_LIB
    else
        export LD_LIBRARY_PATH="$VXL_LIB":"$LD_LIBRARY_PATH"
    fi
    if [ -z "$DYLD_LIBRARY_PATH" ]; then
        export DYLD_LIBRARY_PATH=$VXL_LIB
    else
        export DYLD_LIBRARY_PATH="$VXL_LIB":"$DYLD_LIBRARY_PATH"
    fi
fi
if [[ $output == 1 && -n "$VXL_LIB" ]]; then
    echo "VXL_LIB         $VXL_LIB"
fi


### ----------- FSL ------------ ####
if [ -n "$FSL_DIR" ]; then
    export FSLDIR=$FSL_DIR
    # FSL >= 6.0.6
    if [ -d $FSL_BIN/share/fsl/bin ]; then
        export FSL_BIN=$FSL_DIR/share/fsl/bin
    # FSL <= 6.0.5.2
    else
        export FSL_BIN=$FSL_DIR/bin
    fi
    if [ ! -d $FSL_BIN ]; then
        if [[ $output == 1 ]]; then
            echo "WARNING: $FSL_BIN does not exist.";
        fi
    fi
    if [ -e $FSL_DIR/etc/fslconf/fsl.sh ]; then
        source $FSL_DIR/etc/fslconf/fsl.sh
    fi
    export FSLOUTPUTTYPE=NIFTI_GZ
    # use local ImageMagick stuff
    if [ -e /usr/bin/display ]; then
        export FSLDISPLAY=/usr/bin/display
    fi
    if [ -e /usr/bin/convert ]; then
        export FSLCONVERT=/usr/bin/convert
    fi
fi
if [ -n "$FSL_BIN" ]; then
    # PATH=$FSL_BIN:$PATH
    # avoid PATH conflicts with what might be in FSL_BIN for newer FSL releases
    PATH=$PATH:$FSL_BIN
fi
if [[ $output == 1 && -n "$FSL_DIR" ]]; then
    echo "FSL_DIR           $FSL_DIR"
fi


### ----------- Freesurfer Bin and Lib Paths  ------------ ####
if [ -d $FREESURFER_HOME/tktools ]; then
    # tktools dir could be deleted to remove Cortech license dependency
    PATH=$FREESURFER_HOME/tktools:$PATH
fi
export PATH=$FREESURFER_HOME/bin:$FSFAST_HOME/bin:$PATH

# set FREESURFER to match FREESURFER_HOME
export FREESURFER=$FREESURFER_HOME

# This turns on "fixing" of group surface area. A group subject made
# with make_average_subject will have a surface area smaller than
# the average of the subjects. This makes it appear to have a surface
# area the same as that of the average of the input subjects. This
# affects surface smoothing (fwhm->niterations), computation of cluster 
# sizes, and cluster-wise correction for multiple comparisons. 
# FIX_VERTEX_AREA does not need to be set to anything in particular.
# To turn off, the env variable should not exist at all.
export FIX_VERTEX_AREA=

# cause OS to build new bin path cache:
#rehash;  # not needed for bash!

return 0;
####################################################################
