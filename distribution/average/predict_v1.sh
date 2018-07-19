#!/bin/bash

# predict the location of v1 for a group of subjects
# NOTE: the V1_average subject must be located (or linked in) the SUBJECTS_DIR
#
# Oliver Hinds <ohinds@mit.edu> 2007-06-16

usage() {
    echo "$0 [options] subject1 subject2 ..."
    echo " options:"
    echo "   -t <template> target image for registration "
    echo "                 (exvivo or invivo [default])"

    echo "   -i: dont use inflated surface as initial "
    echo "       registration (backward compatibility)"

    echo "   -h <hemi>  hemisphere (rh or lh) default is both hemis"

    echo "   -p: print mode (do not run commands, just print them)"
    echo "   -u, -?: print usage"
}

if [ ! "$1" ]; then
    usage
    exit
fi

## defaults

# change to exvivo if you are predicting V1 location in dead brains
regtarget=invivo

# hemisphere: lh and/or rh, default is run both
hemis="lh rh"

# switch the comment for old versions of mris_register
#inflated_flag=""
inflated_flag="-inflated"

printmode=

# parse command line arguments
#
# in: $@
parse_args() {
   # get the options
    while getopts ":s:t:h:ip?" Option
    do
	case $Option in
	    t ) regtarget="$OPTARG"
		;;
	    h ) hemis="$OPTARG"
		;;
	    i ) inflated_flag=
		;;
	    p ) printmode=1
		;;
	    * ) usage;
		exit;;
	esac
    done
    shift $(($OPTIND - 1))

    export subjs=$@
}

if [ ! "$1" ]; then
    usage;
    exit 0;
fi

# thresholds probabilities in one file to another at the optimal
# threshold of p=0.8. this threshold is based on comparisons of actual
# and predicted V1 boundary with and independent dataset (zilles
# brains)
threshold_labelfile_at_0p8() {
    infile=$1
    outfile=$2

    echo -n "thresholding the label..."
    cat $infile | sed "/^#/d;/^[0-9]*$/d;/.* 0\\.[0-7].*$/d" \
      > /tmp/tmpthreshout
    numlab=`wc -l /tmp/tmpthreshout | sed "s/ .*$//"`

    echo "# v1 atlaspredict label file auto-generated" > $outfile
    echo $numlab >> $outfile
    cat /tmp/tmpthreshout >> $outfile
    echo "done"
}

# parse command line
parse_args "$@"

# validate args
if [ "$regtarget" == "invivo" ]; then
    vivo=invivo
    template=invivo
elif [ "$regtarget" == "exvivo" ]; then
    vivo=exvivo
    template=exvivo
else 
    echo unknown target $regtarget
    usage;
    exit 1;
fi

# parm values from Hinds, et al. (2008)
if [ $vivo == "invivo" ]; then
    rhdist=10.0 
    rhparea=0.4
    lhdist=10.0 
    lhparea=1.0
elif [ $vivo == "exvivo" ]; then
    rhdist=5.0 
    rhparea=1.0
    lhdist=5.0 
    lhparea=0.1
fi

# register each subject
for subject in $subjs; do

  for hemi in $hemis; do
    if [ $hemi == "lh" ]; then
      dist=$lhdist;
      parea=$lhparea;
    else
      dist=$rhdist;
      parea=$rhparea;
    fi

    ## setup parms
    warp="-dist $dist -parea $parea"
    surf=$SUBJECTS_DIR/$subject/surf/"$hemi".sphere
    target=$SUBJECTS_DIR/V1_average/label/"$hemi".v1."$template".tif
    regfile=$SUBJECTS_DIR/$subject/surf/"$hemi".v1."$vivo".reg

    # generate curvature for inflated surface
    if [ ! -f $SUBJECTS_DIR/$subject/surf/"$hemi".inflated.H ]; then
        cmd="mris_curvature -n -a 5 -w -distances 10 10 \
            $SUBJECTS_DIR/$subject/surf/"$hemi".inflated"
        echo $cmd
        if [ ! $printmode ]; then $cmd; fi
    else
        echo "$SUBJECTS_DIR/$subject/surf/"$hemi".inflated.H already exists, skipping."; echo ""
    fi

    # register the subject to the template
    if [ ! -f $regfile ]; then
        cmd="mris_register $inflated_flag -a 4096 $warp $surf \
            $target $regfile"
        echo $cmd; echo ""
        if [ ! $printmode ]; then $cmd; fi
    else
        echo "file $regfile exists, so already been registered, skipping registration."
    fi

    # map the atlas to the subject
    atlas="$hemi".v1."$template".label
    mkdir -p "$SUBJECTS_DIR"/"$subject"/label
    reglabel="$SUBJECTS_DIR"/"$subject"/label/"$hemi".v1.prob.label
    surf=v1.$vivo.reg
    cmd="mris_spherical_average -osurf $surf -n \
        -o $subject label $atlas $hemi $surf V1_average $reglabel"
    echo $cmd; echo ""
    if [ ! $printmode ]; then $cmd; fi

    # threshold the label
    threshlabel=$SUBJECTS_DIR/$subject/label/"$hemi".v1.predict.label
    cmd="threshold_labelfile_at_0p8 $reglabel $threshlabel"
    echo $cmd; echo ""
    if [ ! $printmode ]; then $cmd; fi
  done
done
