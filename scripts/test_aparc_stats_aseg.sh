#! /bin/bash -p
# This is a test for aparc_stats_aseg
function usage {
    echo "$0 -ref <ref_data>"
    echo " options:"
    echo "   -ref <ref_data> Reference test data"
    echo "   -r Regenerate test reference data"
}

function error_exit {
    >&2 echo "$(tput setaf 1)error:$(tput sgr 0) $@"
    exit 1
}

function cleanup {
    if [ $? -eq 0 ]; then
        if [ "$FSTEST_REGENERATE" = true ]; then
            cd $FSTEST_REGENERATION_DIR
            rm -rf "$FSTEST_REGENERATION_DIR/testdata/bert/scripts/aparc_stats_aseg.log"
            tar -czvf testdata_aparc_stats_aseg.tar.gz testdata
        else
            if [ $START_FLAG ]; then
                echo "$(tput setaf 2)success:$(tput sgr 0) test passed"
                # remove testdata directory
                rm -rf $FSTEST_TESTDATA_DIR
            fi
        fi
    else
        if [ $START_FLAG ]; then
            echo "$(tput setaf 1)error:$(tput sgr 0) test failed"
        fi
    fi
}
trap cleanup EXIT

function abort {
	error_exit "script has been killed"
	exit 1
}
trap abort SIGINT SIGTERM

if [ ! "$1" ]; then
    usage
    exit
fi

FSTEST_REGENERATE=false

while [ $# -gt 0 ]; do
    case $1 in
        -r|--regenerate)
        # regenerate the test reference data
        FSTEST_REGENERATE=true
        shift
        ;;
	-ref)
        REF_TESTDATA=$2
        shift 2
        ;;
        *)
	echo "ERROR: flag $1 unrecognized."
        exit 1
        ;;
    esac
done

FSTEST_CWD="$(pwd)"
FSTEST_SCRIPT_DIR="$(realpath $(dirname $0))"

if [ "$FSTEST_REGENERATE" = true ]; then
    echo "regenerating testdata"
    # make the temporary dir and untar the original testdata into this directory
    FSTEST_REGENERATION_DIR="${FSTEST_CWD}/testdata_regeneration"
    rm -rf $FSTEST_REGENERATION_DIR && mkdir $FSTEST_REGENERATION_DIR
    tar -xzvf "$REF_TESTDATA" -C $FSTEST_REGENERATION_DIR
    export SUBJECTS_DIR="$FSTEST_REGENERATION_DIR/testdata"
else
    START_FLAG=true
    tar -xzvf "$REF_TESTDATA" -C "$FSTEST_CWD"
    FSTEST_TESTDATA_DIR="$FSTEST_CWD/testdata"
    export SUBJECTS_DIR="$FSTEST_TESTDATA_DIR"
fi

cd $SUBJECTS_DIR

if [ "$FSTEST_REGENERATE" = false ]; then
    aparc_stats_aseg -s bert -gcs DKTatlas40

    mris_diff --maxerrs 1000 --s1 bert --s2 bert --hemi lh --aparc aparc.DKTatlas40 --aparc2 aparc.DKTatlas40_ref
    if [ $? -ne 0 ]; then
        error_exit "test failed"
    fi

    diff bert/stats/lh.aparc.DKTatlas40.stats bert/stats/lh.aparc.DKTatlas40_ref.stats -I#
    if [ $? -ne 0 ]; then
        error_exit "test failed"
    fi

    mri_diff bert/mri/aparc.DKTatlas40+aseg.mgz bert/mri/aparc.DKTatlas40_ref+aseg.mgz --debug
    if [ $? -ne 0 ]; then
        error_exit "test failed"
    fi
else
    aparc_stats_aseg -s bert -gcs DKTatlas40 -name DKTatlas40_ref
    if [ $? -ne 0 ]; then
        error_exit "generation failed"
    fi
fi
