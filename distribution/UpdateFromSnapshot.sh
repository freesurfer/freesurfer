#! /bin/sh

ECHO=
DATE=`(set \`date +%y%m%d\`; echo $1)`

if [ "${1}" = "" ] ; then
    echo "Usage: $0 [freesurfer_directory]";
    exit 1;
fi

backup_and_update_file () {
    SRC_FILE=$1
    DEST_FILE=$2

    # If the source file doesn't exist, do nothing.
    if [ ! -e ${SRC_FILE} ] ; then
	return
    fi

    # If the dest file is the same as the source file, do nothing.
    if [ -e ${DEST_FILE} ] ; then
	diff ${SRC_FILE} ${DEST_FILE} > /dev/null
	if [ "${?}" = "0" ] ; then
	    return
	fi
    fi

    # If the dest file exists, back it up.
    if [ -e ${DEST_FILE} ] ; then  
	${ECHO} mv ${DEST_FILE} ${DEST_FILE}.${DATE};
    fi

    # Copy the source file over the dest file.
    ${ECHO} cp ${SRC_FILE} ${DEST_FILE};
    echo "${DEST_FILE} was updated"
}

process_dir () {
    local CUR=$1

    # Make the destination directory path.
    if [ "${CUR}" == "." ] ; then
	DESTDIR=${DEST}
    else
	DESTDIR=${DEST}/${CUR}
    fi

    # Get the listing for this directory.
    local CONTENTS=`/bin/ls ${CUR}`

    # For each file, look for the corresponding file in the
    # destination and if found, update it.
    for file in $CONTENTS; do
	if [ -f ${CUR}/${file} ] ; then
	    backup_and_update_file ${CUR}/${file} ${DESTDIR}/${file}
	fi
    done

    # For each directory, process it as well.
    for dir in $CONTENTS; do
	if [ -d ${CUR}/${dir} ] ; then
	    process_dir ${CUR}/${dir}
	fi
    done
}


echo "${0} ${DATE} BEGIN"

# Our base dest dir is the argument to the script. Start with the
# current dir and go down.
DEST=$1
process_dir .

echo "${0} ${DATE} COMPLETE"
