#!/bin/bash

HELP=" 
  Checks the status of subjects being processed by recon-all
  in the SUBJECTS_DIR, or in the directory passed as an argument.

  Must specify zero or one arguments.

  USAGE:

  check_recons.sh     (Uses SUBJECTS_DIR)
  check_recons.sh     <subject_directory>
"

if [[ $# -eq 0 ]]; then
  [ -z "${SUBJECTS_DIR}" ] && echo "Need to set SUBJECTS_DIR" && exit 1;
elif [[ $# -eq 1 ]]; then
  if [[ "${1}" == \-*help ]]; then
    echo "${HELP}"
    exit 1
  else
    SUBJECTS_DIR="${1}"
  fi
else
  echo "${HELP}"
  exit 1
fi
[ ! -d "${SUBJECTS_DIR}" ] && echo "${SUBJECTS_DIR} not a valid directory." && exit 1;
echo "SUBJECTS_DIR=${SUBJECTS_DIR}"

INACTIVE_LIMIT=60
COMPLETED=(); ERRORED=(); IS_RUNNING=(); IS_INACTIVE=()
for d in `find ${SUBJECTS_DIR}/ -maxdepth 1 -type d`; do
  ##echo "d = ${d}"
  if [[ -e ${d}/scripts/recon-all.done ]] && [[ ! -e ${d}/scripts/recon-all.error ]]; then
    COMPLETED+=(`basename $d`)
  elif [[ -e ${d}/scripts/recon-all.done ]] && [[ -e ${d}/scripts/recon-all.error ]]; then
    ERRORED+=(`basename $d`)
  elif [[ -e ${d}/scripts/IsRunning.lh+rh ]]; then
    find ${d} -type f -mmin -${INACTIVE_LIMIT} | egrep '.*' > /dev/null 2>&1
    if [[ $? -eq 0 ]]; then 
      IS_RUNNING+=(`basename ${d}`)
    else
      IS_INACTIVE+=(`basename ${d}`)
    fi
  fi
done

echo "(${#COMPLETED[@]}) Subjects completed SUCCESSFULLY:"
IFS=$'\n' sorted=($(sort <<<"${COMPLETED[*]}"))
unset IFS
printf '%s\n' "${sorted[@]}"
if [[ ${#sorted[@]} -ne 0 ]]; then echo "";fi

echo "(${#ERRORED[@]}) Subjects completed with ERRORS:"
IFS=$'\n' sorted=($(sort <<<"${ERRORED[*]}"))
unset IFS
printf '%s\n' "${sorted[@]}"
if [[ ${#sorted[@]} -ne 0 ]]; then echo "";fi

echo "(${#IS_RUNNING[@]}) Subjects STILL RUNNING:"
IFS=$'\n' sorted=($(sort <<<"${IS_RUNNING[*]}"))
unset IFS
printf '%s\n' "${sorted[@]}"
if [[ ${#sorted[@]} -ne 0 ]]; then echo "";fi

echo "(${#IS_INACTIVE[@]}) Subjects INACTIVE for more than ${INACTIVE_LIMIT} minutes and may have not exited gracefully:"
IFS=$'\n' sorted=($(sort <<<"${IS_INACTIVE[*]}"))
unset IFS
printf '%s\n' "${sorted[@]}"
if [[ ${#sorted[@]} -ne 0 ]]; then echo "";fi
