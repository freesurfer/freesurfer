#!/usr/bin/env bash
readonly PROGRAM_DIRECTORY=$(cd "$(dirname "$0")"; pwd)
MY_PYTHON="${PROGRAM_DIRECTORY}/python3"
RUN_SAMSEG="${PROGRAM_DIRECTORY}/run_samseg"
${MY_PYTHON} ${RUN_SAMSEG} $@
