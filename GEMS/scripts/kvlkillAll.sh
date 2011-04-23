#!/bin/bash

# qselect prints out a job list based on specific criterions, xargs takes multi line
# input and run the command you give to it repeatedly until it has consumed the input list.


# Delete all your jobs
qselect -u $USER | xargs qdel

# Delete all your running jobs:
#qselect -u $USER -s R | xargs qdel

# Delete all your queued jobs:
#qselect -u $USER -s Q | xargs qdel
