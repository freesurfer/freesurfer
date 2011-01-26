#!/bin/bash

# Simple script to be run after an 'hg push'
sleep 5
/usr/bin/curl http://avebury-vm.nmr.mgh.harvard.edu/redmine/projects/fsgpu/repository -o /dev/null >& /dev/null

