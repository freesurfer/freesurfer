#!/bin/tcsh -f
gunzip -c test_data.tar.gz | tar xvf -
# ignore warning about 'cannot utime'
exit 0

