#!/bin/tcsh -f

# Extract test data

gunzip -c testdata.tar.gz | tar xvf -

set threads=( 8 )

foreach num ($threads)
    
    setenv OMP_NUM_THREADS $num
    echo "OMP_NUM_THREADS: " $OMP_NUM_THREADS

    time ./run_one.tcsh

    time ./run_two.tcsh
end



#
# cleanup
#
cd ..
rm -Rf testdata
