#!/bin/tcsh -f

# Extract test data

gunzip -c testdata.tar.gz | tar xvf -

# Set up performance monitors (comment out to turn off)
setenv PERF "perf stat -B -e cache-references,cache-misses,faults,migrations,dTLB-loads,dTLB-load-misses"

echo "============================"

set threads=( 1 2 4 8 )

foreach num ($threads)
    
    setenv OMP_NUM_THREADS $num
    echo "OMP_NUM_THREADS: " $OMP_NUM_THREADS

    time ./run_one.tcsh

    echo "-- -- -- --"

    time ./run_two.tcsh

    echo "============================"
end



#
# cleanup
#
cd ..
rm -Rf testdata
