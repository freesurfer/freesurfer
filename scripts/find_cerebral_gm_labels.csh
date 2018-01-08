#! /bin/tcsh -f 

set f = $1

# see whether there are any cerebral GM labels in this file 
# labels of interest: 3, 42

set fname = ${f:r:r}.count.txt
set cmd = (mri_binarize --match 3 42 --count $fname --i $f)
eval $cmd
set nums = `cat $fname`
rm -f $fname

if ($nums[1] > 0) then
  exit 1
else
  exit 0 
endif
