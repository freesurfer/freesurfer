#! /bin/tcsh -f 

set f = $1

# see whether there are any cerebral GM labels in this file 
# labels of interest: 3, 42

set fname = ${f:r:r}.count.txt
# set cmd = (/autofs/space/turan_001/users/lzollei/dev/mri_binarize/mri_binarize --noverbose --match 3 42 --count $fname --i $f)
# set cmd = (mri_binarize --match 3 42 9000 9001 9002 9003 9004 9005 9006 9500 9501 9502 9503 9504 9505 9506 --count $fname --i $f)
set cmd = (mri_binarize --match 3 42 --count $fname --i $f)
eval $cmd
set nums = `cat $fname`
rm -f $fname

if ($nums[1] > 0) then
  exit 1
else
  exit 0 
endif
