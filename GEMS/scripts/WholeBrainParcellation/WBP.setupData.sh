#!/bin/tcsh -ef

#
# Set up the directory structure of the Bucker 39 data.
# Note that this needs to be done in CSH rather than BASH
# because the script containing which subjects to use is
# written in CSH...
#

source /autofs/space/dijon_032/users/salat/RBnew/scripts/subjects.csh
#echo $SUBJECTS
set subjectNumber = 0;
foreach subjectID ( $SUBJECTS )
  @ subjectNumber = $subjectNumber + 1;
 
  #echo $subjectID
  #echo $subjectNumber

  # Make directory named subjectXX
  set directoryName = 'subject';
  if ( $subjectNumber < 10 ) then
    set directoryName = $directoryName'0'
  endif
  set directoryName = $directoryName$subjectNumber
  echo $directoryName
  mkdir $directoryName

  # Get inside
  cd $directoryName

  # Put symbolic links
  ln -s /autofs/space/dijon_032/users/salat/RBnew/$subjectID/mri/orig.mgz .
  ln -s /autofs/space/dijon_032/users/salat/RBnew/$subjectID/mri/seg_edited.mgz .

  # Go back
  cd ..
end

