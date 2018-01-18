#! /bin/tcsh -f

# Find k neighbors by age

echo IN FIND_NEIGHBORS_BYAGE

# source subjects.csh
source /autofs/cluster/con_001/users/lilla/BabyProjects/CNYSegmentations/scripts/subjects.csh
set N = $#AGESORTEDSUBJECTS
set writeout = 0

# read input
set cmdline = ($argv);
while( $#argv != 0 )
  set flag = $argv[1]; shift;
  switch($flag)
    case "--s":
      if ( $#argv < 1) goto arg1err;
      set subj = $argv[1]; shift;
      breaksw
    case "--k"
      if ( $#argv < 1) goto arg1err;
      set k = $argv[1]; shift;
      breaksw
    case "--outfile"
      if ( $#argv < 1) goto arg1err;
      set outfile = $argv[1]; shift;
      set writeout = 1
      breaksw
    default:
      echo ERROR: Flag $flag unrecognized.
      echo $cmdline
      exit 1
      breaksw
  endsw
end

## check input
if($#subj == 0) then
  echo "ERROR: must spec a subject id"
  exit 1;
endif
if($#k == 0) then
  set k = 2
  exit 1;
endif
if($#k >= $N) then
  echo "ERROR: k does not make sense as it is greater than then input number of subjects"
  exit 1;
endif

###
### find age of the subjects
###
set counter = 1
foreach s ($AGESORTEDSUBJECTS)
  if ($s == $subj) then
    set ind = $counter
    set subjage = $SORTEDAGES[$counter]
    echo SELECTED SUBJECT AGE $subjage
    break
  endif
  @ counter = $counter + 1
end
set mainind = $ind
###
### set list of training subjects based on nearest ages
###
set trainingsubjects = ()
set trainingages = ()
set trainingagediffs = ()
set counter = 1
@ upperind = $ind + 1
@ lowerind = $ind - 1
while ($counter <= $k)

  if($upperind > $N) then
    set upperdiff = 1000
  else
    set upper = $SORTEDAGES[$upperind]
    @ upperdiff = $upper - $subjage
  endif

  if($lowerind < 1) then
    set lowerdiff = 1000
  else
  set lower = $SORTEDAGES[$lowerind]
  @ lowerdiff =  $subjage - $lower
  endif
 
  # echo DIFFERENCES $lowerdiff  $upperdiff
  if ($lowerdiff < $upperdiff) then
    set trainingsubjects = ($trainingsubjects $AGESORTEDSUBJECTS[$lowerind] )
    set trainingages = ($trainingages $SORTEDAGES[$lowerind] )
    set trainingagediffs = ($trainingagediffs $lowerdiff )
    # set ind = $lowerind
    @ lowerind = $lowerind - 1
  else
    set trainingsubjects = ( $trainingsubjects $AGESORTEDSUBJECTS[$upperind])
    set trainingages = ($trainingages $SORTEDAGES[$upperind]  )
    set trainingagediffs = ($trainingagediffs $upperdiff  )
    # set ind = $upperind
    @ upperind = $upperind + 1
  endif
  @ counter = $counter + 1
end
# echo The list of training subjects ${trainingsubjects}
echo ${trainingsubjects}
# echo The list of training ages ${trainingages}
if ($writeout) then 
  echo ${trainingsubjects} > $outfile
  echo ${trainingagediffs} >> $outfile
endif


