#! /bin/tcsh -f

# Find k neighbors by MI after registration

# run mi avlauation so that that file already exists

# echo FIND NEIGHBORS MI ARG LIST $argv

source subjects.csh
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

# echo DONE WITH INPUT

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

# echo DONE WITH CHECK

###### gather all MI scores for the input subject
set voltype = norm

set miscores = ()
foreach target ($SUBJECTS)
  if !($target == $subj) then
    set DRAMMSmoved = $DRAMMSOUTDIR/$subj/$voltype-${subj}-2-${target}.DRAMMS.nii.gz
    set tmpval = `cat ${DRAMMSmoved:r:r}.finalmival.txt`
    set miscores = ($miscores $tmpval)
  else
    set miscores = ($miscores -1000.00) 
  endif
end

# echo MISCORES $miscores

### sort the scores -- ascending
### and get rid of all commas -- pretty ugly!
set sortedmiscores = `echo "$miscores" | tr ', ' '\n' | sort -k1,1n | paste -s -d' ' -`

# echo SORTEDMISCORES $sortedmiscores

### now find the corresponsing subjects!    
set trainingsubjects = ()
set counter = 0
while ($counter < $k)
  set ind = `echo "$N-$counter" | bc`; #echo $ind
  set currscore = $sortedmiscores[$ind]
  set innercounter = 1
  foreach m ($miscores)
    #echo $m $currscore
    if (`echo "$m == $currscore" | bc`) then
      set trainingsubjects = ($trainingsubjects $SUBJECTS[$innercounter])
      break
    endif
    @ innercounter = $innercounter + 1
  end
  @ counter = $counter + 1 
end

#echo The list of training subjects ${trainingsubjects}
echo ${trainingsubjects}
if ($writeout) then 
  echo ${trainingsubjects} > $outfile
endif
