#! /bin/tcsh -f

# Find k neighbors by MI after registration
# echo FIND NEIGHBORS MI ARG LIST $argv

##   set cmd = (find_neighbors_byMI_extended.csh --f $preprocessedfile --trpool $TEMPLATE_SUBJECTS --k $trainingsetsize --outfile $outfname --checkconflict $avoidtrainingname)
##   set cmd = (find_neighbors_byMI_extended.csh --s $s --trpool $TEMPLATE_SUBJECTS --k $trainingsetsize --outfile $outfname --checkconflict $avoidtrainingname)

set writeout = 0
set N = 0;
set trainingscores = ()
set trainingsubjects = ()
set outfile = ()
set checkconflict = 0
set fixed = ()
set k = 0

# read input
set cmdline = ($argv);
# echo CMDLINE $cmdline
while( $#argv != 0 )
  set flag = $argv[1]; shift;
  switch($flag)
    case "--f":
      if ( $#argv < 1) goto arg1err;
      set fixed = $argv[1]; shift;
      breaksw
    case "--trpool"
      if ( $#argv < 1) goto arg1err;
      while (1) 
        set trainingsubjects = ($trainingsubjects $argv[1])
        if ($#argv < 2) then
          break
        endif
        set tmp = `echo $argv[2] | sed 's/.*\(--\).*/\1/'` 
        if ($tmp == "--") then
          break
        else
          shift;
        endif
      end
      if ($#argv != 0) then 
        shift;
      endif
      # echo The list of all training subjects is: $trainingsubjects
      set N = $#trainingsubjects
      breaksw
    case "--checkconflict"  # for testing purposes -- make sure training dataset will not be equal to test dataset
      if ( $#argv < 1) goto arg1err;
      set conflictsubject = $argv[1]; shift;
      set checkconflict = 1
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
if($#fixed == 0) then
  echo "ERROR: must spec an input file"
  exit 1;
endif
if($k == 0) then
  set k = 2
  exit 1;
endif
if($k > $N) then
  echo "ERROR: k does not make sense as it is greater than input number of subjects"
  exit 1;
endif
# echo DONE WITH CHECK

## (0) identify conflict subject
if ($checkconflict) then
  set found = 0
  set counter = 1
  foreach s ($trainingsubjects)
    if ($s == $conflictsubject) then
      set found = 1
      break
    endif
    @ counter = $counter + 1
  end
  if ($found) then
    echo COUNTER $counter
    switch ($counter)
      case "1"
        echo CASE 1
        set trainingsubjects = ($trainingsubjects[2-$N])
        breaksw
      case "$N"
        echo CASE N
        @ endn = $N - 1
        set trainingsubjects = ($trainingsubjects[1-$endn])
        breaksw
      default:
        echo CASE DEF
        @ a = $counter - 1
        @ z = $counter + 1
        echo A to Z $a $z
        set trainingsubjects = ($trainingsubjects[1-$a] $trainingsubjects[$z-$N])
        echo $trainingsubjects
        breaksw
    endsw
    @ N = $N - 1
  endif
endif
## (1) register subjects and (2) get MI scores
set miscores = ()
foreach ts ($trainingsubjects)
  set mov    = $TEMPLATE_SUBJECTS_DIR/$ts/norm.nii.gz
  set out    = $WORK_DIR/subj-2-$ts.affine.mgz
  set outlta = $WORK_DIR/subj-2-$ts.affine.lta
  if !(-e $out) then
    set cmd = (mri_robust_register --mov $mov --dst $fixed --mapmov $out --lta $outlta --affine --satit) # --cost NMI 
    eval $cmd
  endif
  set cmd = (mri_mi --silent $fixed $out)
  eval $cmd
  set miscores = ($miscores `mri_mi --silent $fixed $out`)
  echo MI scores $miscores
end
## (3) sort the scores -- descending 
echo MI SCORES $miscores
set sortedmiscores = `echo $miscores | fmt -1 | sort -rn`
echo SORTEDMISCORES K $sortedmiscores $k
set selectedsortedmiscores = ($sortedmiscores[1-$k])
echo SELECTED SORTEDMISCORES $selectedsortedmiscores

### now find the corresponsing subjects!    
set selectedtrainingsubjects = ()
set counter = 1
while ($counter <= $k)
  set currscore = $sortedmiscores[$counter]
  set innercounter = 1
  foreach m ($miscores)
    if (`echo "$m == $currscore" | bc`) then
      set selectedtrainingsubjects = ($selectedtrainingsubjects $trainingsubjects[$innercounter])
      break
    endif
    @ innercounter = $innercounter + 1
  end
  @ counter = $counter + 1 
end

#echo The list of selected training subjects ${selectedtrainingsubjects}
echo ${selectedtrainingsubjects}
if ($writeout) then 
  echo ${selectedtrainingsubjects} > $outfile
endif
