#! /bin/tcsh -f

# Find k neighbors by age

# echo IN FIND_NEIGHBORS_BYAGE_EXTENDED

set writeout = 0
set N = 0
set gmwmN = 0;
set subjage = -10
set k = 0
set trainingages = ()
set trainingsubjects = ()
set gmwmtrainingages = ()
set gmwmtrainingsubjects = ()
set outfile = ()
set checkconflict = 0
set forcegmwm = 0

# read input
set cmdline = ($argv);
# echo CMDLINE $cmdline
while( $#argv != 0 )
  set flag = $argv[1]; shift;
  switch($flag)
    case "--age":
      if ( $#argv < 1) goto arg1err;
      set subjage = $argv[1]; shift;
      # echo Age of input subject: $subjage
      breaksw
    case "--k"
      if ( $#argv < 1) goto arg1err;
      set k = $argv[1]; shift;
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
    case "--gmwmtrpool"
      if ( $#argv < 1) goto arg1err;
      while (1) 
        set gmwmtrainingsubjects = ($gmwmtrainingsubjects $argv[1])
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
      # echo The list of all GMWM training subjects is: $gmwmtrainingsubjects
      set gmwmN = $#gmwmtrainingsubjects
      breaksw
    case "--trages"
      if ( $#argv < 1) goto arg1err;
      while (1)
        set trainingages = ($trainingages $argv[1])
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
      # echo The list of ages of all training subjects is: $trainingages
      breaksw
    case "--gmwmtrages"
      if ( $#argv < 1) goto arg1err;
      while (1)
        set gmwmtrainingages = ($gmwmtrainingages $argv[1])
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
      # echo The list of ages of all training subjects is: $trainingages
      breaksw
    case "--outfile"
      if ( $#argv < 1) goto arg1err;
      set outfile = $argv[1]; shift;
      # echo The output file is: $outfile
      set writeout = 1
      breaksw
    case "--checkconflict"  # for testing purposes -- make sure training dataset will not be equal to test dataset
      if ( $#argv < 1) goto arg1err;
      set conflictsubject = $argv[1]; shift;
      set checkconflict = 1
      breaksw
    case "--gmwm"
      set forcegmwm = 1;
      breaksw
    default:
      echo ERROR: Flag $flag unrecognized.
      echo $cmdline
      exit 1
      breaksw
  endsw
end

## check input
if($subjage < 0) then
  echo "ERROR: must spec a subject age"
  exit 1;
endif
if($k == 0) then
  echo Using the default neighborhood size of 2
  set k = 2
endif
if($#trainingages == 0) then
  echo "ERROR: must spec a training dataset ages" 
  exit 1;
endif
if($N == 0) then
  echo "ERROR: must spec a training dataset"
  exit 1;
endif
if($forcegmwm && (($gmwmN == 0) || ($#gmwmtrainingages == 0))) then
  echo "ERROR: must spec a GMWM training dataset and ages if at least one GMWM training data is needed"
  exit 1;
endif
if($k >= $N) then
  echo "USE the full training data set! No selection is needed"
  echo ${trainingsubjects}
  if ($writeout) then 
    echo ${trainingsubjects} > $outfile
  endif
  exit 1;
endif

##
## (0) Check for training conflict
##
if ($checkconflict) then
# echo CHECK CONFLICT
# training subjects
set found = 0
set counter = 1
foreach s ($trainingsubjects)
  if ($s == $conflictsubject) then
    set found = 1
    break
  endif
  @ counter = $counter + 1
end
# echo COUNTER $counter SUBJECT $s FOUND $found
if ($found) then
  switch ($counter)
    case "1"
      set trainingsubjects = ($trainingsubjects[2-$N])
      set trainingages = ($trainingages[2-$N])
      # echo TRAINING $trainingsubjects $trainingages
      breaksw
    case "$N"
      @ endn = $N - 1
      set trainingsubjects = ($trainingsubjects[1-$endn])
      set trainingages = ($trainingages[1-$endn])
      breaksw
    default:
      @ n1 = $counter - 1
      @ n2 = $counter + 1
      set trainingsubjects = ($trainingsubjects[1-$n1] $trainingsubjects[$n2-$N])
      set trainingages     = ($trainingages[1-$n1] $trainingages[$n2-$N])
      breaksw
  endsw
  @ N = $N - 1
endif

if ($forcegmwm) then
# echo FORCE GM
# gmwm trainign subjects
set found = 0
set counter = 1
foreach s ($gmwmtrainingsubjects)
  if ($s == $conflictsubject) then
    set found = 1
    break
  endif
  @ counter = $counter + 1
end
if ($found) then
  switch ($counter)
    case "1"
      set gmwmtrainingsubjects = ($gmwmtrainingsubjects[2-$gmwmN])
      set gmwmtrainingages = ($gmwmtrainingages[2-$gmwmN])
      breaksw
    case "$gmwmN"
      @ endn = $gmwmN - 1
      set gmwmtrainingsubjects = ($gmwmtrainingsubjects[1-$endn])
      set gmwmtrainingages = ($gmwmtrainingages[1-$endn])
      breaksw
    default:
      @ n1 = $counter - 1
      @ n2 = $counter + 1
      set gmwmtrainingsubjects = ($gmwmtrainingsubjects[1-$n1] $gmwmtrainingsubjects[$n2-$gmwmN])
      set gmwmtrainingages     = ($gmwmtrainingages[1-$n1] $gmwmtrainingages[$n2-$gmwmN])
      breaksw
  endsw
  @ gmwmN = $gmwmN - 1
endif
endif
endif

##
## (1) sort the input ages
##
set allages = ($trainingages $subjage)
set sortedages = `echo $allages | fmt -1 | sort -n`
# echo SORTED TRAINING AGES including SUBJ AGE:  $sortedages 

if ($forcegmwm) then
  set allgmwmages = ($gmwmtrainingages $subjage)  
  set sortedgmwmages = `echo $allgmwmages | fmt -1 | sort -n`
# echo SORTED TRAINING AGES including SUBJ AGE:  $sortedgmwmages
endif

##
## (2) find where the input age fits in -- first occurance of the searched age in the sorted age list
##
set subjind = 1
foreach n ($sortedages)
  if ($subjage == $n) then 
    break
  endif
  @ subjind = $subjind + 1
end
# echo SUBJ AGE INDEX $subjind
set ind = $subjind

if ($forcegmwm) then
  set subjind = 1
  foreach n ($sortedgmwmages)
    if ($subjage == $n) then 
      break
    endif
    @ subjind = $subjind + 1
  end
  # echo SUBJ AGE INDEX $subjind
  set gmwmind = $subjind
endif
# echo SUBJ INDEX $subjind

##
## (3) select the desired number of training subjects
##

###
### Find the list of ages
###
set lowerind = 0
set upperind = 0
set lowerdiff = 0
set upperdiff = 0

set selectedages = ()
set selectedagediffs = ()
set counter = 1
@ upperind = $ind + 1
@ lowerind = $ind - 1
# echo COUNTER K $counter $k
@ N = $N + 1 #LZ 02132021
while ($counter <= $k)

  if($upperind > $N) then
    set upperdiff = 1000
  else
    set upper = $sortedages[$upperind]
    @ upperdiff = $upper - $subjage
  endif
  #echo UPPERDIFF $upperdiff

  if($lowerind < 1) then
    set lowerdiff = 1000
  else
  set lower = $sortedages[$lowerind]
  @ lowerdiff =  $subjage - $lower
  endif
  #echo LOWERDIFF $lowerdiff 
 
  # echo DIFFERENCES $lowerdiff  $upperdiff
  if ($lowerdiff < $upperdiff) then
    set selectedages = ($selectedages $sortedages[$lowerind] )
    # echo $selectedages
    set selectedagediffs = ($selectedagediffs $lowerdiff )
    @ lowerind = $lowerind - 1
  else
    set selectedages = ($selectedages $sortedages[$upperind]  )
    # echo $selectedages
    set selectedagediffs = ($selectedagediffs $upperdiff  )
    @ upperind = $upperind + 1
  endif
  #echo SELECTEDAGES $selectedages
  @ counter = $counter + 1
end
#echo THE LIST OF AGES OF THE SELECTED TRAINING SUBJECTS: ${selectedages}

##
## Assemble the list of training subjects
##
set lowerc = 0
set higherc = 0
set selectedsubjects = ()
foreach age ($selectedages)
  set counter = 1
  foreach listage ($trainingages)
    if ($age == $listage) then
      set selectedsubjects = ($selectedsubjects $trainingsubjects[$counter])
      # echo SELECTED SUBJECTS
      @ lowerc  = $counter - 1
      @ higherc = $counter + 1
      # in order to avoid duplications in the selection....
      if ($lowerc == 0) then
        set L = $#trainingages
        set trainingages = ($trainingages[$higherc-$L])
        # echo $trainingages
        set L = $#trainingsubjects
        set trainingsubjects = ($trainingsubjects[$higherc-$L])
        # echo $trainingsubjects
      else 
        if ($higherc > $#trainingages) then
          set trainingages = ($trainingages[1-$lowerc])
          # echo $trainingages
          set trainingsubjects = ($trainingsubjects[1-$lowerc])
          # echo $trainingsubjects
        else
          set L = $#trainingages
          set trainingages = ($trainingages[1-$lowerc] $trainingages[$higherc-$L])
          # echo $trainingages  
          set L = $#trainingsubjects
          set trainingsubjects = ($trainingsubjects[1-$lowerc] $trainingsubjects[$higherc-$L])
          # echo $trainingsubjects
        endif
      endif
      break
    endif
    @ counter = $counter + 1
  end
end

if ($forcegmwm) then
###
### Find the list of gmwm ages
###
set selectedgmwmages = ()
set selectedgmwmagediffs = ()
set counter = 1
@ upperind = $gmwmind + 1
@ lowerind = $gmwmind - 1
# echo COUNTER K $counter $k
while ($counter <= 1) # for now only 1 GMWM subject is required if forcegmwm is set

  if($upperind > $gmwmN) then
    set upperdiff = 1000
  else
    set upper = $sortedgmwmages[$upperind]
    @ upperdiff = $upper - $subjage
  endif

  if($lowerind < 1) then
    set lowerdiff = 1000
  else
  set lower = $sortedgmwmages[$lowerind]
  @ lowerdiff =  $subjage - $lower
  endif
 
  # echo DIFFERENCES $lowerdiff  $upperdiff
  if ($lowerdiff < $upperdiff) then
    set selectedgmwmages = ($selectedgmwmages $sortedgmwmages[$lowerind] )
    # echo $selectedages
    set selectedgmwmagediffs = ($selectedgmwmagediffs $lowerdiff )
    @ lowerind = $lowerind - 1
  else
    set selectedgmwmages = ($selectedgmwmages $sortedgmwmages[$upperind]  )
    # echo $selectedages
    set selectedgmwmagediffs = ($selectedgmwmagediffs $upperdiff  )
    @ upperind = $upperind + 1
  endif
  @ counter = $counter + 1
end
# echo THE LIST OF AGES OF THE SELECTED GMWM TRAINING SUBJECTS: ${selectedgmwmages}

##
## Assemble the list of training subjects
##
set selectedgmwmsubjects = ()
foreach age ($selectedgmwmages)
  set counter = 1
  foreach listage ($gmwmtrainingages)
    if ($age == $listage) then
      set selectedgmwmsubjects = ($selectedgmwmsubjects $gmwmtrainingsubjects[$counter])
      # echo SELECTED SUBJECTS
      @ lowerc  = $counter - 1
      @ higherc = $counter + 1
      # in order to avoid duplications in the selection....
      if ($lowerc == 0) then
        set L = $#gmwmtrainingages
        set gmwmtrainingages = ($gmwmtrainingages[$higherc-$L])
        # echo $gmwmtrainingages
        set L = $#gmwmtrainingsubjects
        set gmwmtrainingsubjects = ($gmwmtrainingsubjects[$higherc-$L])
        # echo $gmwmtrainingsubjects
      else 
        if ($higherc > $#gmwmtrainingages) then
          set gmwmtrainingages = ($gmwmtrainingages[1-$lowerc])
          # echo $gmwmtrainingages
          set gmwmtrainingsubjects = ($gmwmtrainingsubjects[1-$lowerc])
          # echo $gmwmtrainingsubjects
        else
          set L = $#gmwmtrainingages
          set gmwmtrainingages = ($gmwmtrainingages[1-$lowerc] $gmwmtrainingages[$higherc-$L])
          # echo $gmwmtrainingages  
          set L = $#gmwmtrainingsubjects
          set gmwmtrainingsubjects = ($gmwmtrainingsubjects[1-$lowerc] $gmwmtrainingsubjects[$higherc-$L])
          # echo $gmwmtrainingsubjects
        endif
      endif
      break
    endif
    @ counter = $counter + 1
  end
end
endif

########################################################
if ($forcegmwm) then
  set found = 0
  set counter = 1
  foreach s ($selectedsubjects)
    if ($s == $selectedgmwmsubjects[1]) then
      set found = 1
      break
    endif
    @ counter = $counter + 1
  end
  if ($found == 0) then
    @ endn = $k - 1  
    set selectedsubjects = ($selectedsubjects[1-$endn] $selectedgmwmsubjects[1])
  endif
endif

########################################################

# echo The list of training subjects ${selectedsubjects}
echo ${selectedsubjects}
# echo The list of training ages ${selectedages}
if ($writeout) then 
  echo ${selectedsubjects} > $outfile
  # echo ${selectedagediffs} >> $outfile
endif
