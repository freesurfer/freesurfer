#! /bin/tcsh -f

##############
# process input paramters
##############

# echo ARG LIST $argv
## read input
set kneighbor = 0;
set trainingsubjects = ()
set cmdline = ($argv);
while( $#argv != 0 )
  set flag = $argv[1]; shift;
  switch($flag)
    case "--subj":
      if ( $#argv < 1) goto arg1err;
      set testsubj = $argv[1]; shift;
      breaksw
    case "--ts":
      if ( $#argv < 1) goto arg1err;
      set counter = 1
      echo ARG $argv
      set N = $#argv
      echo ARG COUNT $N
      while ($counter <= $N)
        set trainingsubjects = ($trainingsubjects $argv[$counter]);
        @ counter = $counter + 1
      end
      # shift; shift
      set counter = 1
      while ($counter <= $N)
        shift;
        @ counter = $counter + 1
      end
      set selected = 1;
      echo TRAINING SUBJECTS $trainingsubjects
      breaksw
    case "--kn":
      if ( $#argv < 1) goto arg1err;
      set kneighbor = $argv[1]; shift;
      set selected = 0;
      breaksw
    default:
      echo ERROR: Flag $flag unrecognized.
      echo $cmdline
      exit 1
      breaksw
  endsw
end

echo CHECKING THE INPUTS...
## check input
if($#testsubj == 0) then
  echo "ERROR: must spec a subject to be segmented"
  exit 1;
endif
if(($kneighbor == 0) && ($#trainingsubjects == 0)) then
  echo "ERROR: must identify the number of neighbors or the trainingsubjects to be used for segmentation"
  exit 1;
endif
if(($kneighbor > 0) && ($#trainingsubjects > 0)) then
  echo "ERROR: must identify trainingsubjects in a unique way (either kn or ts)!"
  exit 1;
endif

##############
# assemble inputs: training data, seg options
##############
if ($kneighbor == 0) then
  set kneighbor = $#trainingsubjects
endif

echo NUMBER OF NEIGHBORS $kneighbor

set testnorm  = $MASKEDINPUT
set outputdir = $SEGOUTDIR

set outfile = $outputdir/$testsubj.$kneighbor.neighborlist.txt
if !($selected) then
  set cmd = ($FSSCRIPTSDIR/find_neighbors_byAge.csh --s $testsubj --k $kneighbor --outfile $outfile)
  echo $cmd; 
  $cmd
  set trainingsubjects = (`cat $outfile`)
else 
  echo $trainingsubjects > $outfile
endif
echo TRAINING SUBJECTS $trainingsubjects

set trainingdata = ()  
foreach subj ($trainingsubjects)
  set DRAMMSmovedmanseg = $outputdir/manseg-${subj}-2-${testsubj}.DRAMMS.nii.gz
  echo MOVED SEG FILE $DRAMMSmovedmanseg
  if (-e $DRAMMSmovedmanseg) then
    set trainingdata = ($trainingdata $DRAMMSmovedmanseg)
  endif
end
echo TRAININGDATA $trainingdata
set N = $#trainingdata 

set rho = `echo "0.5 + $N * 0.05" | bc`
#set rho = `echo "0.2 + $N * 0.05" | bc`
set options = ($MAXLAB $BIASFIELDORDER $rho $BETA)
echo OPTIONS $options

##############
# assemble hack options
##############

## decide whether hacks are needed for WM labels: (a) cerebral, (b) cerebellar and (c) basal ganglia
set counter = 1
set val = 0
# while ( $val == 0 && $counter <= $N) 
while ( $counter <= $N) 
  set d = $trainingdata[$counter]
  echo DATA $d
  #source $FSSCRIPTSDIR/find_cerebral_wm_labels.csh $d
  ## set val = $status
  #@ val = $val + $status
  @ val = $val + `$FSSCRIPTSDIR/find_cerebral_wm_labels.csh $d` 
  echo VAL $val
  @ counter = $counter + 1
end
#if (($val > 0) && ($N > 1)) then
if (($val >= 1) && ($N > $val) && ($N > 1)) then
  set lfoptions1 = (1 $LFOPTIONS1)
else
  set lfoptions1 = (0 $NEGLFOPTIONS1)
endif
## 
set counter = 1
set val = 0
# while ( $val == 0 && $counter <= $N) 
while ( $counter <= $N) 
  set d = $trainingdata[$counter]
  source $FSSCRIPTSDIR/find_cerebellar_wm_labels.csh $d
  # set val = $status
  @ val = $val + $status
  @ counter = $counter + 1
end
#if ($val > 0) then
if (($val >= 1) && ($N > $val) && ($N > 1)) then
  set lfoptions2 = (1 $LFOPTIONS2)
else
  set lfoptions2 = (0 $NEGLFOPTIONS2)
endif
##

set counter = 1
set val = 0
# while ( $val == 0 && $counter <= $N) 
while ( $counter <= $N)
  set d = $trainingdata[$counter]
  source $FSSCRIPTSDIR/find_putamen_labels.csh $d
  # set val = $status
  @ val = $val + $status
  @ counter = $counter + 1
end
# if ($val > 0) then
if (($val >= 1) && ($N > $val) && ($N > 1)) then
  set lfoptions3 = (1 $LFOPTIONS3)
else
  set lfoptions3 = (0 $NEGLFOPTIONS3)
endif

##############
# run process
##############
set postfix = ()
#foreach s ($trainingsubjects)
#  set postfix = ($postfix$s)
#end
echo POSTFIX $postfix 
#set outputprefix = $outputdir/${testsubj}-with$postfix-BETA$BETA-RHO$rho
set outputprefix = $outputdir/${testsubj}-BETA$BETA-RHO$rho
echo OUTPUT PREFIX $outputprefix
if !(-e ${outputprefix}.labels.mrf.nii.gz) then 
  set cmd = ($CMDDIR/run_LabelFusion.sh $MLAB $lfoptions1 $lfoptions2 $lfoptions3 $testnorm $trainingdata $options $outputprefix $MRFparams)
  echo "\n"
  echo LABELFUSION CMD $cmd; 
  echo "\n"
  $cmd;
endif
