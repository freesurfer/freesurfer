#! /bin/tcsh -f

##############
# process input paramters
##############

# echo ARG LIST $argv
## read input
set kneighbor = 0;
set regtype   = NIFTYREG;
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
    case "--regtype":
      if ( $#argv < 1) goto arg1err;
      set regtype = $argv[1]; shift;
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
  set movedmanseg = $outputdir/manseg-${subj}-2-${testsubj}.$regtype.nii.gz
  echo MOVED SEG FILE $movedmanseg
  if (-e $movedmanseg) then
    set trainingdata = ($trainingdata $movedmanseg)
  endif
end
echo TRAININGDATA $trainingdata
set N = $#trainingdata 

# set rho = `echo "0.2 + $N * 0.05" | bc`
set rho = `echo "0.5 + $N * 0.05" | bc`
set beta = 0.3
set maxlab = 3
set biasfieldorder = 4

##############
# run process
##############
set outputfile = $outputdir/${testsubj}-BETA${beta}-RHO${rho}.labels.mrf.nii.gz
echo OUTPUT FILE ${outputfile}
if !(-e ${outputfile}) then
  set cmd = (mri_label_fusion -i $testnorm -s $trainingdata -o $outputfile)
  # options
  set cmd = ($cmd --smooth --rho $rho --beta $beta --bf-order $biasfieldorder --max-lab $maxlab --unary-weight 5 --verbose)
  # cerebral WM hack
  set cmd = ($cmd -e 2 41 3 42)
  # cerebellar WM hack
  set cmd = ($cmd -e 7 46 8 47)
  # basal ganglia hack
  set cmd = ($cmd -e 12 51 13 52)
  echo "\n"
  echo LABELFUSION CMD $cmd; 
  echo "\n"
  $cmd;
endif
