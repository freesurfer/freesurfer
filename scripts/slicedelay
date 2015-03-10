#!/usr/bin/python
import sys;
import os;
import string;

#---------------------------------------------------------
def print_usage():
  print "USAGE: slicedelay --help";
  print "  --o slicedelayfile";
  print "  --nslices nslices : total number of slices in the volume";
  print "  --order order (up,down,odd,even,siemens)";
  print "  --ngroups ngroups (number of slice groups for SMS)";
  return 0;
#end def print_usage:

#---------------------------------------------------------
def print_help():
  print_usage();
  print "\nCreates an FSL custom slice delay file for use with slicetimer (--tcustom=sdfile).\n"\
      "It has a single column of values, one for each slice. Each value is the slice delay\n" \
      "measured as a fraction of the TR and range from +0.5 (beginning of the TR) to \n"\
      "-0.5 (end of the TR). Used for slice-time correction of fMRI.\n"
#end def print_help:

#---------------------------------------------------------
def argnerr(narg,flag):
  print "ERROR: flag %s requires %d arguments" % (flag,narg);
  sys.exit(1);
#end def parse_args(argv)

#---------------------------------------------------------
def parse_args(argv):
  global sdf;
  global nslices;
  global order;
  global ngroups;
  global debug;

  del argv[0]; # get past program name (like shift)

  while(len(argv) != 0):
    flag = argv[0];
    del argv[0];
    if(debug): print "flag = %s" % flag;

    if(flag == "--o"):
      if(len(argv) < 1): argnerr(1,flag);
      sdf = argv[0]; del argv[0];
    elif(flag == "--nslices"):
      if(len(argv) < 1): argnerr(1,flag);
      nslices = int(argv[0]); del argv[0];
    elif(flag == "--order"):
      if(len(argv) < 1): argnerr(1,flag);
      order = argv[0]; del argv[0];
    elif(flag == "--ngroups"):
      if(len(argv) < 1): argnerr(1,flag);
      ngroups = int(argv[0]); del argv[0];
    elif(flag == "--up"):      order = "up";
    elif(flag == "--down"):    order = "down";
    elif(flag == "--odd"):     order = "odd";
    elif(flag == "--even"):    order = "even";
    elif(flag == "--siemens"): order = "siemens";
    elif(flag == "--help"):
      print_help();
      sys.exit(1);
    elif(flag == "--debug"):
      debug = 1;
    else:
      print "ERROR: flag %s not recognized" % flag; 
      sys.exit(1);
    #endif
  #endwhile
  return 0;

#end def parse_args(argv)

#---------------------------------------------------------
def isflag(arg):
  if(len(arg) < 3): return 0;
  if(arg[0] == "-" and arg[1] == "-"): return 1;
  return 0;
# end def isflag(arg)

#---------------------------------------------------------
def check_args():
  global sdf;
  global nslices;
  global order;
  global ngroups;
  global debug;

  if(len(sdf) == 0):
    print "ERROR: output file needed";
    sys.exit(1);
  #endif    
  if(nslices == 0):
    print "ERROR: nslices needed";
    sys.exit(1);
  #endif    
  if(len(order) == 0):
    print "ERROR: order needed";
    sys.exit(1);
  #endif  

  return 0;
#end check_args()

#-----------------------------------------------------------
# ------ main -----------------------------------------------
#-----------------------------------------------------------

sdf = ();
nslices = 0;
order = ();
ngroups = 1;
debug = 0;

nargs = len(sys.argv) - 1;
if(nargs == 0):
  print_usage();
  sys.exit(0);
#end
parse_args(sys.argv);
check_args();

err = nslices%ngroups;
if(err):
  print "ERROR: cannot divide nslices=%d by ngroups=%d " % (nslices,ngroups);
  sys.exit(0);
#endif
nslicespg = nslices/ngroups;

print "sdf %s" % sdf;
print "nslices %d" % nslices;
print "order %s" % order;
print "ngroups %d" % ngroups;
print "nslicespg %d" % nslicespg;

AcqSliceDelay = [];
for sno in range(1,nslicespg+1):
  D = ((nslicespg-1.0)/2.0-(sno-1.0))/nslicespg;
  AcqSliceDelay.append(D);
#end
print  AcqSliceDelay

AnatSliceDelay = [];
AnatSliceOrder = ();
for sg in range(1,ngroups+1):
  print "sg %d" % (sg);
  if(order == "up"):   AnatSliceOrder = range(1,nslicespg+1,+1);
  if(order == "down"): AnatSliceOrder = range(nslicespg,0,-1);
  if(order == "odd"):  AnatSliceOrder = range(1,nslicespg+1,+2) + range(2,nslicespg+1,+2);
  if(order == "even"): AnatSliceOrder = range(2,nslicespg+1,+2) + range(1,nslicespg+1,+2);
  if(order == "siemens"):  
    if(nslicespg%2==1): # siemens-odd
      AnatSliceOrder = range(1,nslicespg+1,+2) + range(2,nslicespg+1,+2);
    else: # siemens-even
      AnatSliceOrder = range(2,nslicespg+1,+2) + range(1,nslicespg+1,+2);
    #endif
  #endif
  if(len(AnatSliceOrder) == 0):
    print "ERROR: slice order %s not recognized " % (order);
    sys.exit(0);
  #endif

  AnatAcqSliceOrder0 =  sorted(range(len(AnatSliceOrder)),key=lambda x:AnatSliceOrder[x]);
  for s in AnatAcqSliceOrder0:
    AnatSliceDelay.append(AcqSliceDelay[s]);
  #endfor

#end for sg

fp = open(sdf,'w');
for d in AnatSliceDelay:
  fp.write('%15.13f\n' % d);
#end
err = fp.close();

sys.exit(0);
#-------------------------------------------------#






