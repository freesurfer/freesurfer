#!/usr/bin/env python3

###########################
#
# Takes inputs required for `kvlBuildAtlasMesh` and a training schedule
# and makes repeated calls to `kvlBuildAtlasMesh` according to the training
# schedule, using the output of one call to `kvlBuildAtlasMesh` as input to the
# next (via `explicitStartCollection`)
#
# The file `training-schedule-example.txt` contains an example training schedule
# every line should contain 3 numbers
#   - The number of interations to run for this epoch (unsigned int)
#   - The 'stiffness' factor to pass to kvlBuildAtlasMesh (float)
#   - The 'edgeCollapseEncouragementFactor' to pass to kvlBuildAtlasMesh (float)
#
# Example:
#  /autofs/cluster/gerenuk/pwighton/fs/freesurfer/samseg/gems_train_mesh \
#    --num-upsamples 1 \
#    --mesh-size 3 3 3 \
#    --schedule-file /autofs/cluster/gerenuk/pwighton/fs/freesurfer/samseg/training-schedule-example.txt \
#    --work-dir /autofs/cluster/gerenuk/pwighton/samseg/test-gems-train-mesh-02 \
#    --binary /autofs/cluster/gerenuk/pwighton/samseg/install/gems/bin/kvlBuildAtlasMesh \
#    --label-files \
#        seg_1.mgz \
#        seg_2.mgz \
#        seg_3.mgz
###########################

import os
import sys
import argparse
import tempfile
import errno
import shutil

def parse_args(args):
  parser = argparse.ArgumentParser()
  parser.add_argument('-n','--num-upsamples', required=True, \
                      help='The number of upsapmling steps for `kvlBuildAtlasMesh` to perform.')
  parser.add_argument('-m','--mesh-size', nargs=3, required=True, \
                      help='The mesh size (x, y, z) to pass to `kvlBuildAtlasMesh`.')
  parser.add_argument('-s','--schedule-file', required=False, default=None, \
                      help='Filename containing the training schedule for sucessive calls to `kvlBuildAtlasMesh`')
  parser.add_argument('-w', '--work-dir', required=False, default=None, \
                      help='Directory under which to keep output and intermediary files (will be created)')
  parser.add_argument('-b', '--binary', required=False, default=None, \
                      help='Location of `kvlBuildAtlasMesh` binary')
  parser.add_argument('-l', '--label-files', nargs='+', required=True, \
                      help='The ground truth lables to train against')
  #parser.add_argument('-t', '--lookup-table', required=False, default=None, \
  #                    help='The FreeSurfer-like lookup table')
  args = parser.parse_args()

  if args.schedule_file is None:
    args.schedule = default_training_schedule()
  else:
    args.schedule = read_schedule_file(args.schedule_file)

  # Make sure we can find the kvlBuildAtlasMesh binary    
  if args.binary is None:
    if os.environ.get('FREESURFER_HOME'):
      args.binary = os.path.join(os.environ.get('FREESURFER_HOME'),'gems/bin/kvlBuildAtlasMesh')
  if args.binary is None or not os.path.exists(args.binary):
    print("gems_train_mesh:  ERROR: Can't find kvlBuildAtlasMesh, either set the FREESURFER_HOME env var or use -b")
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.binary)

  # Make sure we can find the FreeSurfer lookup table
  #if args.lookup_table is None:
  #  if os.environ.get('FREESURFER_HOME'):
  #    args.lookup_table = os.path.join(os.environ.get('FREESURFER_HOME'),'FreeSurferColorLUT.txt')
  #if args.lookup_table is None or not os.path.exists(args.lookup_table)
  #  print("gems_train_mesh:  ERROR: Can't find FreeSurfer-like lookup table, either set the FREESURFER_HOME env var or use -t")
  #  raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.binary)
  
  # Make sure all input label files exist
  for file in args.label_files:
    if not os.path.exists(file):
      print("gems_train_mesh:  ERROR: Can't find label file "+file)
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

  # Make the working dir if it doesn't exist
  if args.work_dir is None:
    args.work_dir = tempfile.mkdtemp()
    print("gems_train_mesh:  No workdir specified so using "+args.work_dir)
  else:
    os.makedirs(args.work_dir, exist_ok=True)
  
  return args

# The training shedule format is a list of lists, eg:
# [ [num_itr, stiffness, edgeCollaseFactor],
#   [num_itr, stiffness, edgeCollaseFactor] ]
def default_training_schedule():
  return [ [1, 1.0, 1.0],
           [10, 0.1, 1.25] ]

def read_schedule_file(filename):
  with open(filename) as f:
    schedule = [[float(x) for x in line.split()] for line in f]
  # Every line of the training schedule should have exactly three values:
  #   - num itr (uint)
  #   - stiffness (float)
  #   - edge collapse factor (float)
  for epoch in schedule:
    assert(len(epoch) == 3)  
  return schedule

# Runs kvl_BuildAtlasMesh and returns a filename string with the last meshCollection created   
def run_kvlBuildAtlasMesh(binary, work_dir, epoch_num, num_upsamples, mesh_size, label_files, num_itr=5, stiffness=1.0, edgeCollapse=1.0, startCollection=None):
  launch_dir = os.path.join(work_dir, 'epoch_'+'{:04d}'.format(epoch_num)+'_launch_dir')
  out_dir = os.path.join(work_dir, 'epoch_'+'{:04d}'.format(epoch_num)+'_out_dir')

  os.makedirs(launch_dir, exist_ok=True)
  os.makedirs(out_dir, exist_ok=True)
  
  print(f"gems_train_mesh: epoch: {epoch_num}\n  launch dir: {launch_dir}\n  out dir: {out_dir}")

  # Copy explicitStartCollection.gz to launch dir if specified
  if startCollection is not None:
    startCollection_dest = os.path.join(launch_dir, 'explicitStartCollection.gz')
    shutil.copy(startCollection, startCollection_dest)

  # launch  
  orig_dir = os.getcwd()
  os.chdir(launch_dir)
  launch_cmd = f"{binary} {num_upsamples} {mesh_size[0]} {mesh_size[1]} {mesh_size[2]} {stiffness} {num_itr} {edgeCollapse} {out_dir} "+' '.join(label_files)
  print("  launch command: "+launch_cmd)
  os.system(launch_cmd)
  os.chdir(orig_dir)

  # kvlBuildAtlasMesh will create meshCollections in out_dir up to and including num_itr.
  # If we can't find the last meshCollection, then something went wrong with kvlBuildAtlasMesh 
  lastMeshCollection_file = os.path.join(out_dir, f"CurrentMeshCollection{num_itr}.gz")
  if not os.path.exists(lastMeshCollection_file):
    print(f"gems_train_mesh:  ERROR: Can't find the file {lastMeshCollection_file} that we were expecting kvlBuildAtlasMesh to create")
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), lastMeshCollection_file)
  return lastMeshCollection_file

def main(argv):
  args = parse_args(argv)
  orig_cwd = os.getcwd()
  startCollection = None
  try:
    for epoch_num, epoch in enumerate(args.schedule):
      kvl_num_itr = int(epoch[0])
      kvl_stiffness = float(epoch[1])
      kvl_edge_collapse = float(epoch[2])
      startCollection = run_kvlBuildAtlasMesh(\
                              binary=args.binary, \
                              work_dir=args.work_dir, \
                              epoch_num=epoch_num, \
                              num_upsamples=args.num_upsamples, \
                              mesh_size=args.mesh_size, \
                              label_files=args.label_files, \
                              num_itr=kvl_num_itr, \
                              stiffness=kvl_stiffness, \
                              edgeCollapse=kvl_edge_collapse, \
                              startCollection=startCollection)
  except Exception as e:
    os.chdir(orig_cwd)
    raise e

if __name__ == "__main__":
  sys.exit(main(sys.argv))
