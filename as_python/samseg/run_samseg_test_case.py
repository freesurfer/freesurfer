import os
import sys

from as_python.samseg.dev_utils.debug_client import CheckpointManager
from as_python.samseg.run_samseg_ported import run_samseg

def run_single_case(case_file_folder, savePath):
    image_file_path = os.path.join(case_file_folder, 'orig.mgz')
    if os.path.isfile(image_file_path):
        print('testing with file {0}'.format(image_file_path))
        run_samseg(
            imageFileNames=[image_file_path],
            savePath=savePath,
            checkpoint_manager=create_checkpoint_manager(case_file_folder)
        )
    else:
        print("{0} is not a file".format(image_file_path))

def run_samseg_test_cases(case_names=None, testing_dir=None, action=run_single_case):
    if testing_dir is None:
        testing_dir = os.getenv('TESTING_DIR')
    if not testing_dir or not os.path.isdir(testing_dir):
        raise ValueError("testing_dir={0}: must be a valid directory".format(testing_dir))
    if not case_names:
        case_names = sorted([f for f in os.listdir(testing_dir)])
    for case_name in case_names:
        if not 'temp_data' in case_name:
            case_file_folder = os.path.join(testing_dir, case_name)
            savePath = os.path.join(testing_dir, 'python_temp_data', case_name)
            if os.path.isdir(case_file_folder):
                action(case_file_folder, savePath)
            else:
                print("{0} is not a folder".format(case_file_folder))

def create_checkpoint_manager(case_file_folder):
    matlab_dump_dir = make_checkpoint_dir(case_file_folder, 'matlab_checkpoints')
    python_dump_dir = make_checkpoint_dir(case_file_folder, 'python_checkpoints')
    return CheckpointManager(matlab_dump_dir=matlab_dump_dir, python_dump_dir=python_dump_dir)

def make_checkpoint_dir(case_file_folder, subdir_name):
    checkpoint_dir = os.path.join(case_file_folder, subdir_name)
    os.makedirs(checkpoint_dir, exist_ok=True)
    return checkpoint_dir

if __name__ == '__main__':
    run_samseg_test_cases(sys.argv[1:])
