import os
import sys

from as_python.samseg.dev_utils.debug_client import create_checkpoint_manager, run_test_cases, \
    create_part3_inspection_team
from as_python.samseg.run_samseg_ported import run_samseg


def run_single_full_case(case_file_folder, savePath):
    image_file_path = os.path.join(case_file_folder, 'orig.mgz')
    if os.path.isfile(image_file_path):
        print('testing with file {0}'.format(image_file_path))
        checkpoint_manager = create_checkpoint_manager(case_file_folder)
        run_samseg(
            imageFileNames=[image_file_path],
            savePath=savePath,
            checkpoint_manager=checkpoint_manager
        )
        create_part3_inspection_team().inspect_all(checkpoint_manager)
    else:
        print("{0} is not a file".format(image_file_path))


if __name__ == '__main__':
    run_test_cases(action=run_single_full_case)
