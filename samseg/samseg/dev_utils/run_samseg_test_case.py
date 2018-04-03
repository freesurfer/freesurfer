import os

from samseg.dev_utils.debug_client import create_checkpoint_manager, run_test_cases
from samseg.run_samseg_ported import run_samseg
from samseg.show_figures import ShowFigures

USE_CHECKPOINT_MANAGER = False
THREAD_COUNT = 6

def run_single_full_case(case_name, case_file_folder, savePath):
    image_file_path = os.path.join(case_file_folder, 'orig.mgz')
    if os.path.isfile(image_file_path):
        print('testing with file {0}'.format(image_file_path))
        checkpoint_manager = create_checkpoint_manager(case_file_folder) if  USE_CHECKPOINT_MANAGER else None
        run_samseg(
            imageFileNames=[image_file_path],
            savePath=savePath,
            checkpoint_manager=checkpoint_manager,
            numberOfThreads=THREAD_COUNT,
            atlas_only=False,
            showFigures=ShowFigures(interactive=False, movie_flag=True, show_flag=False),
        )
    else:
        print("{0} is not a file".format(image_file_path))


if __name__ == '__main__':
    run_test_cases(action=run_single_full_case)