import sys

import GEMS2Python

from as_python.samseg.dev_utils.debug_client import create_part1_inspection_team, \
    create_part2_inspection_team, create_part3_inspection_team, create_reduced_alphas_inspection_team, \
    create_optimizer_inspection_team, create_multiresWarp_inspection_team, create_optimizer_exit_inspection_team, \
    create_optimizer_em_exit_inspection_team, create_bias_correction_inspection_team
from as_python.samseg.run_samseg_test_case import run_samseg_test_cases, create_checkpoint_manager


def compare_single_case(case_file_folder, savePath):
    checkpoint_manager = create_checkpoint_manager(case_file_folder)
    # matlab_checkpoint_0 = checkpoint_manager.load('part0', 1)
    # register_atlas_checkpoint = checkpoint_manager.load('register_atlas_fixture', 1)
    for team in [
        create_part1_inspection_team(),
        create_bias_correction_inspection_team(),
        create_reduced_alphas_inspection_team(),
        create_optimizer_inspection_team(),
        create_multiresWarp_inspection_team(),
        create_optimizer_exit_inspection_team(),
        create_optimizer_em_exit_inspection_team(),
        create_part2_inspection_team(),
        create_part3_inspection_team(),
    ]:
        team.inspect_all(checkpoint_manager)


def compare_template_files():
    matlab_file_name = '/home/willy/work/cm/innolitics_testing/buckner40/matlab_temp_data/004/template_coregistered.nii'
    python_file_name = '/home/willy/work/cm/innolitics_testing/buckner40/python_temp_data/004/template_coregistered.nii'
    matlab_image = GEMS2Python.KvlImage(matlab_file_name)
    python_image = GEMS2Python.KvlImage(python_file_name)
    matlab_transform_matrix = matlab_image.transform_matrix.as_numpy_array
    python_transform_matrix = python_image.transform_matrix.as_numpy_array
    matlab_image_buffer = matlab_image.getImageBuffer()
    python_image_buffer = python_image.getImageBuffer()
    print('hello')


if __name__ == '__main__':
    # run_samseg_test_cases(['004'], action=compare_single_case)
    # compare_template_files()
    run_samseg_test_cases(sys.argv[1:], action=compare_single_case)
