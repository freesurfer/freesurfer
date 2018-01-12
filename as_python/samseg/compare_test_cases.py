import sys
import GEMS2Python

from as_python.samseg.dev_utils.debug_client import measure_closeness, measure_scalar, measure_list
from as_python.samseg.run_samseg_test_case import run_samseg_test_cases, create_checkpoint_manager


# 'biasFieldCoefficients',
# 'colors',
# 'croppingOffset',
# 'FreeSurferLabels',
# 'imageBuffers',
# 'imageSize',
# 'imageToWorldTransformMatrix',
# 'kroneckerProductBasisFunctions',
# 'mask',
# 'names',
# 'nonCroppedImageSize',
# 'numberOfBasisFunctions',
# 'numberOfClasses',
# 'numberOfContrasts',
# 'numberOfGaussians',
# 'numberOfGaussiansPerClass',
# 'reducingLookupTable',
# 'savePath',
# 'transformMatrix',
# 'voxelSpacing',


def compare_single_case(case_file_folder, savePath):
    checkpoint_manager = create_checkpoint_manager(case_file_folder)
    # matlab_checkpoint_0 = checkpoint_manager.load('part0', 1)
    # register_atlas_checkpoint = checkpoint_manager.load('register_atlas_fixture', 1)
    things_to_compare = [
        ['part1', 1, {
            'measure_closeness': [
                'biasFieldCoefficients',
                'colors',
                'croppingOffset',
                'FreeSurferLabels',
                'imageBuffers',
                'imageSize',
                'imageToWorldTransformMatrix',
                # 'kroneckerProductBasisFunctions',
                'mask',
                # 'names',
                'nonCroppedImageSize',
                'reducingLookupTable',
                # 'savePath',
                'transformMatrix',
                'voxelSpacing',
            ],
            'measure_scalar': [
                'numberOfClasses',
                'numberOfContrasts',
                'numberOfGaussians',
            ],
            'measure_list': [
                # 'numberOfBasisFunctions',
                'numberOfGaussiansPerClass',
            ]
        }
         ],
        ['part2', 1, {
            'measure_closeness': [
            ],
        }
         ]
    ]
    for checkpoint, index, action_list in things_to_compare:
        try:
            matlab_checkpoint = checkpoint_manager.load(checkpoint, index)
        except Exception as flaw:
            print('no matlab checkpoint {0} for {1}'.format(checkpoint, case_file_folder))
            continue
        try:
            python_checkpoint = checkpoint_manager.load_python(checkpoint, index)
        except Exception as flaw:
            print('no python checkpoint {0} for {1}'.format(checkpoint, case_file_folder))
            continue
        for key in sorted(python_checkpoint.keys()):
            if key in matlab_checkpoint:
                print('found matlab {0}'.format(key))
                test_message =':'.join((checkpoint, str(index), key))
                python_value = python_checkpoint[key]
                matlab_value = matlab_checkpoint[key]
                if key in action_list['measure_closeness']:
                    measure_closeness(matlab_value, python_value, test_message)
                if key in action_list['measure_scalar']:
                    measure_scalar(matlab_value, python_value, test_message)
                if key in action_list['measure_list']:
                    measure_list(matlab_value, python_value, test_message)
            else:
                print('could not find {0} in matlab'.format(key))


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
