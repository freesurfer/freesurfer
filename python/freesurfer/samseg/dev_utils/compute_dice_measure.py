import os
import freesurfer.gems as gems

from samseg.dev_utils.debug_client import run_test_cases, compare_ndarray_dice, create_checkpoint_manager

RESULT_LEAF_NAME = 'crispSegmentation.nii'
RESULT_REGISTRATION_LEAF_NAME = 'template_coregistered.nii'


def dice_compare_test_case(case_name, case_file_folder, save_path, cross_compare):
    print("")
    gold_image_buffer, gold_file_path = find_gold_image(case_name)
    if gold_image_buffer is None:
        print('no gold file at {0}'.format(gold_file_path))
        return
    python_image_buffer, python_file_path = find_python_image(save_path)
    if python_image_buffer is None:
        print('no python file at {0}'.format(python_file_path))
        return

    compare_ndarray_dice(gold_image_buffer, python_image_buffer,
                         '/'.join((case_name, 'gold vs python')))
    if cross_compare:
        checkpoint_manager = create_checkpoint_manager(case_file_folder)
        python_checkpoint_image = checkpoint_manager.load_python('part3', 1).get('uncroppedFreeSurferSegmentation')
        if python_checkpoint_image is not None:
            compare_ndarray_dice(gold_image_buffer, python_checkpoint_image,
                                 '/'.join((case_name, 'gold vs python checkpoint')))
            compare_ndarray_dice(python_image_buffer, python_checkpoint_image,
                                 '/'.join((case_name, 'python vs python checkpoint')))

            matlab_image_buffer, matlab_file_path = find_matlab_image(save_path, case_name)
            if matlab_image_buffer is not None:
                compare_ndarray_dice(gold_image_buffer, matlab_image_buffer,
                                     '/'.join((case_name, 'gold vs matlab')))
                compare_ndarray_dice(matlab_image_buffer, python_image_buffer,
                                     '/'.join((case_name, 'matlab vs python')))
            else:
                print('no data at {0}'.format(matlab_file_path))

def find_gold_image(case_name):
    return find_image_and_path(os.getenv('GOLD_REFERENCE_DIR'), case_name)

def find_python_image(save_path):
    return find_image_and_path(save_path)

def find_matlab_image(save_path, case_name):
    matlab_save_path = os.path.join(os.path.dirname(os.path.dirname(save_path)), 'matlab_temp_data')
    return find_image_and_path(matlab_save_path, case_name)

def find_image_and_path(case_dir, case_name=None, leaf_name=RESULT_LEAF_NAME):
    if case_name:
        case_dir = os.path.join(case_dir, case_name)
    file_path = os.path.join(case_dir, leaf_name)
    if os.path.isfile(file_path):
        image_buffer = gems.KvlImage(file_path).getImageBuffer()
        return image_buffer, file_path
    else:
        return None, file_path

def find_registration_image_and_path(case_dir, case_name=None, leaf_name=RESULT_REGISTRATION_LEAF_NAME):
    return find_image_and_path(case_dir, case_name=case_name, leaf_name=leaf_name)


def simple_dice_compare(case_name, case_file_folder, savePath):
    dice_compare_test_case(case_name, case_file_folder, savePath, cross_compare=False)


if __name__ == '__main__':
    run_test_cases(action=simple_dice_compare)
