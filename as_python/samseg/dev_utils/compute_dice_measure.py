import os

from gems2python import GEMS2Python


from samseg.dev_utils.debug_client import run_test_cases, compare_ndarray_dice, create_checkpoint_manager

RESULT_LEAF_NAME = 'crispSegmentation.nii'


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
    gold_dir = os.getenv('GOLD_REFERENCE_DIR')
    gold_case_dir = os.path.join(gold_dir, case_name)
    gold_file_path = os.path.join(gold_case_dir, RESULT_LEAF_NAME)
    if os.path.isfile(gold_file_path):
        gold_image_buffer = GEMS2Python.KvlImage(gold_file_path).getImageBuffer()
        return gold_image_buffer, gold_file_path
    else:
        return None, gold_file_path

def find_python_image(save_path):
    python_file_path = os.path.join(save_path, RESULT_LEAF_NAME)
    if os.path.isfile(python_file_path):
        python_image_buffer = GEMS2Python.KvlImage(python_file_path).getImageBuffer()
        return python_image_buffer, python_file_path
    else:
        return None, python_file_path

def find_matlab_image(save_path, case_name):
    matlab_save_path = os.path.join(os.path.dirname(os.path.dirname(save_path)), 'matlab_temp_data')
    matlab_case_dir = os.path.join(matlab_save_path, case_name)
    matlab_file_path = os.path.join(matlab_case_dir, RESULT_LEAF_NAME)
    if os.path.isfile(matlab_file_path):
        matlab_image_buffer = GEMS2Python.KvlImage(matlab_file_path).getImageBuffer()
        return matlab_image_buffer, matlab_file_path
    else:
        return None, matlab_file_path


def simple_dice_compare(case_name, case_file_folder, savePath):
    dice_compare_test_case(case_name, case_file_folder, savePath, cross_compare=False)


if __name__ == '__main__':
    run_test_cases(action=simple_dice_compare)
