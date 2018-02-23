from samseg.dev_utils.compute_dice_measure import find_image_and_path
from samseg.dev_utils.debug_client import compare_ndarray_dice


def compare_results(prior_result_path, current_result_path, case_name):
    prior_image_buffer, prior_path = find_image_and_path(prior_result_path, case_name)
    if prior_image_buffer is None:
        print('no prior file at {0}'.format(prior_path))
        return
    current_image_buffer, current_path = find_image_and_path(current_result_path, case_name)
    if current_image_buffer is None:
        print('no current file at {0}'.format(current_path))
        return
    compare_ndarray_dice(current_image_buffer, prior_image_buffer,
                         '/'.join((case_name, 'current vs prior')))


if __name__ == '__main__':
    prior_result_path = '/home/willy/work/cm_p/python_temp_data'
    current_result_path = '/home/willy/work/cm_p/innolitics_testing/buckner40/python_temp_data'
    case_name = '004'
    compare_results(prior_result_path, current_result_path, case_name)
