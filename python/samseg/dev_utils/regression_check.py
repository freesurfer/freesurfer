from samseg.dev_utils.compute_dice_measure import find_image_and_path, find_registration_image_and_path
from samseg.dev_utils.debug_client import compare_ndarray_dice


def compare_results(prior_result_path, current_result_path, case_name, finder, test_type='final'):
    prior_image_buffer, prior_path = finder(prior_result_path, case_name)
    if prior_image_buffer is None:
        print('no prior file at {0}'.format(prior_path))
        return
    current_image_buffer, current_path = finder(current_result_path, case_name)
    if current_image_buffer is None:
        print('no current file at {0}'.format(current_path))
        return
    compare_ndarray_dice(current_image_buffer, prior_image_buffer,
                         '/'.join((case_name, test_type, 'current vs prior')))

def compare_final_image_results(prior_result_path, current_result_path, case_name):
    return compare_results(
        prior_result_path,
        current_result_path,
        case_name,
        finder=find_image_and_path,
        test_type='final'
    )

def compare_registration_image_results(prior_result_path, current_result_path, case_name):
    return compare_results(
        prior_result_path,
        current_result_path,
        case_name,
        finder=find_registration_image_and_path,
        test_type='register'
    )

if __name__ == '__main__':
    prior_result_path = '/home/willy/work/cm/python_temp_data'
    current_result_path = '/home/willy/work/cm/innolitics_testing/python_temp_data'
    case_name = '004'
    compare_final_image_results(prior_result_path, current_result_path, case_name)
    compare_registration_image_results(prior_result_path, current_result_path, case_name)
