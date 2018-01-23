import json

import os

from as_python.samseg.compute_dice_measure import find_gold_image, find_python_image, \
    find_matlab_image
from as_python.samseg.dev_utils.debug_client import valid_case_folders_and_save_paths, \
    measure_label_differences, find_testing_dir


def measure_and_report(report_save_path=None):
    if report_save_path is None:
        report_save_path = os.path.join(find_testing_dir(), 'test_data.json')
    save_as_json(report(measure()), report_save_path)

def save_as_json(measurements, report_save_path):
    with open(report_save_path, mode='w') as outfile:
        json.dump(measurements, outfile, indent=4)


def measure():
    return {case_name: case_measurements(case_name, case_file_folder, save_path) for
            case_name, case_file_folder, save_path in valid_case_folders_and_save_paths()}


def case_measurements(case_name, case_file_folder, save_path):
    measures = {}
    named_images = find_named_images(case_name, case_file_folder, save_path)
    for i, first_named_image in enumerate(named_images):
        for j, second_named_image in enumerate(named_images):
            if j > i:
                first_name, first_image = first_named_image
                second_name, second_image = second_named_image
                key = " versus ".join((first_name, second_name))
                measures[key] = measure_label_differences(first_image, second_image)
    return measures


def find_named_images(case_name, case_file_folder, save_path):
    possible_named_images = []
    gold_image_buffer, _ = find_gold_image(case_name)
    possible_named_images.append(('gold', gold_image_buffer))
    python_image_buffer, _ = find_python_image(save_path)
    possible_named_images.append(('python', python_image_buffer))
    matlab_image_buffer, _ = find_matlab_image(save_path, case_name)
    possible_named_images.append(('matlab', matlab_image_buffer))
    return [(name, buffer) for name, buffer in possible_named_images if buffer is not None]


def report(measurements):
    sorted_case_names = sorted(measurements.keys())
    for case_name in sorted_case_names:
        print(case_name)
        report_case_measurements(measurements[case_name])
        print()
    return measurements


def report_case_measurements(case_measurements):
    sorted_comparisons = sorted(case_measurements.keys())
    for comparison in sorted_comparisons:
        print('  {}'.format(comparison))
        report_comparison_measurments(case_measurements[comparison])


def report_comparison_measurments(comparison_measurments):
    def show(message):
        print('    {}'.format(message))

    show(match_message(comparison_measurments))
    show(dice_message(comparison_measurments))
    show(jaccard_message(comparison_measurments))


def match_message(comparison_measurments):
    return 'Match   {0} = {1}/{2}'.format(
        comparison_measurments['Match'],
        comparison_measurments['total_matches'],
        comparison_measurments['total_size'],
    )


def dice_message(comparison_measurments):
    return 'Dice    {0} = 2*{1}/({2}+{3})'.format(
        comparison_measurments['Dice'],
        comparison_measurments['interior_match_count'],
        comparison_measurments['expected_interior_count'],
        comparison_measurments['actual_interior_count'],
    )


def jaccard_message(comparison_measurments):
    return 'Jaccard {0} = {1}/{2}'.format(
        comparison_measurments['Jaccard'],
        comparison_measurments['interior_match_count'],
        comparison_measurments['interior_union_count'],
    )


if __name__ == '__main__':
    measure_and_report()
